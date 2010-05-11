#!/usr/local/bin/perl

use strict;
use warnings;
use OBO::Core::Ontology;
use OBO::Core::Term;
use OBO::Core::Def;
use OBO::Core::Relationship;
use OBO::Core::RelationshipType;
use Data::Dumper;

###################################
my $id_prefix = "EC";
my $default_namespace = "EC";
my $xref_db = "EC";
###################################

##############################
my $id_num = "1";
my @transferred = ();
##############################

my $enzclass = $ARGV[0] or die("enzclass");
my $enzyme_dat = $ARGV[1] or die("Please pass in the enzyme.dat file");

my $ontology = &setup_ontology;

my %enzyme_class = &parse_enzclass( $enzclass );
&enzyme_classes_to_ontology( \%enzyme_class, $ontology );

print STDERR "Parsing dat file and add to ontology\n";
&parse_enzyme_dat( $enzyme_dat );

print STDERR "Adding replaced by\n";
&add_replaced_by( \@transferred, $ontology );

$ontology->export(\*STDOUT, 'obo');

sub enzyme_classes_to_ontology {
    my ($ec, $ont) = @_;
    my @relationships;
    
    foreach my $ec_class ( keys %{$ec} ) {
        my ($new_def, $is_a) = &concat_parent_defs( $ec_class, $ec );
        my $term = new OBO::Core::Term();
        $term->id( &next_id );
        $term->name( $new_def );
        $term->xref_set_as_string("[$xref_db:$ec_class]");
        if( $is_a ) {
            push(@relationships, { 'parent' => $is_a, 'child' => $ec_class } );
        }
        $ont->add_term($term);
    }

    ## now add the relationships
    foreach my $rel ( @relationships ) {
        my $parent_term = $ont->get_term_by_xref( $xref_db, $rel->{'parent'} );
        my $child_term  = $ont->get_term_by_xref( $xref_db, $rel->{'child'}  );
        my $relationship = new OBO::Core::Relationship();
        $relationship->type('is_a');
        $relationship->link($child_term, $parent_term);
        $relationship->id( $child_term->id()."_is_a_".$parent_term->id() );
        $ont->add_relationship( $relationship );
    }
    
}

sub concat_parent_defs {
    my ($ec_num, $defs) = @_;
    my $new_def = $defs->{$ec_num};
    my $parent;
    unless( $ec_num =~ /\d+(\.-){3}/ ) {
        $parent = &get_direct_parent_ec( $ec_num );
        my ($parent_def) = &concat_parent_defs($parent, $defs);
        die("Could not get parent_definition for $parent") unless( defined( $parent_def ) );
        $new_def = $parent_def." ".$defs->{$ec_num};
    }
    return $new_def, $parent;
}

sub get_direct_parent_ec {
    my ($ec_num) = @_;
    my @nums = split(/\./, $ec_num );
    foreach ( reverse @nums ) {
        if( /\d+/ ) {
            $_ = "-";
            last;
        }
    }
    join(".", @nums);
}

sub parse_enzclass {
    my ($file) = @_;
    my %results = ();

    open(my $in, "<$file") or die("Could not open $file ($!)");
    while( <$in> ) {
        chomp;
        #only parse lines that start with ec numbers
        if( /^(([\d\s-]+\.){3}[\d\s-]+)\s+(.*)$/ ) {
            my $num = $1;
            my $def = $3;
            $num =~ s/\s//g;
            $results{$num} = $def;
        }
    }
    close($in);
    return %results;
}


sub setup_ontology {
    my $ontology = new OBO::Core::Ontology();
    $ontology->default_namespace('EC');
    $ontology->add_relationship_type_as_string("is_a", "is_a");
    return $ontology;
}

sub parse_enzyme_dat {
    my ($file) = @_;
    my $enzyme = {};
    my $flag = 0;
    my $action;
    open(my $in, "<$file") or die("Could not open $file ($!)");
    while( <$in> ) {
        chomp;
        if( !$flag && m|^//| ) {
            $flag = 1;
        } elsif( $flag ) {
            if( m|^//| ) {
                #this is where we should add the info
                &add_enzyme_to_ontology( $enzyme, $ontology, $action );
                undef $action;
                $enzyme = {};
            } elsif( /^ID\s+(.*)/ ) {
                $enzyme->{'xref'} = $1;
            } elsif( /^DE\s+(.*)/ ) {
                my $name = $1;
                if( $name =~ /Transferred entry:\s+(.*)/ ) {
                    $action = "transferred";
                    my $transfer_to = $1;
                    $transfer_to =~ s/\.$//;
                    my @ids;
                    while( $transfer_to =~ /(([\d-]+\.){3}[\d-]+)/g ) {
                        push(@ids, $1)
                    }
                    $enzyme->{'transfer_to'} = \@ids;
                } elsif( $name =~ /Deleted entry/ ) {
                    $action = "deleted";
                }
                $enzyme->{'name'} = $name;
            } elsif( /^AN\s+(.*)/ ) {
                $enzyme->{'synonyms'} = [] unless( exists( $enzyme->{'synonyms'} ) );
                push(@{$enzyme->{'synonyms'}}, $1 );
            } elsif( /^CA\s+(.*)/ ) {
                $enzyme->{'definition'} .= " " if( exists( $enzyme->{'definition'} ) );
                $enzyme->{'definition'} .= $1;
            } elsif( /^CC\s+(.*)/ ) {
                $enzyme->{'comment'} .= " " if( exists( $enzyme->{'comment'} ) );
                $enzyme->{'comment'} .= $1;
            }
        } 
    }
    close($in);
}

sub add_enzyme_to_ontology {
    my ($einfo, $ont, $action) = @_;

    #create the term
    my $term = new OBO::Core::Term;
    $term->id( &next_id );
    $term->name( $einfo->{'name'} );
    $term->xref_set_as_string("[$xref_db:".$einfo->{'xref'}."]");
    $term->def_as_string($einfo->{'definition'}, "[]" ) if( exists( $einfo->{'definition'} ) );
    if( exists( $einfo->{'synonyms'} ) ) {
        map { $term->synonym_as_string( $_, "[]", "RELATED" ); } @{$einfo->{'synonyms'}};
    }
    $term->comment($einfo->{'comment'} ) if( exists( $einfo->{'comment'} ) );

    if( defined( $action ) && ( $action eq 'transferred' || $action eq 'deleted' ) ) {
        $term->is_obsolete( 1 );
    }

    if( defined( $action ) && $action eq 'transferred' ) {
        push( @transferred, { 'from' => $term, 'to' => $einfo->{'transfer_to'} } );
    }

    $ont->add_term($term);

    my $parent_ec = &get_direct_parent_ec( $einfo->{'xref'} );
    die("Could not find parent_ec for $einfo->{'xref'}") unless(defined( $parent_ec));
    my $parent_term = $ont->get_term_by_xref( $xref_db, $parent_ec );

    if( defined( $parent_term ) ) {
        my $relationship = new OBO::Core::Relationship();
        $relationship->type('is_a');
        $relationship->link($term, $parent_term);
        $relationship->id( $term->id()."_is_a_".$parent_term->id() );
        $ont->add_relationship( $relationship );
    }
}

sub add_replaced_by {
    my ($transferred, $ont) = @_;

    foreach my $t ( @{$transferred} ) {
        my $from_term = $t->{'from'};
        my @ids;
        foreach my $xref ( @{$t->{'to'}} ) {
            my $to_term = $ont->get_term_by_xref( $xref_db, $xref );
            die("Could not find term which replaces another term: [$t->{'to'}]") 
                unless( defined( $to_term ) );
            push( @ids, $to_term->id );
        }
        $from_term->replaced_by( @ids );
    }
}

sub next_id {
    my $id = $id_prefix.":".sprintf("%08d", $id_num++);
    return $id;
}
