package PFunc::EvidenceParser;

=head1 NAME

EvidenceParser.pm - used to parse evidence files for use in Annotation

=head1 SYNOPSIS

    my $parser = new EvidenceParser( 'evidence_file' => "/path/to/evidence.file",);

    or

    my $parser = new EvidenceParser( 'evidence_list' => "/path/to/evidence_files.list",
                                     'bsml_file' => "/path/to/file.bsml",
                                     'bsml_list' => "/path/to/bsml.list",
                                     'annotate_on' => 'transcript', );

=head1 DESCRIPTION

=head2 Overview

    Parent class for specific EvidenceParsers.  

=head2 Constructor and initialization

    At the least, an evidence file is required.  The format of the evidence
    bsml_file option specifies a bsml file to gather feature relationship information from. This
    is used when the annotate_on option is used.  For example, an input evidence set might contain
    information linked to polypeptide ids when this information linked to transcript models is preferred.
    In this case, the bsml would be used to lookup which polypeptide ids are related to transcript ids.

    Contructor options:
    evidence_file | evidence_list :: File containing evidence to be parsed. One of the two options required
    bsml_file | bsml_list         :: File containing evidence to be parsed. Required if annotate_on option is used
    annotate_on                   :: The feature type to apply annotations to. 

    

=head2 Class and object methods

=over 4

=cut


use strict;
use warnings;
use Carp;
use PFunc::Annotation;
use File::OpenFile qw( open_file );
use Data::Dumper;
use PFunc::EvidenceParser::ValidEvidenceTypes;
use BSML::FeatureRelationshipLookup;
$|++;
our %valid_evidence_types;

sub new {
    my ($class, %args) = @_;

    my $self = {
        '_evidence' => undef,
        '_bsml' => undef,
        '_lookup' => undef,
        '_evidence_type' => undef,
        '_features' => undef,
        '_annotate_on' => undef,
    };

    bless($self, $class);
    $self->_init(%args);
    return $self;
}

## adds evidence files to the _evidence arrayref
sub evidence {
    my ($self, @evidence) = @_;
    push( @{$self->{'_evidence'}}, @evidence ) if( @evidence  );
    return $self->{'_evidence'};
}

## will clear all input
sub clear_evidence {
    my ($self) = @_;
    delete $self->{'_evidence'};
    $self->{'_evidence'} = undef;
}

## bsml
sub bsml {
    my ($self, @bsml) = @_;
    push( @{$self->{'_bsml'}}, @bsml ) if( @bsml );
    return $self->{'_bsml'};
}

sub get_feature_annotation {
    my ($self, $feature_id) = @_;
    my $annotation;

    if( exists( $self->{'_features'}->{$feature_id}->{'_annotation'} ) ) {

        $annotation = $self->{'_features'}->{$feature_id}->{'_annotation'};

    } else {

        ## if the feature doesn't exist in the lookup, don't create it. If the lookup doesn't exist, 
        ## make whatever gets passed in.
        if( ( defined( $self->{'_lookup'} ) && $self->lookup_feature_id( $feature_id, $self->annotate_on ) ) ||
            !defined( $self->{'_lookup'} ) ) {
            $annotation = $self->store_annotation( new PFunc::Annotation( 'feature_id' => $feature_id ) );
        }

    }

    return $annotation;
}

sub store_annotation {
    my ($self, $annotation) = @_;
    my $feature_id = $annotation->get_feature_id;

    $self->{'_features'}->{$feature_id}->{'_annotation'} = 
        $annotation;

    return $annotation;
}

sub get_annotations {
    my ($self) = @_;
    my @annotations = ();
    map { push( @annotations, $self->{'_features'}->{$_}->{'_annotation'} ) 
              if( exists( $self->{'_features'}->{$_}->{'_annotation'} ) ) } keys %{$self->{'_features'}};
    return @annotations;
}

sub get_feature_ids {
    my ($self) = @_;
    return keys %{$self->{'_features'}};
}

sub annotate_on {
    my ($self, $feature_type) = @_;
    if( $feature_type ) {
        $self->{'_annotate_on'} = $feature_type;
    }
    return $self->{'_annotate_on'};
}

sub lookup_feature_id {
    my ($self, $id, $class) = @_;

    my $fl = $self->{'_lookup'};
    my $retval;

    if( defined( $fl ) && defined( $class ) ) {
        eval {
            $retval = $fl->lookup( $id, $class );
        };
        if( $@ ) {
            $retval = undef;
        }
    } else {
        $retval = $id;
    }

    return $retval;
}


sub to_file {
    my ($self, $file, $append) = @_;
    
    my $out = &open_file( $file, 'out' );

    foreach my $feature_id ( $self->get_feature_ids ) {
        foreach my $annotation ( $self->get_feature_annotation( $feature_id ) ) {
            next unless( $annotation->has_annotation );
            print $out $annotation->to_string()."\n";
        }
    }

    close($out);
    
}


## parses the evidence
sub parse {
    my ($self) = @_;
    my $evidence = $self->evidence();

    ## do we need to parse the feature_lookup bsml?
    if( !defined($self->{'_lookup'}) && $self->bsml ) {
        ## parse the bsml files to create id lookup
        print "Parsing bsml feature lookup\n";
        my $fr = new BSML::FeatureRelationshipLookup( 'bsml' => $self->bsml );
        $self->{'_lookup'} = $fr;
    }

    print "Pre parse\n";
    $self->_pre_parse;

    my $total = scalar( @{$evidence} );
    my $count = 0;

    print "starting evidence parse\n";
    foreach my $ev ( @{$evidence} ) {
        my $in = &open_file( $ev, 'in' );
        $self->_parse( $in );
        $count++;
        print "\r$count/$total";
        close( $in );
    }
}

############ private subroutines ################
sub _parse {
    my ($self) = @_;
    return undef;
}

sub _pre_parse {
    my ($self) = @_;
    return undef;
}

sub _init {
    my ($self, %args) = @_;

    if( $args{'evidence_file'} ) {
        croak("File $args{'evidence_file'} does not exist")
            unless( -e $args{'evidence_file'} );
        $self->evidence( $args{'evidence_file'} );
    }

    if( $args{'evidence_list'} ) {
        my $in = &open_file( $args{'evidence_list'}, 'in' );
        chomp( my @files = <$in> );
        close($in);
        $self->evidence( @files );
    }

    if( $args{'evidence'} ) {
        my $ev = $args{'evidence'};
        if( ref($ev) eq 'ARRAY' ) {
            $self->evidence( @{$ev} );
        }
    }

    if( $args{'bsml_list'} ) {
        my $in = &open_file( $args{'bsml_list'}, 'in' );
        chomp( my @bsml_files = <$in> );
        close( $in );
        $self->bsml( @bsml_files );
    }
    
    if( $args{'bsml_file'} ) {
        croak("File $args{'bsml_file'} does not exist")
            unless( -e $args{'bsml_file'} );
        $self->bsml( $args{'bsml_file'} );
    }
            
    if( $args{'annotate_on'} ) {
        $self->annotate_on( $args{'annotate_on'} );
    } else {
        $self->annotate_on( 'transcript' );
    }
}

1;
