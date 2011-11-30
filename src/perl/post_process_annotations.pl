#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1 NAME

post_process_annotations.pl - Contains functionality to refine final annotations.

=head1 SYNOPSIS

 USAGE: post_process_annotations.pl
       --input_file=/path/to/file.tab
       --output=/path/to/output.tab
       --tigr_roles_db_dir=/path/to/db_dir
     [ --input_ec_dat=/path/to/enzyme.dat       
       --hypo="conserved hypothetical protein" or "hypothetical protein"
       --log=/path/to/file.log
       --debug=4
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>
    The input tab file to parse.  Should be output from assign_annotations.pl. In general, a file containing
    Annotation.pm objects ->to_string representations

B<--output,-o>
    Output tab file.

B<--tigr_roles_db_dir,-b>
    The path to the database directory. The directory should contain the file:
    tigr_roles/tigr_roles_keywords.txt.

B<--input_ec_dat,-e>
    Path to the enzyme.dat containing EC numbers file downloaded from 
    [ftp://ftp.expasy.org/databases/enzyme/enzyme.dat]

B<--hypo,-h>
    Default is hypothetical protein. Option is conserved hypothetical protein

B<--log,-l>
    Logfile.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

    Will process various annotations assigned to the features including:

    1. Name clean up as described in CommonNameAnnotation.pm
    2. TIGR role assignment by common name keywords
         DESCRIPTION or link to module doing the work
 
=head1  INPUT

    A file based on the Annotation->to_string implementation.  See Annotation.pm for details.  

=head1 OUTPUT

    Same format file, with the annotations changed where described in description

=head1  CONTACT

    Kevin Galens
    kgalens@som.umaryland.edu

=cut


use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::OpenFile qw(open_file);
use PFunc::Annotation;
use PFunc::EvidenceParser::Hypothetical;
use PFunc::CommonNameAnnotation qw(clean_common_name clean_gene_symbol);
use PFunc::TIGRRolesAnnotation;
use Data::Dumper;

##### Globals #####
my $input_file;
my $output;
my $tigr_roles_annot;
my $input_ec_dat;
my $chp;
###################

&check_options();
my $enzymes = &parse_dat_file($input_ec_dat) if (defined $input_ec_dat);
my $in = open_file( $input_file, 'in' ); 
my $out = open_file( $output, 'out' ); 

while( my $line = <$in> ) {
    chomp( $line );

    my $annotation = PFunc::Annotation::to_object( $line );

    &post_process( $annotation );
}

close($in);
close($out);

sub post_process {
    my ($annotation) = @_;

    ## check for ec numbers in the common name and add them as annotations
    &check_for_ec_numbers( $annotation );

    ## format ec numbers
    &format_ec_numbers( $annotation );

    ## clean up the common name  with the CommonNameAnnotation module
    my ($gene_product_name, $gpn_source, $gpn_source_type) = $annotation->get_gene_product_name;
    $annotation->set_gene_product_name( clean_common_name( $gene_product_name->[0],$chp ),
                                        $gpn_source, $gpn_source_type );
    my ($gene_symbol, $gene_source, $gene_source_type) = $annotation->get_gene_symbol;
    if(defined($gene_symbol->[0])) {
    	$annotation->set_gene_symbol(clean_gene_symbol($gene_symbol->[0]), $gene_source, $gene_source_type);
    }
    
    ($gene_product_name, $gpn_source, $gpn_source_type) = $annotation->get_gene_product_name;

    &correct_gene_symbol($annotation, $gene_product_name->[0]);
    ## Assign changed GO terms to conserved hypothetical proteins
    &correct_go_terms($annotation, $gene_product_name->[0]);
    &correct_ec_number($annotation, $gene_product_name->[0]);
    ## can we assign any TIGR roles based on common name keyword? It helps if this is after the
    ## clean up function
    my ($tigr_roles, $tr_source, $tr_source_type) = $annotation->get_TIGR_Role;
    my $new_tigr_roles;
    if( $new_tigr_roles = $tigr_roles_annot->assign_tigr_roles_by_keyword( $gene_product_name->[0], $tigr_roles ) ) {
        $annotation->set_TIGR_Role( $tigr_roles_annot->assign_tigr_roles_by_keyword( $gene_product_name->[0], $tigr_roles ),
                                    'by_keyword', 'by_keyword' );
    }

    print $out $annotation->to_string."\n";
}

sub format_ec_numbers {
    my ($annotation) = @_;
    my ($ecs, $source, $type) = $annotation->get_EC;
    my @new_ecs;
    my $change = 0;

    foreach my $ec ( @{$ecs} ) {
        if( $ec =~ /(([\d\-]+\.){3}[\d\-]+)/ ) {
            my $tmp_ec = $1;
            if( $ec ne $tmp_ec ) {
                $change = 1;
                $ec = $tmp_ec;
            }
	    if(defined $input_ec_dat) { 
	    	next if( $ec =~ /-/ );
	   	if( exists( $enzymes->{'obsolete'}->{$ec} ) ) {
			$change = 1;
#			print  "Change $ec to $enzymes->{'obsolete'}->{$ec}\n";
			$ec = $enzymes->{'obsolete'}->{$ec};
            	} elsif( !exists( $enzymes->{'current'}->{$ec} ) ) {
            		$change = 1;
#			print "Remove from $ec\n";
         		$ec = "";
	    	}
	    }
	    push( @new_ecs, $ec );
        }
    }

    if( $change ) {
        $annotation->set_EC(\@new_ecs, $source, $type);
    }
}

sub check_for_ec_numbers {
    my ($annotation) = @_;

    my ($gene_product_name_array, $source, $type) = $annotation->get_gene_product_name( 'gene_product_name' );
    my $gene_product_name = $gene_product_name_array->[0];

    if( !defined( $gene_product_name ) ) {
        &PFunc::EvidenceParser::Hypothetical::_set_annotation( $annotation );
    } else {

        my $feature_id = $annotation->get_feature_id;
    
        while( $gene_product_name =~ /ec\s+([\d\.\-]+)/gi ) {
            
            my $ec = $1;
            if( $annotation->has_annotation( 'EC' ) ) {
                $annotation->add_EC( $ec );
            } else {
                $annotation->set_EC( $ec, $source, $type );
            }
        }
    }
}

sub correct_go_terms {
    my ($annotation, $gene_product) = @_;
## In case where protein is 'hypothetical protein' assign 'GO:0008150', 'GO:0003674', 'GO:0005575 GO terms    
    if ($gene_product eq $chp) {
    	my ($go_terms, $go_source, $go_source_type) = $annotation->get_GO;
    	my @new_go = ('GO:0008150', 'GO:0003674', 'GO:0005575');
    	$annotation->set_GO(\@new_go, $go_source, $go_source_type);
     }
}

sub correct_gene_symbol {
	my ($annotation, $gene_product) = @_;
	my ($gene_symbol, $gene_source, $gene_source_type) = $annotation->get_gene_symbol;
# Gene symbol should be blank if gene product name is "conserved hypothetical protein" or begins with putative 	
	if (($gene_product eq $chp) || ($gene_product eq "hypothetical protein")) {
		$annotation->set_gene_symbol("", $gene_source, $gene_source_type);		
	} 
	if($gene_product =~ /^putative/i) {
		$annotation->set_gene_symbol("", $gene_source, $gene_source_type);
	}
}

sub correct_ec_number {
	my ($annotation, $gene_product) = @_;
	my ($ecs, $source, $type) = $annotation->get_EC;
# EC number should be blank if gene product name is "conserved hypothetical protein" or begins with putative 	
	if (($gene_product eq $chp) || ($gene_product eq "hypothetical protein")) {
		$annotation->set_EC("", $source, $type);
	} 
	if($gene_product =~ /^putative/i) {
		$annotation->set_EC("", $source, $type);
	}
}

sub parse_dat_file {
    my ($file) = @_; 
    open(IN, "< $file") or die("Could not open $file: $!");

    my %data = ('current'=>{},'obsolete'=>{});
    my ($current, $transfer, $def, $name);
    while( my $line = <IN> ) { 
        chomp($line);
        if( $line =~ m|^//| ) { 

            if( defined( $current )  ) { 
                ## we should store current
                if( $transfer ) { 
                    $data{'obsolete'}->{$current} = $transfer;
                } else {
                    $data{'current'}->{$current} = { 
                        'name' => $name,
                        'def' => $def
                        };
                }
            }

            ## reset the variables
            map { undef($_); } ($current, $transfer, $def, $name);

        } elsif( $line =~ /^ID\s+(.*)/ ) { 
            $current = $1; 
        } elsif( $line =~ /^DE\s+(.*)/ ) { 
            my $deline = $1; 
            if( $deline =~ /Transferred entry:\s+(\S*)/ ) {
		if($deline =~ /and|\,/) {
			map { undef($_); } ($current, $transfer, $def, $name);
		} else {
                	$transfer = $1; 
	   		$transfer =~ s/\.$//;
		}
            } elsif($deline =~ /Deleted entry/) {
		map { undef($_); } ($current, $transfer, $def, $name);
	    } else {
                if(defined($current)) {
			$name = $deline." [$current]";
		}
            }
        } elsif( $line =~ /^CA\s+(.*)/ ) { 
            if(defined($def)) {
                $def.= " ".$1;
            } else {
                $def = $1; 
            }
        }
    }   
    close(IN);
    return \%data;
}

sub check_options {
     my %options;
     my $results = GetOptions (\%options,
                              'input_file|i=s',
                              'output|o=s',
                              'tigr_roles_db_dir|d=s',
			      'input_ec_dat|e=s',
			      'hypo|h=s',
                              'log|l=s',
                              'debug|d=s',
                              'help|h',
                              );

    if( $options{'help'} ) {
        &_pod;
    }

    if( $options{'tigr_roles_db_dir'} ) {
        $tigr_roles_annot = new PFunc::TIGRRolesAnnotation( 'tigr_roles_db_dir' => $options{'tigr_roles_db_dir'} );
    } else {
        die("Option --tigr_roles_db_dir is required");
    }

    if( $options{'input_file'} ) {
        $input_file = $options{'input_file'};
    } else {
        die("Option --input_file is required");
    }

    if( $options{'output'} ) {
        $output = $options{'output'};
    } else {
        die("Option --output is required");
    }
    
    if( $options{'input_ec_dat'} ) { 
        $input_ec_dat = $options{'input_ec_dat'};
    } 
 
    if($options{'hypo'}) {
	$chp = $options{'hypo'};
    } else {
	$chp = 'conserved hypothetical protein';
    }	

}
