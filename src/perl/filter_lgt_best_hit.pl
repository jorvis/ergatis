#!/usr/local/bin/perl

=head1 NAME

filter_lgt_best_hit.pl

=head1 SYNOPSIS

 USAGE: filter_lgt_best_hit.pl
       --input=/path/to/intput_btab
       --output1=/path/to/output_btab1
       --lineage1=a;lineage;to;filter;on
       --lineage1tophit
       --output2=/path/to/output/btab2
       --lineage2=another;lineage;to;filter;on
       --lineage2tophit
      [--help]
     
=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use GiTaxon;
use List::Util qw[min max];
use File::Basename;
$| =1;

my %options;
my $results = GetOptions (\%options,
                          'input=s',
                          'trace_mapping=s',
                          'output1=s',
                          'output2=s',
                          'lineage1=s',
                          'lineage1tophit',
                          'lineage2=s',
                          'lineage2tophit',
                          'overalloutput=s',
                          'filter_min_overlap=s',
                          'filter_lineage=s',
                          'ncbitax=s',
                          'gitax=s',
                          'dbhost=s',
                          'taxondb=s',
                          'taxoncoll=s',
                          'idx_dir=s',
                          'output_dir=s',
                          'help|h'
                          );

my $filter_lineage = $options{filter_lineage};

my $FILTER_MIN_OVERLAP = $options{filter_min_overlap} ? $options{filter_min_overlap} : 50;
open IN, "<$options{input}" or die "Unable to open input $options{input}\n";
my($name,$directories,$suffix) = fileparse($options{input},qr/\.[^.]*/);
my $overallfile = $options{overalloutput};
if(!$options{overalloutput} && $options{output_dir}) {
    $overallfile = $options{output_dir}."/$name\_overall.out";
    print STDERR "Printing to $overallfile\n";
}

open my $out3, ">$overallfile" or die "Unable to open output $options{overalloutput}\n";
my $map_file = $options{trace_mapping};
my $trace_lookup = {};

my $clones_by_clone_id = {};

if($map_file) {
    print STDERR "Reading the map file\n";
    open MAP, "<$map_file" or die "Unable to open $map_file\n";
    while(<MAP>) {
        my @fields = split;
        $trace_lookup->{$fields[0]} = 
        {'template_id' => $fields[2],
         'trace_end' => $fields[1]
        };
        if(!$clones_by_clone_id->{$fields[2]}) {
            $clones_by_clone_id->{$fields[2]} = {
                'forward'=>{},
                'reverse'=> {}
        };
        }
    }

    close MAP;
    print STDERR "Done reading the map file\n";
}
else {
    print STDERR "No map file, will try and assume mates.\n";

}

my $gi2tax = &get_gi2taxon();

my $id;
my $best_e1 = 100;
my $best_e2 = 100;
my $best_rows1 = [];
my $best_rows2 = [];
my $best_e = 100;
my $best_overall = [];
my $filter_hits = [];

my $lineages = [];
if($options{lineage1} && $options{lineage2}) {

    my $out1file = $options{output1};
    if(!$options{output1} && $options{output_dir}) {
        $out1file = $options{output_dir}."/$name\_lineage1.out";
    }
    open my $out1, ">$out1file" or die "Unable to open output $options{output1}\n";

    my $out2file = $options{output2};
    if(!$options{output1} && $options{output_dir}) {
        $out1file = $options{output_dir}."/$name\_lineage2.out";
    }
    open my $out2, ">$out2file" or die "Unable to open output $options{output2}\n";

    push(@$lineages,({
        'lineage' => $options{lineage1},
        'tophit' => 1,#$options{lineage1tophit},
        'best_e' => 100,
        'id' => '',
        'handle' => $out1,
        'best_rows' => [],
        'name' => 'lineage1'
       },
       {'lineage' => $options{lineage2},
        'tophit' => 1,#$options{lineage2tophit},
        'best_e' => 100,
        'id' => '',
        'handle' => $out2,
        'best_rows' => [],
        'name' => 'lineage2'
       }));
}
push(@$lineages,{'lineage' => '',
     'tophit' => 1,
     'best_e' => 100,
     'id' => '',
     'handle' => $out3,
     'best_rows' => [],
     'name' => 'overall'
     });


&process_file($options{input});

sub process_file {
    my $file = shift;
    open IN, "<$file" or die "Unable to open $file\n";


    while(<IN>) {
        chomp;
        my @new_fields = split(/\t/, $_);

        my $tax;

        my $found_tax = 0;
        # If we already have lineage info in here we'll not append it again
        if($new_fields[14]) {
            $found_tax = 1;
            $tax = {'taxon_id' => $new_fields[12],
                    'lineage' => $new_fields[14],
                    'name' => $new_fields[13]
            };
        }
        else {
            $tax = $gi2tax->getTaxon($new_fields[1]);
        }

   

        if($tax->{'taxon_id'}) {
            if(!$found_tax) {
                push(@new_fields,$tax->{'taxon_id'});
            }
            if($tax->{'name'}) {
                if(!$found_tax) {
                    push(@new_fields,($tax->{'name'},$tax->{'lineage'}));
                }
                my $fields = &add_trace_info(\@new_fields);

                my $done = 0;
                foreach my $lineage (@$lineages) {
#                    print STDERR "$fields->[10]\n";
                    $done = &process_line($fields,$tax,$lineage);
                    if($done ==1) {
                        &append_hits($lineage->{best_rows},$lineage->{name});
                    }
                    elsif($done eq 'filter') {
                        print "Filtered $tax->{name} $fields->[0]\n";
                    }
                    
#                if(&process_line(\@fields,$best_e1,$best_rows1,$out1,$options{lineage1tophit},$options{lineage1})) {
#                    $best_e1=100;
#                    &append_hits($best_rows1,'lineage1');
#                    print STDERR "$fields[0] found ".scalar @$best_rows1." lineage1\n";
#                }
#                if(&process_line(\@fields,$best_e2,$best_rows2,$out2,$options{lineage2tophit},$options{lineage2})) {
#                    $best_e2 =100;
#                    &append_hits($best_rows2,'lineage2');
#                    print STDERR "found ".scalar @$best_rows2." lineage2\n";
#                }
#                if(&process_line(\@fields,$best_e,$best_overall,*STDOUT,1,'')) {
#                    $best_e = 100;
#                    &append_hits($best_overall,'best');
#                    print STDERR "found ".scalar @$best_overall." best overall\n";
#                    $id = $fields[0];
#                }
                    
                }
                if($done) {
                    $filter_hits = [];
                }
            }
            else {
                print STDERR "Unable to find name or lineage for taxon_id $tax->{'taxon_id'} in trace $new_fields[1]\n";
            }
        }
        else {
            print STDERR "Unable to find taxon info for $new_fields[1]\n";
        }
    }
    close IN;
    
    # here we'll take care of the last trace in the file.
    foreach my $lineage (@$lineages) {
        &append_hits($lineage->{best_rows},$lineage->{name});
        &print_hits($lineage);
    }
}

sub get_gi2taxon {
    my $ncbitax = $options{ncbitax} ? $options{ncbitax} : '/local/db/repository/ncbi/blast/20110831_122720/taxonomy/taxdump/';
    my $gi2tax = $options{gitax} ? $options{gitax} : '/local/db/repository/ncbi/blast/20110831_122720/taxonomy/gi_taxid_nucl.dmp';
    my $dbhost = $options{dbhost} ? $options{dbhost} : 'tettelin-lx.igs.umaryland.edu';
    my $taxondb = $options{taxondb} ? $options{taxondb} : 'gi2taxon';
    
    my $taxoncoll = $options{taxoncoll};
    if(!$taxoncoll) {
        $ncbitax =~ /(\d+\_\d+)/;
        my $date = $1;
    if($gi2tax =~ /nuc/) {
        $taxoncoll = "gi2taxonnuc_$date";
    }
        else {
            $taxoncoll = "gi2taxonprot_$date";
        }
        
    }
    
    my $idx_dir = $options{idx_dir};
    
    if(!$idx_dir && -e "$ncbitax/names") {
        $idx_dir = $ncbitax;
    }
    elsif(!$idx_dir) {
        $idx_dir='/tmp/';
}
    
    my $gi2tax = GiTaxon->new(
        {'nodes' => "$ncbitax/nodes.dmp",
         'names' => "$ncbitax/names.dmp",
         'gi2tax' => $gi2tax,
         'chunk_size' => 10000,
         'idx_dir' => $idx_dir,
         'host' => $dbhost,
     'gi_db' => $taxondb,
         'gi_coll' => $taxoncoll
        });
    
    return $gi2tax;
}
sub append_hits {
    my $hits = shift;
    my $list_name = shift;


    my $template_id = $hits->[0]->[12];
    
    my $clone = $clones_by_clone_id->{$template_id};

    map {
        # Check what strand we're on
        if($_->[13] eq 'F') {

            # Create the trace entry
            if(!$clone->{'forward'}->{$_->[0]}) {
                $clone->{'forward'}->{$_->[0]} = {
                    'trace_id' => $_->[0],
                    "$list_name" => []
                };
            }
            if(!$clone->{'forward'}->{$_->[0]}->{$list_name}) {
                $clone->{'forward'}->{$_->[0]}->{$list_name} = [];
            }
            push(@{$clone->{'forward'}->{$_->[0]}->{$list_name}},{
                'accession' => $_->[1],
                'pid' => $_->[2],
                
                 });
        }
        elsif($_->[13] eq 'R') {

            # Create the trace entry
            if(!$clone->{'reverse'}->{$_->[0]}) {
                $clone->{'reverse'}->{$_->[0]} = {
                    'trace_id' => $_->[0],
                    "$list_name" => []
                };
            }
            if(!$clone->{'reverse'}->{$_->[0]}->{$list_name}) {
                $clone->{'reverse'}->{$_->[0]}->{$list_name} = [];
            }
            push(@{$clone->{'reverse'}->{$_->[0]}->{$list_name}},{
                'accession' => $_->[1],
                'pid' => $_->[2],
                
                 });
        }
    } @$hits;
}

sub process_line {
    my ($fields,$tax,$lineage) = @_;

#    print STDERR "Looking for $lineage in $fields->[16]\n";
    my $finished_id = 0;
    # See if we are on the last id still
    if($lineage->{id} eq $fields->[0]) {
        if($lineage->{tophit}) {
#            print STDERR "here with $fields->[16]\n";
            # If this is one of the best
            if($fields->[10] == $lineage->{best_e}) {
                # Check that we are one of the lineage
                if($tax->{lineage} =~ /$lineage->{lineage}/) {
                    push(@{$lineage->{best_rows}}, $fields);
                }
            }
            # If we find a lower e-value
            elsif($fields->[10] < $lineage->{best_e}) {
                # Check that we are one of the lineage
                if($tax->{lineage} =~ /$lineage->{lineage}/) {
                    $lineage->{best_e} = $fields->[10];
                    $lineage->{best_rows} = [$fields];
                }
            }
        }
        # If we aren't filtering for best hits...
        else {
            # Not doing it for now.
        }
    }
    # If we are no longer on the last id
    else {
        if($lineage->{tophit}) {
#            if($lineage->{id}) {
                # We have finished a hit
                &print_hits($lineage);
                $finished_id = 1;
#                if(scalar @{$lineage->{best_rows}}) {
#                    print STDERR "Should have printed to $lineage->{name}\n";
#                    map {
#                        print {$lineage->{handle}} join("\t", @$_);
#                        print {$lineage->{handle}} "\n";
#                    }@{$lineage->{best_rows}};
#                }
#            }
            
#            $lineage->{best_e} = $fields->[10];
#            $id = $fields->[0];
            # Check that we are one of the lineages we want
            if($tax->{lineage} =~ /$lineage->{lineage}/) {
#                print "Resetting the lineage to $fields->[0] $lineage->{name}\n";
                $lineage->{id} = $fields->[0];
                $lineage->{best_e} = $fields->[10];
                $lineage->{best_rows} = [$fields];
            }
            else {
                $lineage->{best_rows} = [];
                $lineage->{id} = undef;
                $lineage->{best_e} = 100;
            }
        }
        else {
            # Not doing it for now
        }
    }

    if(&filter_hit($tax->{name})) {
        push(@$filter_hits, $fields);
    }
    return $finished_id;
}

sub print_hits {
    my $lineage = shift;
    if($lineage->{id}) {
        # We have finished a hit
        
        if(scalar @{$lineage->{best_rows}}) {
            if(!&filter_best_hits($lineage->{best_rows})) {
                
#                print STDERR "Should have printed to $lineage->{name}\n";
                map {
                    print {$lineage->{handle}} join("\t", @$_);
                    print {$lineage->{handle}} "\n";
                }@{$lineage->{best_rows}};
            }
            else {
                print STDERR "Filtered $lineage->{name}\n";
            }
        }
    }
}

sub add_trace_info {
    my $list = shift;

    my $id = $list->[0];
    if($trace_lookup->{$list->[0]}) {
        push(@$list,($trace_lookup->{$id}->{'template_id'},$trace_lookup->{$id}->{'trace_end'}));
    }
    elsif(!$options{trace_mapping}) {
        # Ghetto way of checking for directionality.
        $list->[0] =~ /(.*)\_(\d)/;
        if($2 == 1) {
            push(@$list, ($1,'F'));
        }
        elsif($2 ==2) {
            push(@$list, ($1,'R'));            
        }
        else {
            print STDERR "Couldn't figure out the clone name from $list->[0]\n";
        }
    }
    else {
        print STDERR "couldn't find trace info for $list->[0]\n";
    }
    return $list;
}

sub filter_best_hits {
    my $hits = shift;
    my $filter = 0;
    foreach my $fhit (@$filter_hits) {
        foreach my $hit (@$hits) {
            my $overlap = &get_overlap_length([$fhit->[6],$fhit->[7]],[$hit->[6],$hit->[7]]);
            if($overlap >=$FILTER_MIN_OVERLAP) {
                print STDERR "Here to filter out $hit->[0] $fhit->[14] with $overlap\n";
                $filter = 1;
                last;
            }
        }
        last if $filter;
    }
    $filter;
    
}

sub filter_hit {
    my $lineage = shift;

    my $retval = 0;
    if($filter_lineage && $lineage =~ /$filter_lineage/i) {
        $retval = 1;
    }
    return $retval;
}

sub get_overlap_length {
    my ($q, $s) = @_;

    my $qlen =  $q->[1] - $q->[0];
    my $slen = $s->[1] - $s->[0];
    my $min = min($q->[1],$q->[0],$s->[1],$s->[0]);
    my $max = max($q->[1],$q->[0],$s->[1],$s->[0]);

    my $len = ($qlen + $slen) - ($max - $min);

    return $len;
}
