#!/usr/local/bin/perl

=head1 NAME

lgt_finder.pl

=head1 SYNOPSIS

 USAGE: lgt_finder.pl
       --input=/path/to/file/list
      [--help]
     
=cut

use strict;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use List::Util qw[min max];
$| =1;

my %options;
my $results = GetOptions (\%options,
                          'input=s', # Comma separated list of files
                          'input_file_list=s',
                          'output_dir=s',
                          'output_prefix=s',
                          'max_overlap=s',
                          'min_length=s',
                          'append_links=s',
                          'help|h'
                          );

my $traces_by_template = {};
my $filtered_clones = {};
my $traces_by_trace_id = {};

my $MAX_OVERLAP = $options{max_overlap} ? $options{max_overlap}: 0;
my $MIN_LENGTH = $options{min_length} ? $options{min_length}: 0;

my $output_dir = $options{output_dir};
my $filename = 'lgt_finder';

if($options{output_prefix}) {
    $filename = $options{output_prefix};
    
}
elsif($options{input_file_list}) {
    my($name,$directories,$suffix) = fileparse($options{input_file_list},qr/\.[^.]*/);
    $filename = $name;
    $output_dir = $options{output_dir} ? $options{output_dir} : $directories;
}
if($options{input}) {
    my($name,$directories,$suffix) = fileparse($options{input});
    $output_dir = $options{output_dir} ? $options{output_dir} : $directories;
}

print STDERR "$output_dir/$filename\_by_clone.txt\n";
open OUTCLONE, ">$output_dir/$filename\_by_clone.txt" or die "Couldn't open $output_dir";
open OUTTRACE, ">$output_dir/$filename\_by_trace.txt" or die "Couldn't open $output_dir";

if($options{input}) {

    my $old_prefix;
    foreach my $file (split(/,/,$options{input})) {
        chomp;
        my($fn,$directories,$suffix) = fileparse($file,qr/\.[^.]*/);
        my $prefix =$options{output_prefix};
        if(!$prefix) {
            $fn =~ /^(.*)\_[^\_]+/;
            $prefix = $1;
        }
        if($prefix && ($prefix ne $old_prefix)) {
            print STDERR "$prefix is the new prefix\n";
            &_find_lgt_in_trace();
            $traces_by_trace_id = {};
            $old_prefix = $prefix;
        }
        &process_file($file);
    }
    
    &find_lgt();
    close OUTCLONE;
    close OUTTRACE;        
}
elsif($options{input_file_list}) {
    open IN, "<$options{input_file_list}" or die "Unable to open input $options{input_file_list}\n";
    my($filename,$directories,$suffix) = fileparse($options{input_file_list},qr/\.[^.]*/);
    
    my $old_prefix;
    while(<IN>) {
        chomp;
        my($filename,$directories,$suffix) = fileparse($_,qr/\.[^.]*/);
        $filename =~ /^([^.]*)/;
        my $prefix = $1;
        if($prefix && ($prefix ne $old_prefix)) {
            print STDERR "$prefix is the new prefix\n";
            &_find_lgt_in_trace();
            $traces_by_trace_id = {};
            $old_prefix = $prefix;
        }
        &process_file($_);
    }
    
    &find_lgt();
    close OUTCLONE;
    close OUTTRACE;
}

sub process_file {
    my $file = shift;
    $file =~ /_([^_]+).out$/;
    my $type = $1;
    if($type eq 'overall') {
        &_process_overall($file);
    }
    else {
        &_process_within_trace($file,$type);
    }
}


sub _process_overall {
    my $file = shift;
    open IN2, "<$file" or die "Unable to open $file\n";
    while(<IN2>) {
        chomp;
        my @fields = split(/\t/, $_);

        my $hit = {
            trace_id => $fields[0],
            accession => $fields[1],
            pid => $fields[2],
            align_len => $fields[3],
            mismatches => $fields[4],
            gaps => $fields[5],
            query_start => $fields[6],
            query_end => $fields[7],
            subj_start => $fields[8],
            subj_end => $fields[9],
            evalue => $fields[10],
            score => $fields[11],
            taxon_id => $fields[12],
            scientific_name => $fields[13],
            lineage => $fields[14],
            template_id => $fields[15],
            direction => $fields[16]


        };
        
        my $dir = 'forward';
        if($hit->{direction} eq 'R') {
            $dir = 'reverse';
        }

        $hit->{scientific_name} =~ /^(\w+) /;
        my $genera = $1;
        if(!$traces_by_template->{$hit->{template_id}}) {
            $traces_by_template->{$hit->{template_id}} = {
                'forward' => {},
                'reverse' => {}
            };
        }
        if(!$traces_by_template->{$hit->{template_id}}->{$dir}->{$genera}) {
            $traces_by_template->{$hit->{template_id}}->{$dir}->{$genera} = [];
        }
        if(&filter_clone($hit->{scientific_name})) {
            $filtered_clones->{$hit->{template_id}} = 1;
        }
        if(!$filtered_clones->{$hit->{template_id}}) {
            if(!$traces_by_template->{$hit->{template_id}}->{$dir}->{hits}) {
                $traces_by_template->{$hit->{template_id}}->{$dir}->{hits} = [];
                }
            $traces_by_template->{$hit->{template_id}}->{$dir}->{traces}->{$hit->{trace_id}} = 1;
            $hit->{genera} = $genera;
            $traces_by_template->{$hit->{template_id}}->{$dir}->{hit} = {evalue => $hit->{evalue},
                                                                         align_len => $hit->{align_len},
                                                                         trace_id => $hit->{trace_id},
                                                                         genera => $hit->{genera}};
            
            push(@{$traces_by_template->{$hit->{template_id}}->{$dir}->{genera}->{$genera}},$hit->{trace_id});
        }
    }
    # This strategy was too memory intensive
    
#        if(!$traces_by_template->{$hit->{template_id}}) {
#            $traces_by_template->{$hit->{template_id}} = {
#                forward_hits => {},
#                reverse_hits => {}
#            };
#        }
#        my $trace = $traces_by_template->{$hit->{template_id}}->{$dir}->{$hit->{trace_id}};
#        if(!$trace) {
#            $traces_by_template->{$hit->{template_id}}->{$dir}->{$hit->{trace_id}} = {};
#            $trace = $traces_by_template->{$hit->{template_id}}->{$dir}->{$hit->{trace_id}};
#        }
#        if(!$trace->{$type}) {
#            $trace->{$type} = [];
#        }
#        push(@{$trace->{$type}}, $hit);
#    }
    close IN2;
}
sub find_lgt {

    if(keys %$traces_by_template) {
        &_find_lgt_in_clone();
    }
    if(keys %$traces_by_trace_id) {
        &_find_lgt_in_trace();
    }

}

sub filter_clone {
    my $scientific = shift;

    my $retval = 0;

    # Currently NOT doing this here
    #if($scientific =~ /vector/i) {
    #    $retval = 1;
    #}
    return $retval;
}
sub filter_genera {
    my $genera = shift;

    my $genera_to_ignore = {
        'Pan' => 1,
        'Papio' => 1,
        'Macaca' =>1
    };

    my $retval = 0;
    if($genera_to_ignore->{$genera}) {
        $retval = 1;
    }
    return $retval;
}

sub _find_lgt_in_trace {
    
    foreach my $trace (keys %$traces_by_trace_id) {
        
        my $found_lgt = 0;
        my $current_clone = '';
        my $genera = {};
        my $evalue = '';
        my $length = '';
        # Look at all combinations of genera to see if any have LGT evidence.
        foreach my $g1 (keys %{$traces_by_trace_id->{$trace}}) {
            foreach my $h1 (@{$traces_by_trace_id->{$trace}->{$g1}}) {
                foreach my $g2 (keys %{$traces_by_trace_id->{$trace}}) {
                    foreach my $h2 (@{$traces_by_trace_id->{$trace}->{$g2}}) {
                        if(($g1 ne $g2) && &check_overlap($h1,$h2) &&
                           ($h1->{align_len} >= $MIN_LENGTH) &&
                           ($h2->{align_len} >= $MIN_LENGTH) &&
                            !$filtered_clones->{$h1->{template_id}}) {
                            $current_clone = $h1->{template_id};
                            if($g1 =~ /Homo/) {
                                $genera->{"$g2:$g1"} = 1;
                                $found_lgt="$g2\t$g1";
                                $evalue = $h2->{evalue};
                                $length = $h2->{align_len};
                            }
                            else {
                                $genera->{"$g1:$g2"} = 1;
                                $found_lgt="$g1\t$g2";
                                $evalue=$h1->{evalue};
                                $length = $h1->{align_len};
                                    
                            }
                        }
#                        last if($found_lgt);
                    }
#                    last if($found_lgt);
                }
#                last if($found_lgt);
            }
#            last if($found_lgt);
        }
        if($found_lgt) {
#            print STDERR "Found LGT in trace $trace between ".join(',',keys %$genera)."\n";
            print OUTTRACE join("\t",($trace,$current_clone,$evalue,$length,join(',',sort(keys %$genera)),"http://driley-lx:8080/blast/$trace\_$current_clone.html","http://trace.ncbi.nlm.nih.gov/Traces/trace.cgi?&cmd=retrieve&val=TRACE_NAME=\%22$trace\%22\%20and\%20TEMPLATE_ID=\%22$current_clone\%22&dopt=trace&size=1&seeas=Show"));
            print OUTTRACE "\n";
        }
        $found_lgt = 0;
    }
}

sub check_overlap {
    my($hit1,$hit2) = @_;
    my $overlap = &get_overlap_length([$hit1->{query_start},$hit1->{query_end}],[$hit2->{query_start},$hit2->{query_end}]);
    
#    print STDERR "here with $hit1->{query_start} $hit1->{query_end} $hit2->{query_start} $hit2->{query_end} $overlap\n";
    
    my $retval = 0;
    if($overlap < $MAX_OVERLAP) {
        $retval =1;
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

sub _find_lgt_in_clone {
    my $lgt_by_genera = {};
    foreach my $clone (keys %$traces_by_template) {
        
#        my $forward_genera = {};
#        my $reverse_genera = {};

#        my $filter = 0;
#        foreach my $trace_id (keys %{$traces_by_template->{$clone}->{'forward_hits'}}) {
            
#            my $best_hits = $traces_by_template->{$clone}->{'forward_hits'}->{$trace_id}->{'overall'};

#            map {
#                $_->{'scientific_name'} =~ /^(\w+) /;
#                if(!&filter_genera($1)) {
#                    $forward_genera->{$1} = 1;
#                    if(&filter_clone($_->{'scientific_name'})) {
#                        $filter = 1;
#                    }
#                }
#            } @$best_hits;
#        }

#        foreach my $trace_id (keys %{$traces_by_template->{$clone}->{'reverse_hits'}}) {
            
#            my $best_hits = $traces_by_template->{$clone}->{'reverse_hits'}->{$trace_id}->{'overall'};

#            map {
#                $_->{'scientific_name'} =~ /^(\w+) /;
#                if(!&filter_genera($1)) {
#                    $reverse_genera->{$1} = 1;
#                    if(&filter_clone($_->{'scientific_name'})) {
#                        $filter = 1;
#                    }
#                }
#            } @$best_hits;
#        }

        if(!$filtered_clones->{$clone}) {
            my $forward = join(",", keys %{$traces_by_template->{$clone}->{forward}->{genera}});
            my $reverse = join(",", keys %{$traces_by_template->{$clone}->{reverse}->{genera}});
#            print $clone;
#            print "\tForward:";
#            print $forward;
#            print "\tReverse:";
#            print $reverse;
#            print "\n";
#            print STDERR "Here with $forward##$reverse\n";
            if($forward && $reverse) {

                my $fwd_traces = [];
                my $rev_traces = [];
                my $fhit = $traces_by_template->{$clone}->{forward}->{hit};
                my $nft = scalar keys %{$traces_by_template->{$clone}->{forward}->{traces}};
                push(@$fwd_traces, join("\t",($fhit->{trace_id},$fhit->{evalue},$fhit->{align_len})));
                my @fields = ($clone,'F',$nft,$forward,join(",",@$fwd_traces));
                if($options{append_links}) {
                    push(@fields,("driley-lx:8080/blast/$clone\_F.html","http://www.ncbi.nlm.nih.gov/Traces/trace.cgi?&cmd=retrieve&val=TEMPLATE_ID%20%3D%20%22".$clone."%22%20and%20SPECIES_CODE%20%3D%20%22HOMO%20SAPIENS%22%20and%20TRACE_END%20%3D%20%22FORWARD%22&dopt=trace&size=2&dispmax=5&seeas=Show"));
                }
                my $fline =  join("\t",@fields);

                my $rhit = $traces_by_template->{$clone}->{reverse}->{hit};
                my $nrt = scalar keys %{$traces_by_template->{$clone}->{forward}->{traces}};
                push(@$rev_traces, join("\t",($rhit->{trace_id},$rhit->{evalue},$rhit->{align_len})));
                my @fields = ($clone,'R',$nft,$reverse,join(",",@$rev_traces));
                if($options{append_links}) {
                    push(@fields,("driley-lx:8080/blast/$clone\_R.html","http://www.ncbi.nlm.nih.gov/Traces/trace.cgi?&cmd=retrieve&val=TEMPLATE_ID%20%3D%20%22".$clone."%22%20and%20SPECIES_CODE%20%3D%20%22HOMO%20SAPIENS%22%20and%20TRACE_END%20%3D%20%22REVERSE%22&dopt=trace&size=2&dispmax=5&seeas=Show"));
                }
                my $rline = join("\t",@fields);
                
                my $output_line = 0;
                if($forward =~ /Homo/ && $reverse !~ /Homo/) {
                    my @outfields = ($clone,$reverse,$forward,join(",",@$rev_traces),join(",",@$fwd_traces));
                    if($options{append_links}) {
                        push(@outfields,("driley-lx:8080/blast/$clone\_R.html","driley-lx:8080/blast/$clone\_F.html"));
                    }
                    $output_line = join("\t",@outfields);
                    print "$fline\n$rline\n";
                    push(@{$lgt_by_genera->{$reverse}},$clone);
                }
                elsif($reverse =~ /Homo/ && $forward !~ /Homo/) {
                    print "$fline\n$rline\n";
                    my @outfields = ($clone,$forward,$reverse,join(",",@$fwd_traces),join(",",@$rev_traces));
                    if($options{append_links}) {
                        push(@outfields, ("driley-lx:8080/blast/$clone\_F.html","driley-lx:8080/blast/$clone\_R.html"));
                    }
                    $output_line = join("\t",@outfields);
                    push(@{$lgt_by_genera->{$forward}},$clone);
                }
                elsif($reverse !~ /Homo/ && $forward !~ /Homo/) {
                    print "$fline\n$rline\n";
                    my @outfields = ($clone,$forward,$reverse,join(",",@$fwd_traces),join(",",@$rev_traces));
                    if($options{append_links}) {
                        push(@outfields,("driley-lx:8080/blast/$clone\_R.html","driley-lx:8080/blast/$clone\_F.html"));
                    }                    
                    $output_line = join("\t",@outfields);
                    push(@{$lgt_by_genera->{$forward.":".$reverse}},$clone);
                }
                if($output_line) {
                    print OUTCLONE $output_line;
                    print OUTCLONE "\n";
                }
            }
        }
        
    }
    foreach my $genera (keys %$lgt_by_genera) {
#        print join("\t",($genera,scalar @{$lgt_by_genera->{$genera}}));
#        print "\n";
    }
}
sub _process_within_trace {
    my $file = shift;
    my $type = shift;
    open IN3, "<$file" or die "Unable to open $file\n";

    while(<IN3>) {
        my @fields = split(/\t/, $_);
        
        my $hit = {
            trace_id => $fields[0],
            accession => $fields[1],
            pid => $fields[2],
            align_len => $fields[3],
            mismatches => $fields[4],
            gaps => $fields[5],
            query_start => $fields[6],
            query_end => $fields[7],
            subj_start => $fields[8],
            subj_end => $fields[9],
            evalue => $fields[10],
            score => $fields[11],
            template_id => $fields[12],
            direction => $fields[13],
            taxon_id => $fields[14],
            scientific_name => $fields[15],
            lineage => $fields[16]
        };
        
        $hit->{scientific_name} =~ /^(\w+) /;
        my $genera = $1;

        if(!$traces_by_trace_id->{$hit->{trace_id}}) {
            $traces_by_trace_id->{$hit->{trace_id}} = {};
        }
        push(@{$traces_by_trace_id->{$hit->{trace_id}}->{$genera}}, $hit);
    }
    close IN3;

}
