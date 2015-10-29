=head1 NAME

LGTsam2lca.pm - Run bwa and generate lca's

=head1 SYNOPSIS

Need to put something useful here

=head1 DESCRIPTION

A module to run bwa

=head1 AUTHOR - David R. Riley

e-mail: driley@som.umaryland.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _
my $results = GetOptions (
    \%options,
    'file_list=s',
    'out_file=s',
    'ncbitax=s',
    'gitax=s',
    'dbhost=s',
    'taxondb=s',
    'taxoncoll=s',
    'overwrite=s',
    'idx_dir=s',
    'chunk_size:s',
    'taxon_dir:s',
    'samtools_bin:s',
    'step:s',
    'check_mates:s',
    'independent_lca',
    'se_lca',
    'help|h') || pod2usage();

=cut

package LGTsam2lca;
use strict;
use version;
use File::Basename;
# Dependencies
use GiTaxon;
use LGTBestBlast;
use LGTFinder;
use LGTbwa;
use LGTSeek;

$| = 1;

my $check_mates;
my $CHUNK_SIZE;
#my $reads_by_read_id = {};
#my $reads_by_mate_id = {};
#my $seen_mates = {};
#my $seen_reads = {};
my $chunk = [];
my $counter = 0;
my $step = 0;
my $out2;
my $samtools;
my $outf;

my %options;

=head2 new

 Title   : new
 Usage   : my $lgtsam2lca = LGTsam2lca->new({out_file => $outfile,...})
 Function: Creates a new LGTsam2lca object.
 Returns : An instance of LGTsam2lca
 Args    : A hash containing potentially several config options:

           out_file - Path to the output file.
           samtools_bin - path to the samtools binary
           gi2tax - GiTaxon object
=cut
sub new {
    my ($class, $args) = @_;

    my $step = $args->{step} ? $args->{step} : 0;
    my $self = {
        reads_by_read_id => {},
        reads_by_mate_id => {},
        seen_mates => {},
        seen_reads => {},
        out_file => $args->{out_file},
        samtools_bin => $args->{samtools_bin},
        gi2tax => $args->{gi2tax},
        step => $step,
        check_mates => $args->{check_mates},
        complete_bam => $args->{complete_bam}
    };
    bless $self;
    
    if($self->{complete_bam}) {
        $self->prime_hash($self->{complete_bam});
    }
    return $self;
}

=head2 new

 Title   : runSam2Lca
 Usage   : my $output_files = $lgtsam2lca->runSam2lca()
 Function: run sam2lca on the input sam/bam list
 Returns : An object with the output files in it
 Args    : 
=cut
sub runSam2Lca {
    my $self = shift;
    my $config = shift;

    %options = %$config;

    # If the user passes in a list of files, then we'll process them here.
    if($self->{file_list}) {

#    if(defined($self->{step})) {
        my $count = `wc -l $self->{file_list}`;
        $step = $self->{step};
        print STDERR "Working on ".($step*$CHUNK_SIZE)." to ".($step*$CHUNK_SIZE+$CHUNK_SIZE)."\n";
        open IN, "<$self->{file_list}" or die "Couldn't open $self->{file_list}\n";
        my $cnter = 0;
        while(<IN>) {
            chomp;
            $self->process_file({file => $_});
        }
        print STDERR "\nFinished looping through step $step. Now writing output\n";
        my $out_obj = $self->writeOutput();
#    }
    }
    else {
        print STDERR "Called runSaml2Lca with no file list, nothing to do!\n";
    }
}

=head2 writeOutput

 Title   : writeOutput
 Usage   : my $output_files = $lgtsam2lca->writeOutput()
 Function: Write the output of the currently compiled LCAs
 Returns : An object with the output files in it.
 Args    : 
=cut
sub writeOutput {
    my $self = shift;
    my $out2;
    my $outf = $self->{out_file} ? ">$self->{out_file}": '>-';

#    if($self->{independent_lca}){
        if($self->{out_file}){
           $outf=~/^>(.+?)\.(\w+?)$/;
           $out2="$1\_independent_lca.$2";
        } else {
            $out2=">-";
        }
        open(OUT2,">","$out2") or die "Couldn't open $out2\n";
#    }

    open my $out, $outf or die "Couldn't open output\n";
    foreach my $key (keys %{$self->{reads_by_mate_id}}) {
            #print STDERR "Mate seen:$key processing with ...";
#            if($self->{independent_lca}){
             #   print STDERR "\tindependent_lca ...";
                print OUT2 join("\t",("$key\_1", $self->{reads_by_read_id}->{"$key\_1"}));
                print OUT2 "\n";
                print OUT2 join("\t",("$key\_2", $self->{reads_by_read_id}->{"$key\_2"}));
                print OUT2 "\n";
#            }
            #if($self->{se_lca}){
            if(!defined($self->{reads_by_mate_id}->{$key})) {
              #  print STDERR "\tFound no hits for $key";
                $self->{reads_by_mate_id}->{$key} = '';
            }
#            print STDERR "\tSingleEnd_lca ...";
            my $new_conservative_se_lca = &find_lca([$self->{reads_by_read_id}->{"$key\_1"},$self->{reads_by_read_id}->{"$key\_2"}]);
            if($self->{reads_by_read_id}->{"$key\_1"} =~ $self->{reads_by_read_id}->{"$key\_2"} || $self->{reads_by_read_id}->{"$key\_2"} =~ $self->{reads_by_read_id}->{"$key\_1"}){
                ## Liberal
                if(length($self->{reads_by_read_id}->{"$key\_1"}) >= length($self->{reads_by_read_id}->{"$key\_2"})){
                    print $out join("\t",($key, $self->{reads_by_mate_id}->{$key}, $new_conservative_se_lca,$self->{reads_by_read_id}->{"$key\_1"}));
                    print $out "\n";
                } else {
                    print $out join("\t",($key, $self->{reads_by_mate_id}->{$key}, $new_conservative_se_lca, $self->{reads_by_read_id}->{"$key\_2"}));
                    print $out "\n";
                }
            } else {
                print $out join("\t",($key, $self->{reads_by_mate_id}->{$key}, $new_conservative_se_lca, $new_conservative_se_lca));
                print $out "\n";
            }
 #       } else {
 #           if($self->{reads_by_mate_id}->{$key}) {
 #               print STDERR "\ttraditional lca ...";
 #               print $out "$key\t$self->{reads_by_mate_id}->{$key}\n";
#            } else {
#                print STDERR "\tFound no hits for $key";
#                print $out join("\t",($key, ""));
#                print $out "\n";
#            }
        #print STDERR "\n";
    }
    
    return {
        independent => $out2,
        normal => $outf
    };
}

sub prime_hash {
    my $self = shift;
    my $f = shift;
    open(my $in, "-|", "samtools view $f") or die "Unable to open $f for priming the hash\n";
    while(<$in>) {
        my @fields = split(/\t/);
        $self->{reads_by_mate_id}->{$fields[0]} = undef;
   }
}

sub jump_to_line {
    my $handle = shift;
    my $count = 0;
    my $start = $step*$CHUNK_SIZE;
    while ($count < $start) {
        my $line = <$handle>;
        if ($line !~ /^@/) {
            $count++;
        }
    }
    return $count;
}

sub process_sam_line {
    my($self,$l) = @_;
    chomp $l;
    # Don't count @seq lines
    if($l =~ /^@/) {
        next;
    }
    #$count++;
    my @fields = split(/\t/,$l);
    my $flag = &parseFlag($fields[1]);
    my $read_name = "$fields[0]\_1";
    if(!$flag->{'first'}) {
        $read_name = "$fields[0]\_2";
    }
    $self->{seen_mates}->{$fields[0]} = 1;
    if($flag->{query_mapped}) {
        
        my $tax = $self->{gi2tax}->getTaxon($fields[2]);
        # Here we'll deal with keeping track of things by read
        if($self->{reads_by_read_id}->{$read_name}) {
            my $lca = &find_lca([$self->{reads_by_read_id}->{$read_name},$tax->{lineage}]);
            $self->{reads_by_read_id}->{$read_name} = $lca;
            # print "$read_name\t$tax->{lineage}\t$lca\n";
            
        }
        else {
            $self->{reads_by_read_id}->{$read_name} = $tax->{lineage};
        }
        
        if(!$self->{check_mates} || $self->{check_mates} && $flag->{mate_mapped}) {
            # Here we'll keep track of things by mate
            if($self->{reads_by_mate_id}->{$fields[0]}) {
                # print STDERR "$fields[0]\n";
                my $lca = &find_lca([$self->{reads_by_mate_id}->{$fields[0]},$tax->{lineage}]);
                $self->{reads_by_mate_id}->{$fields[0]} = $lca;
            }
            else {
                # print STDERR "$fields[0]\n";
                $self->{reads_by_mate_id}->{$fields[0]} = $tax->{lineage};
            }
        }
    }
}

sub process_file {
    my ($self, $config) = @_;
    my $line = $config->{file};
    $self->{seen_mates} = {};
    $self->{seen_reads} = {};
    my $handle;
    my $presplit = 0;
    my $file = $line;
        
#    my($name,$path,$suff) = fileparse($line,('.sam','.bam'));


#    print STDERR "Checking for $path/split/$name\_$step$suff\n"; 
#    if( -e "$path/split/$name\_$step$suff") { 
#        print STDERR "Looks like we are pre-split\n";
#        $file = "$path/split/$name\_$step$suff";
#        ($name,$path,$suff) = fileparse($line,('.sam','.bam'));
#        $presplit = 1;
#    }
    if($config->{handle}) {
        $handle = $config->{handle};
    }
    elsif($file =~ /.bam$/) {
        open($handle, "-|", "$samtools view $file") or die "Unable to open $file\n";
    }
    else {
        open $handle, "<$file" or die "Unable to open $file\n";
    }

    # Skip to the start
#    if(!$presplit) {
#        print STDERR "Jumping to line $step\n";
#        &jump_to_line($handle,$step);
#    }
    # Loop till we're done.
    my $end = $CHUNK_SIZE;
    my $count = 0;
    my $l;
    my $hit = 0;

    # If we have presplit the files we'll process the whole thing
#    if($presplit) {
        while ($l = <$handle>) {
            chomp $l;
            # Don't count @seq lines
            if($l =~ /^@/) {
                next;
            }
            $count++;
            $self->process_sam_line($l);
        }
#    }

    # If we have not presplit the files we'll only process a chumk
#    else {
#        while ($count < $end && ($l = <$handle>)) {
#            chomp $l;
#            # Don't count @seq lines
#            if($l =~ /^@/) {
#                next;
#            }
#            $count++;
#            &process_sam_line($l);
#        }
#    }
#    print STDERR scalar keys (%{$self->{seen_mates}}) . " reads seen in $file\n";
#    print STDERR scalar keys (%{$self->{reads_by_mate_id}}) . " reads with hits\n";
}


sub find_lca {
    my $lineages = shift;

    # prime it
    my @lca = split(';', $lineages->[0]);

    foreach my $l (@$lineages) {
        my $newlca = [];
        my @lineage = split(';',$l);
        for( my $i = 0; $i < @lineage; $i++) {
            if($lca[$i] eq $lineage[$i]) {
                push(@$newlca, $lineage[$i]);
            }
            else {
                last;
            }   
        }
        @lca = @$newlca;
    }
    return join(';',@lca);
}

sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    $str =~ s/^0+(?=\d)//;   # otherwise you'll get leading zeros
    return $str;
}

sub parseFlag {
    my $flag = shift;
    my $rawbin = dec2bin($flag);
    my $rev = scalar $rawbin;
    if($rev eq $rawbin) {
        #    print "ERROR $rev $rawbin\n";
    }
    my $bin = sprintf("%011d", $rev);
    my $final_bin = reverse $bin;
    my $prop = substr($final_bin, 1, 1);
    my $qmap = substr($final_bin, 2, 1);
    my $mmap = substr($final_bin, 3, 1);
    my $qstrand = substr($final_bin, 4, 1);
    my $mstrand = substr($final_bin, 5, 1);
    my $first = substr($final_bin, 6, 1);
    my $last = substr($final_bin, 7, 1);
    return {
        'query_mapped' => !$qmap,
        'mate_mapped' => !$mmap,
        'first' => $first,
        'last' => $last
    };
}
1;