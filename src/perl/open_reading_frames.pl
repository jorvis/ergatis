#!/usr/bin/perl

=head1 NAME

open_reading_frames.pl - orf finding script for GOS sequence data

=head1 SYNOPSIS

USAGE: open_reading_frames.pl 
            --input_file=/path/to/some_file.fsa || --stdin
            --translation_table=11
          [ --output_dir=/some/dir
            --gzip_output=1
            --beginning_as_start=1
            --end_as_stop=1
            --assume_stops=0 
            --full_orfs=0
            --min_orf_size=180
            --max_orf_size=9999
            --min_unmasked_size=150
            --frames=1,2,3,4,5,6
            --force_methionine=0
            --header_additions='foo=bar,ergatis_id=12345'
            --unknown_aa=X
          ]

=head1 OPTIONS

B<--input_file,-i>
    Provide full path to fasta format sequence file (may contain multiple records).
    Alternatively, input can be accepted from STDIN using --stdin flag.

B<--stdin>
    Optional. Instead of providing --input_file, fasta file input can be read from

B<--output_dir>
    Optional. Directory to write output sequence files to. (default=.)

B<--gzip_output>
    Optional. Compress output sequence files using gzip.
    
B<--translation_table,-t>
    Translation table to use for ORF finding and amino acid translations.
    Provide an integer value specifying a script-supported translation table from
    the following set: (0-6, 9-16, 21-23), or provide the full path to an EMBOSS
    formatted translation table file.
    STDIN using this flag.

B<--frames>
    Optional. Translation frames for which to output ORFs. Values can be:
    'F' (forward frames, 1-3), 'R' (reverse frames, 4-6), a comma-delimited string 
    composed of integer values 1-6 such as '1,4,5,6', or '0' for all six frames (default).
    
B<--beginning_as_start>
    Optional. Treat the beginning of the sequence as a start codon and mark ORFs as partial.
    (default = 1)

B<--force_methionine>
    Optional. Force translation of first codon in any ORF as 'M'.
    
B<--full_orfs>
    Optional. Predicts ORFs from stop codon to stop codon, without requiring start codons.
    
B<--end_as_stop>
    Optional. Treat the end of the sequence as a stop codon and mark ORFs as partial.
    (default = 1)

B<--assume_stops>
    Optional. Any codon that could be translated as a stop should be. 
    Eg: TAN could be a stop codon (TAG) or tyrosine (Y) but should be marked as a stop codon.
    (default = 1)
   
B<--min_orf_size>
    Optional. Set the minimum size of ORFs to output. (default = 180)
 
B<--max_orf_size>
    Optional. Set the maximum size of ORFs to output. (default = 999999999)

B<--min_unmasked_size>
    Optional. The minimum number of non-softmasked bases and ORF must contain to be valid. (default = 150)
    
B<--unknown_aa>
    Optional. Character to use for aa translation of ambiguous or partial non-degenerate codons. (default = 'X')
    
B<--debug> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

This script is a re-implementation of Doug Rusch's GOS ORF finding script.

=head1  INPUT

Input should be cleaned na sequence reads.

=head1  OUTPUT

Output will be na and aa sequence of predicted ORFs. Two output files will be created:
    INPUT_FILE_NAME_BASE.fna[.gz]
    INPUT_FILE_NAME_BASE.faa[.gz]

=head1  CONTACT

    Brett Whitty
    bwhitty@tigr.org

=cut

use lib '/usr/local/annotation/CAMERA/lib';
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Basename qw(basename fileparse);
use GUID;

my %options = ();
my $results = GetOptions (\%options, 
                          'input_file|i=s',
                          'output_dir|o=s',
                          'gzip_output=i',
                          'stdin',
                          'frames=s',
                          'translation_table|t=s',
                          'force_methionine=i',
                          'beginning_as_start=i',
                          'full_orfs=i',
                          'end_as_stop=i',
                          'assume_stops=i',
                          'min_orf_size=i',
                          'max_orf_size=i',
                          'min_unmasked_size=i',
                          'header_additions=s',
                          'unknown_aa=s',
                          'log|l=s',
                          'debug=i',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "can't create log file: $!";
}

## id generator for globally unique orf and peptide ids
my $idgen = new GUID( 100 );

## initialize the hash that stores the codon table
## $codon_table_ref->{start|stop|aa}->{$codon}
## hash values are AA character or '*'
my $codon_table_ref = initialize_codon_table($options{'translation_table'});

## initialize array of frames to process
my $frames_ref = parse_frame_string($options{'frames'});

## initialize the hash of fasta header additional attributes 
my $header_additions_ref = {};
if ($options{'header_additions'}) {
    $header_additions_ref = parse_header_additions($options{'header_additions'});
}

my $outfh_na;
my $outfh_aa;
if ($options{'input_file'}) {
    my $output_basename = $options{'output_dir'}."/".fileparse(basename($options{'input_file'}), qr/\.[^.]*/);
    $outfh_na = open_file_write("$output_basename.fna", $options{'gzip_output'});
    $outfh_aa = open_file_write("$output_basename.faa", $options{'gzip_output'});
} elsif ($options{'stdin'}) {
    open ($outfh_na, ">&STDOUT") || die "Couldn't open STDOUT for writing: $!";
    open ($outfh_aa, ">&STDOUT") || die "Couldn't open STDOUT for writing: $!";
}

my $seq_fh;
if ($options{'stdin'}) {
    open ($seq_fh, "<&STDIN") || die "Couldn't open STDIN for reading: $!";
} else {
    $seq_fh = open_file_read($options{'input_file'});
}
my ($header, $sequence);
while (<$seq_fh>) {
    chomp;
    if (/^#/) { 
        next;
    } elsif (/^>(.*)/) {
        if ($sequence) {
            my $header_hash_ref = process_header($header);
            process_sequence(\$sequence, $header_hash_ref);
            $header = '';
            $sequence = '';
        }
        $header = $1;
    } else {
        $sequence .= $_;
    }
}
if ($sequence) {
    my $header_hash_ref = process_header($header);
    process_sequence(\$sequence, $header_hash_ref);
}


## processes an individual fasta sequence record
sub process_sequence {
    my ($sequence_ref, $header_ref) = @_;

    ## clean whitespace from sequence
    ${$sequence_ref} =~ s/\s+//g;
    
    ## create reverse complement sequence for frames 4,5,6
    my $sequence_revcomp = reverse_complement_dna(${$sequence_ref}); 

    ## store references to the forward and reverse sequences
    ## so we can iterate over them easier
    my $sequence_array_ref = [
                                $sequence_ref, 
                                $sequence_ref, 
                                $sequence_ref, 
                                \$sequence_revcomp,
                                \$sequence_revcomp,
                                \$sequence_revcomp,
                             ];
    
    my $seqlen = length(${$sequence_ref});

    my $orfs = [];
    foreach my $frame (@{$frames_ref}) {
        $orfs->[$frame] = [];
    }
    my $orf_buffer = [];
    my $orf_flags = [ 0, 0, 0, 0, 0, 0];
     
    if ($options{'beginning_as_start'}) {
        ## initialize values for option --beginning_as_start
        $orf_flags = [ 1, 1, 1, 1, 1, 1 ];
        foreach my $frame (@{$frames_ref}) {
            $orf_buffer->[$frame]->{'begin'} = 0 + $frame % 3;
            $orf_buffer->[$frame]->{'5_prime_stop'} = 0;
        }
    }
    
    ## codon startpos
    for (my $i = 0; $i < $seqlen - 2; $i += 3) {
        
        ## predict orfs in specified frames 
        foreach my $frame (@{$frames_ref}) {
            
            ## codon start position
            my $startpos = $i + $frame % 3;

            ## extract a codon
            my $codon = substr(${$sequence_array_ref->[$frame]}, $startpos, 3);
           
            ## translate the codon
            my $translated_codon = translate_codon($codon);
           
            ## if we're not in an ORF, look for start codons
            if (! $orf_flags->[$frame] && is_start_codon($codon)) {
                $orf_flags->[$frame] = 1;
                $orf_buffer->[$frame]->{'begin'} = $startpos;
            }
           
            ## if we're in an ORF
            if ($orf_flags->[$frame]) {
                
                ## store codon in orf buffer
                push(@{$orf_buffer->[$frame]->{'codon'}}, $codon);
            
                ## store aa in orf buffer
                push(@{$orf_buffer->[$frame]->{'aa'}}, $translated_codon);
            
                ## if we've hit a stop codon, store the orf details
                if (is_stop_codon($codon)) {
                    
                    ## reset orf flag (end of orf)
                    $orf_flags->[$frame] = 0;

                    ## save some typing
                    my $orf_ref = $orf_buffer->[$frame];
                    
                    ## attach parent header ref
                    $orf_ref->{'_header'} = $header_ref;
                    ## attach header additions ref
                    $orf_ref->{'_header_additions'} = $header_additions_ref;
                    
                    ## set orf translation table attribute
                    $orf_ref->{'ttable'} = get_translation_table_string();
                    ## set three_prime_stop to this stop codon
                    $orf_ref->{'3_prime_stop'} = $codon;
                    ## orf end is current codon startpos + 3
                    $orf_ref->{'end'} = $startpos + 3;
                    ## set orientation attribute
                    $orf_ref->{'orientation'} = ($frame < 3) ? 1 : -1;
                    $orf_ref->{'length'} = $orf_ref->{'end'} - $orf_ref->{'begin'};
                    ## adjust begin and end coordinates if on reverse strand
                    if ($frame > 2) {
                        my $begin = $seqlen - $orf_ref->{'end'};
                        my $end   = $seqlen - $orf_ref->{'begin'};
                        $orf_ref->{'begin'} = $begin;
                        $orf_ref->{'end'}   = $end;
                    }
  
                    process_orf($orfs->[$frame], $orf_ref);
                    
                    ## reset the ORF buffer for this frame
                    $orf_buffer->[$frame] = {};
                    ## mark current stop as next ORF's five_prime_stop
                    $orf_buffer->[$frame]->{'5_prime_stop'} = $codon;
                }
            }
        }
    }
    ## if we reached the end of the sequence without hitting a stop codon
    ## then any partial orfs will be processed if --end_as_stop is enabled
    if ($options{'end_as_stop'}) {
        foreach my $frame (@{$frames_ref}) {
            ## save some typing
            my $orf_ref = $orf_buffer->[$frame];
            if (defined($orf_ref->{'codon'})) {

                ## no 3' stop
                $orf_ref->{'3_prime_stop'} = '0';
                
                ## attach parent header ref
                $orf_ref->{'_header'} = $header_ref;
                ## attach header additions ref
                $orf_ref->{'_header_additions'} = $header_additions_ref;
                    
                ## set orf translation table attribute
                $orf_ref->{'ttable'} = get_translation_table_string();
                $orf_ref->{'end'} = $seqlen;
                $orf_ref->{'orientation'} = ($frame < 3) ? 1 : -1;
                $orf_ref->{'length'} = $orf_ref->{'end'} - $orf_ref->{'begin'};
                ## adjust begin and end coordinates if on reverse strand
                if ($frame > 2) {
                    my $begin = $seqlen - $orf_ref->{'end'};
                    my $end   = $seqlen - $orf_ref->{'begin'};
                    $orf_ref->{'begin'} = $begin;
                    $orf_ref->{'end'}   = $end;
                }
 
                process_orf($orfs->[$frame], $orf_buffer->[$frame]);
                
            }
        } 
    }

    ## commenting this out here and in process_orf to save memory, to enable -> SEARCH_ON_THIS
    # return $orfs;
    
}

## does whatever we're going to do with an orf once we've got one
sub process_orf {
    my ($orf_array_ref, $orf_ref) = @_;

    ## enforce limits on softmasked regions
    if ($options{'min_unmasked_size'}) {
        my $unmasked_size = 0;
        foreach my $codon(@{$orf_ref->{'codon'}}) {
            $codon =~ s/[a-z]+//g;
            $unmasked_size += length($codon);
        }
        if ($unmasked_size < $options{'min_unmasked_size'}) {
            return;
        } 
    }
  
    ## enforce orf size limits
    if ($orf_ref->{'length'} <= $options{'max_orf_size'} && $orf_ref->{'length'} >= $options{'min_orf_size'}) {
        
        ## force first codon to translate as methionine
        if ($options{'force_methionine'}) {
            $orf_ref->{'aa'}->[0] = 'M';
        }
       
        ## create fasta headers for the orf
        create_fasta_headers($orf_ref);
      
        ## write the orf to fasta output
        write_orf_to_fasta($orf_ref);
        
        ## not currently storing orfs to save memory -> SEARCH_ON_THIS
        # push (@{$orf_array_ref}, $orf_ref);
        
    }

}

## outputs an orf hash ref to na and aa fasta files
sub write_orf_to_fasta {
    my ($orf) = @_;

    ## write na sequence
    write_seq_to_fasta($outfh_na, $orf->{'orf_header'}, join("",@{$orf->{'codon'}}));
    ## write aa sequence
    write_seq_to_fasta($outfh_aa, $orf->{'pep_header'}, join("",@{$orf->{'aa'}}));
}

## writes a record to a fasta file
sub write_seq_to_fasta {
    my ($fh, $header, $sequence) = @_;

    $sequence =~ s/(.{1,60})/$1\n/g;

    print $fh ">$header\n$sequence";
}


## creates strings for the orf and peptide fasta headers and stores them in the orf hash
sub create_fasta_headers {
    my ($orf_ref) = @_;
   
    ## create a temporary duplicate of the orf hash, minus codons and aas
    my %orf = %{$orf_ref};
    delete @orf{'codon', 'aa'};
   
    ## move attributes onto orf hash from read header hash ref
    $orf{'read_id'} = $orf{'_header'}->{'read_id'};
    $orf{'read_defline'} = $orf{'_header'}->{'read_defline'};
    delete $orf{'_header'};
    
    my ($orf_header, $pep_header);
    
    ## create new orf / pep ids
    my $orf_id = "JCVI_ORF_".$idgen->getGUID();
    my $pep_id = "JCVI_PEP_".$idgen->getGUID();
   
    $orf_header = "$orf_id /pep_id=$pep_id";
    $pep_header = "$pep_id /orf_id=$orf_id";
    
    my @ordered_attributes = (
                                'read_id',
                                'begin',
                                'end',
                                'orientation',
                                '5_prime_stop',
                                '3_prime_stop',
                                'ttable',
                             );
    
    foreach my $attribute (@ordered_attributes) {
        $orf_header .= " /".$attribute."=".$orf{$attribute};
        $pep_header .= " /".$attribute."=".$orf{$attribute};
        delete $orf{$attribute};
    }
   
    ## add remaining attributes in no particular order 
    foreach my $attribute (keys(%orf)) {
        ## deal with attributes that are hash references of other attributes
        if ($attribute =~ /^_/) {
            foreach my $sub_att (keys(%{$orf{$attribute}})) {
                $orf_header .= " /".$sub_att."=".$orf{$attribute}->{$sub_att};
                $pep_header .= " /".$sub_att."=".$orf{$attribute}->{$sub_att};
            }
        ## deal with normal attributes
        } else {
                $orf_header .= " /".$attribute."=".$orf{$attribute};
                $pep_header .= " /".$attribute."=".$orf{$attribute};
        }
    }    
    
    $orf_ref->{'orf_header'} = $orf_header;
    $orf_ref->{'pep_header'} = $pep_header;
    
    return $orf_ref;
}

## parses the input sequence header and returns a hash of attributes
sub process_header {
    my ($header) = @_;
   
    my $header_hash_ref = {};
   
    $header =~ /^(\S+)\s*(.*)/ || die "Couldn't parse read ID from header:\n$header";
    
    $header_hash_ref->{'read_id'} = $1;
    $header_hash_ref->{'read_defline'} = "\"$2\"";
    
    ## return a hash of whatever values are obtained from the header
    return $header_hash_ref;
}

## uses the codon table and whatever other rules we're applying
## to translate the codon into an amino acid character
sub translate_codon {
    my ($codon) = @_;

    ## support softmasked sequence
    $codon = uc($codon);
    
    my $aa;
   
    ## append N's to partial codons 
    if (length($codon) != 3) {
        $codon = substr($codon."NNN", 0, 3);
    }
    
    if (defined($codon_table_ref->{'aa'}->{$codon})) {
        $aa = $codon_table_ref->{'aa'}->{$codon};
    } else {
        $aa = $options{'unknown_aa'};
#        die "no codon lookup for -> $codon";
    }
    
    ## return the translated aa
    return $aa;
}

sub reverse_complement_dna {
        my ($r_seq) = @_;
        $r_seq =~ tr/AaCcGgTtMmRrWwSsYyKkVvHhDdBb/TtGgCcAaKkYyWwSsRrMmBbDdHhVv/;
        $r_seq = reverse($r_seq);
        return $r_seq;
}

sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}

## read a codon table and initialize a hash structure for storing it
sub initialize_codon_table {
    my ($table) = @_;

    my $table_ref = [];
    my $aa_ref = {};
    my $start_ref = {};
    my $stop_ref = {};
    my $table_string = '';
    
    my $parsed_code;
 
    my %supported_tables;
    
    ## initialize a flag hash for the internally supported codon tables
    foreach my $table_id ((0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23)) {
        $supported_tables{$table_id} = 1;
    } 
    
    my $found_table_flag = 0;
    
    my $table_fh;
    if ($supported_tables{$table}) {
        open ($table_fh, "<&DATA") || die "Couldn't open DATA for reading";
    } else {
        no warnings;
        if (int($table) eq $table) { ## stupid that this raises a warning
        use warnings;
            die "Codon table '$table' is not supported, must provide path to an EMBOSS format codon table file";
        } else {
            open ($table_fh, $table) || die "Failed to open codon table file '$table': $!";
            $found_table_flag = 1;
            $table_string = 'custom';
        }
    }
    
    while (<$table_fh>) {
        chomp;
        if (/^#/ || /^\s+$/) {
            next;
        }
        if (/Genetic Code \[(\d+)\]/) {
            $parsed_code = $1;
            if ($parsed_code eq $table) {
                $found_table_flag = 1;
                $table_string = $table;
            }
        }
        if ($found_table_flag) {
            if (/AAs\s+=\s+(.*)/) {
                $table_ref->[0] = [split('', $1)];
            } elsif (/Starts\s+=\s+(.*)/) {
                $table_ref->[1] = [split('', $1)];
            } elsif (/Base1\s+=\s+(.*)/) {
                $table_ref->[2] = [split('', $1)];
            } elsif (/Base2\s+=\s+(.*)/) {
                $table_ref->[3] = [split('', $1)];
            } elsif (/Base3\s+=\s+(.*)/) {
                $table_ref->[4] = [split('', $1)];
                last;
            }
        }
    }
    close $table_fh;
     
    unless ($found_table_flag) {
        die "Couldn't initialize codon table '$table'";
    }
   
    ## temp hash for finding degenerate codons
    my $ncodons_ref = {};
    while (scalar(@{$table_ref->[0]})) {
        my $aa    = shift(@{$table_ref->[0]});
        my $start = shift(@{$table_ref->[1]});
        my $base1 = shift(@{$table_ref->[2]}); 
        my $base2 = shift(@{$table_ref->[3]}); 
        my $base3 = shift(@{$table_ref->[4]}); 
        my $codon = $base1 . $base2 . $base3;

        $aa_ref->{$codon} = $aa;
        
        $ncodons_ref->{$base1.$base2.'N'}->{$aa} = 1;
        $ncodons_ref->{$base1.'N'.$base3}->{$aa} = 1;
        $ncodons_ref->{'N'.$base2.$base3}->{$aa} = 1;
        $ncodons_ref->{$base1.'N'.'N'}->{$aa}    = 1;
        
        ## store if not '-' --- may want to put other conditions on here
        if ($start ne '-') {
            $start_ref->{$codon} = $start;
        }

        ## might as well store a hash of stop codons as 'stop'
        if ($aa eq '*') {
            $stop_ref->{$codon} = $aa;
        }
        
    }
   
    ## add codons for unknown bases
    foreach my $codon(keys(%{$ncodons_ref})) {
        
        ## if it could be a stop, we'll assume it is
        if ($options{'assume_stops'}) {
            if ($ncodons_ref->{$codon}->{'*'}) {
                $stop_ref->{$codon} = '*';
            }
        }
        
        ## @aas = the set of aas coded by the ambiguous codon
        my @aas = keys(%{$ncodons_ref->{$codon}});
        
        ## add codon lookups for unknowns at degenerate positions
        if (scalar(@aas) == 1) {
            $aa_ref->{$codon} = shift(@aas);
        }
    }

    ## NNN is to be treated as a hard stop
    $stop_ref->{'NNN'} = '*';
     
    return { 'aa' => $aa_ref, 'start' => $start_ref, 'stop' => $stop_ref, 'table' => $table_string};
}

## boolean conditional based on whether codon is a start codon or not
sub is_start_codon {
    my ($codon) = @_;
   
    ## support softmasked sequence
    $codon = uc($codon);
    
    ## full_orfs flag treats all non-stop-codons as starts
    if ($options{'full_orfs'}) {
        if (is_stop_codon($codon)) {
            return 0;
        } else {
            return 1;
        }
    }
    
    if ($codon_table_ref->{'start'}->{$codon}) {
        return 1;
    } else {
        return 0;
    }
}

## boolean conditional based on whether codon is a stop codon or not
sub is_stop_codon {
    my ($codon) = @_;
   
    ## support softmasked sequence
    $codon = uc($codon);
    
    if ($codon_table_ref->{'stop'}->{$codon}) {
        return 1;
    } else {
        return 0;
    }
}

## finds coordinates of masked regions in a sequence
## default is lowercase masking, but any characters can be provided to the regex
sub find_masked_regions {
    my ($seq_ref, $mask_char) = shift @_;
    
    unless (defined($mask_char)) {
        $mask_char = 'a-z';
    }
    
    my @regions = ();
    
    while (${$seq_ref} =~ /[$mask_char]{1,}/g) {
        push (@regions, [$-[0],$+[0]]);
    }
    return \@regions;
}

sub get_translation_table_string {
    return $codon_table_ref->{'table'};
}

## parses the string provided to the frames option
sub parse_frame_string {
    my ($frame_string) = @_;

    my $f_string = $frame_string;
    
    my @frames;
   
    if (! defined($frame_string)) {
        $frame_string = '';
    }
   
    $frame_string =~ s/\s+/,/g;
    $frame_string =~ s/[^0-9^F^R^,]+//g;
    
    if ($frame_string eq 'F') {
        @frames = (0, 1, 2);
    } elsif ($frame_string eq 'R') {
        @frames = (3, 4, 5);
    } elsif ($frame_string eq '' || $frame_string eq '0') {
        @frames = (0, 1, 2, 3, 4, 5);
    } else {
        foreach my $frame (split(",", $frame_string)) {
            if ($frame >= 1 && $frame <= 6) {
                $frame--;
                push (@frames, $frame);
            }
        }
    }
    
    if (! @frames) {
        die "Frames flag passed bad value, --frames '$f_string'";
    }
   
    return \@frames;
}

## parses the header_additions flag value, eg: 'attrib1=some_value1,attrib2=some_value2'
sub parse_header_additions {
    my ($header_additions) = @_;

    my $header_additions_hash_ref = {};
    
    my @attribs = split(",", $header_additions);

    foreach my $attrib (@attribs) {
        my ($key, $value) = split("=", $attrib);
        unless ($key && $value) {
            die "Couldn't parse a key and a value from '$attrib', from --header_additions '$header_additions'";
        }
        $header_additions_hash_ref->{$key} = $value;
    }

    return $header_additions_hash_ref;
}

## opens a filehandle for reading
sub open_file_read {
    my ($filename) = @_;

    my $fh;
    
    if (! -e $filename && -e $filename.'.gz') {
        $filename .= '.gz';
    }
    if ($filename =~ /\.gz$|\.gzip$/) {
        open ($fh, "<:gzip", $filename) || die "Couldn't open '$filename' for reading: $!";
    } else {
        open ($fh, "<$filename") || die "Couldn't open '$filename' for reading: $!";
    }
   
    return $fh 
}

## opens a filehandle for writing
sub open_file_write {
    my ($filename, $gzip_mode) = @_;

    my $fh;
    
    if ($gzip_mode) {
        open ($fh, ">:gzip", $filename.'.gz') || die "Couldn't open '$filename' for writing: $!";
    } else {
        open ($fh, ">$filename") || die "Couldn't open '$filename' for writing: $!";
    }
   
    return $fh 
}

## check options provided to the script and set defaults
sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( translation_table );
    for my $option ( @required ) {
        unless  ( defined $options->{$option} ) {
            die "--$option is a required option";
        }
    }
  
    unless ($options->{'input_file'} || $options->{'stdin'}) {
        die "either --stdin or --input_file must be specified";
    }
    
    ## handle some defaults
    $options->{'min_unmasked_size'}     = 150      unless (defined($options->{'min_unmasked_size'}));
    $options->{'min_orf_size'}          = 180       unless (defined($options->{'min_orf_size'}));
    $options->{'max_orf_size'}          = 999999999 if (! defined($options->{'max_orf_size'}) 
                                                        || $options->{'max_orf_size'} == 0);
    $options->{'beginning_as_start'}    = 1         unless (defined($options->{'beginning_as_start'}));
    $options->{'end_as_stop'}           = 1         unless (defined($options->{'end_as_stop'}));
    $options->{'output_dir'}            = '.'       unless (defined($options->{'output_dir'}));
    $options->{'gzip_output'}           = 1         unless (defined($options->{'gzip_output'}));
    $options->{'unknown_aa'}            = 'X'         unless (defined($options->{'unknown_aa'}));

    if (length($options->{'unknown_aa'}) < 1 || length($options->{'unknown_aa'}) > 1) {
        print STDERR "value provided to flag --unknown_aa should be one character in length, setting to default 'X'\n";
        $options->{'unknown_aa'} = 'X';
    }
}
__END__

## some codon tables are attached here
__DATA__
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#   Version 3.1 - 1995
#      Initial base data set from Andrzej Elzanowski (PIR/NCBI)
#
# Differs from Genetic Code [1] only in that the initiation sites have been
# changed to only 'AUG'.
#
# This is intended to be the default genetic code where the use of
# the rare initiation sites (CUG, UUG) is not intended.

Genetic Code [0]

Standard with AUG start only
 
AAs  =   FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#  Version 3.4
#     Added CTG,TTG as allowed alternate start codons in Standard code.
#        Prats et al. 1989, Hann et al. 1992
#
# Initiation Codon:
#
# AUG 
# 
# Alternative Initiation Codons
#
# In rare cases, translation in eukaryotes can be initiated from codons
# other than AUG.  A well documented case (including direct protein
# sequencing) is the GUG start of a ribosomal P protein of the fungus
# Candida albicans (Abramczyk et al.).  Other examples can be found in the
# following references: Peabody 1989; Prats et al.  1989; Hann et al. 
# 1992; Sugihara et al.  1990. 
#
# GUG, CUG, UUG 

Genetic Code [1]

Standard
 
AAs  =   FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ---M---------------M---------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#   Version 3.1 - 1995
#      Initial base data set from Andrzej Elzanowski (PIR/NCBI)
#
# Differences from the Standard Code: 
# 
#         Code 2          Standard
# 
#  AGA    Ter  *          Arg  R
#  AGG    Ter  *          Arg  R
#  AUA    Met  M          Ile  I
#  UGA    Trp  W          Ter  *
# 
# 
# Alternative Initiation Codon:
# 
# Bos: AUA 
# Homo: AUA, AUU
# Mus: AUA, AUU, AUC
# Coturnix, Gallus: also GUG (Desjardins and Morais, 1991)
# 
# Systematic Range:
# 
# Vertebrata
# 
# Comment: 
# 
# The transcripts of several vertebrate mitochondrial genes end in U or
# UA, which become termination codons (UAA) upon subsequent
# polyadenylation. 

Genetic Code [2]
 
Vertebrate Mitochondrial

AAs  =   FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG
Starts = --------------------------------MMMM---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#   Version 3.1 - 1995
#      Initial base data set from Andrzej Elzanowski (PIR/NCBI)
#
# Differences from the Standard Code: 
# 
#         Code 3          Standard
# 
#  AUA    Met  M          Ile  I
#  CUU    Thr  T          Leu  L
#  CUC    Thr  T          Leu  L
#  CUA    Thr  T          Leu  L
#  CUG    Thr  T          Leu  L
#  UGA    Trp  W          Ter  *
# 
#  CGA    absent          Arg  R
#  CGC    absent          Arg  R
# 
# 
# Systematic Range: 
# 
# Saccharomyces cerevisiae, Candida glabrata, Hansenula saturnus,
# and Kluyveromyces thermotolerans
# (Clark-Walker and Weiller, 1994)
# 
# Comments:
# 
# The remaining CGN codons are rare in Saccharomyces cerevisiae and
# absent in Candida glabrata (= Torulopsis glabrata). 
#
# The AUA codon is common in the gene var1 coding for the single
# mitochonLIial ribosomal protein, but rare in genes encoding the enzymes. 
#
# The coding assignments of the AUA (Met or Ile) and CUU (possibly Leu,
# not Thr) are uncertain in Hansenula saturnus. 
#
# The coding assignment of Thr to CUN is uncertain in Kluyveromyces
# thermotolerans (Clark-Walker and Weiller, 1994). 

Genetic Code [3]
 
Yeast Mitochondrial
 
AAs  =   FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ----------------------------------MM----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#   Version 3.1 - 1995
#      Initial base data set from Andrzej Elzanowski (PIR/NCBI)
#
# Differences from the Standard Code: 
# 
#         Code 4         Standard
# 
#  UGA    Trp  W          Ter  *
# 
# 
# Alternative Initiation Codons: 
# 
# Trypanosoma: UUA, UUG, CUG
# Leishmania: AUU, AUA 
# Tertrahymena: AUU, AUA, AUG 
# Paramecium: AUU, AUA, AUG, AUC, GUG, GUA(?) 
# (Pritchard et al., 1990)
# 
# Systematic Range: 
# 
# Mycoplasmatales: Mycoplasma, Spiroplasma (Bove et al., 1989); 
# 
# Fungi: Emericella nidulans, Neurospora crassa, Podospora anserina,
# Acremonium (Fox, 1987), Candida parapsilosis (Guelin et al., 1991),
# Trichophyton rubrum (de Bievre and Dujon, 1992), Dekkera/Brettanomyces,
# Eeniella (Hoeben et al., 1993), and probably Ascobolus immersus,
# Aspergillus amstelodami, Claviceps purpurea, and Cochliobolus
# heterostrophus. 
# 
# Protozoa: Trypanosoma brucei, Leishmania tarentolae, Paramecium
# tetraurelia, Tetrahymena pyriformis and probably Plasmodium gallinaceum
# (Aldritt et al., 1989)]. 
# 
# Metazoa: Coelenterata (Ctenophora and Cnidaria) 
# 
# Comments: 
# 
# This code is also used for the kinetoplast DNA (maxicircles,
# minicircles).  Kinetoplasts are modified mitochondria (or their parts). 
# 
# This code is not used in the Acholeplasmataceae and plant-pathogenic
# mycoplasma-like organisms (MLO) (Lim and Sears, 1992)

Genetic Code [4]
  
Mold, Protozoan, Coelenterate Mitochondrial and Mycoplasma/Spiroplasma
  
AAs  =   FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = --MM---------------M------------MMMM---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#  Version 3.3 - 10/13/95
#     Added alternate intiation codon ATC to code 5
#        based on complete mitochondrial genome of honeybee
#        Crozier and Crozier (1993)
#
#  Version 3.2 - 6/24/95
#     GTG allowed as alternate initiator
#
# Comment: 
# 
# The codon AGG is absent in Drosophila. 
# 
# Differences from the Standard Code: 
# 
# 
#         Code 5          Standard
# 
#  AGA    Ser  S          Arg  R
#  AGG    Ser  S          Arg  R
#  AUA    Met  M          Ile  I
#  UGA    Trp  W          Ter  *
# 
# 
# Alternative Initiation Codons: 
# 
# AUA, AUU
# AUC: Apis (Crozier and Crozier, 1993)
# GUG: Polyplacophora (Boore and Brown, 1994 GenBank Accession Number: U09810)
# UUG: Ascaris, Caenorhabditis
# 
# Systematic Range:
# 
# Nematoda: Ascaris, Caenorhabditis;
# Mollusca: Bivalvia (Hoffmann et al., 1992); Polyplacophora (Boore and
# Brown, 1994)
# Arthropoda/Crustacea: Artemia (Batuecas et al., 1988);
# Arthropoda/Insecta: Drosophila [Locusta migratoria (migratory locust),
# Apis mellifera (honeybee)]
# 
# Comments: 
# 
# GUG may possibly function as an initiator in Drosophila (Clary and
# Wolstenholme, 1985; Gadaleta et al., 1988).  AUU is not used as an
# initiator in Mytilus (Hoffmann et al., 1992). 
# 
# "An exceptional mechanism must operate for initiation of translation of
# the cytochrome oxidase subunit I mRNA in both D.  melanogaster (de
# Bruijn, 1983) and D.  yakuba (Clary and Wolstenholme 1983), since its
# only plausible initiation codon, AUA, is out of frame with the rest of
# the gene.  Initiation appears to require the "reading" of of an AUAA
# quadruplet, which would be equivalent to initiation at AUA followed
# immediately by a specific ribosomal frameshift.  Another possible
# mechanism ...  is that the mRNA is "edited" to bring the AUA initiation
# into frame." (Fox, 1987)

Genetic Code [5]
    
Invertebrate Mitochondrial
  
AAs  =   FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG
Starts = ---M----------------------------MMMM---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG


# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#   Version 3.1 - 1995
#      Initial base data set from Andrzej Elzanowski (PIR/NCBI)
#
# Differences from the Standard Code: 
# 
#           Code 6       Standard
# 
#  UAA      Gln  Q        Ter  *
#  UAG      Gln  Q        Ter  *
# 
# 
# Systematic Range: 
# 
# Ciliata: Oxytricha and Stylonychia (Hoffman et al.  1995), Paramecium,
# Tetrahymena, Oxytrichidae and probably Glaucoma chattoni. 
# 
# Dasycladaceae: Acetabularia (Schneider et al., 1989) and Batophora
# (Schneider and de Groot, 1991). 
# 
# Diplomonadida: 
# Scope: Hexamita inflata, Diplomonadida ATCC50330, and ATCC50380. 
# Ref.: Keeling, P.J.  and Doolittle, W.F.  1996.  A non-canonical genetic
# code in an early diverging eukaryotic lineage.  The EMBO Journal 15,
# 2285-2290. 
# 
# Comment: 
# 
# The ciliate macronuclear code has not been determined completely.  The
# codon UAA is known to code for Gln only in the Oxytrichidae. 

Genetic Code [6]
 
Ciliate, Dasycladacean and Hexamita Nuclear
    
AAs  =   FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
 
# Genetic Code Table 
# 
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
# 
#  Version 3.8
#     Added GTG start to Echinoderm mitochondrial code, code 9
#
# Differences from the Standard Code: 
# 
# 
#           Code 9        Standard
# 
#  AAA      Asn  N        Lys K
#  AGA      Ser  S        Arg R
#  AGG      Ser  S        Arg R
#  UGA      Trp  W        Ter *
# 
# 
# Systematic Range: 
# 
# Asterozoa (starfishes) (Himeno et al., 1987) Echinozoa (sea urchins)
# (Jacobs et al., 1988; Cantatore et al., 1989)

Genetic Code [9]
   
Echinoderm and Flatworm Mitochondrial
       
AAs  =   FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG
Starts = -----------------------------------M---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#  Version 3.2 - 6/24/95
#     Alternative Ciliate Macronuclear renamed to Euplotid Macronuclear
#
# Differences from the Standard Code:
# 	Code 10     Standard
# UGA 	Cys  C        Ter  *
# 
# Systematic Range: 
# Ciliata: Euplotidae (Hoffman et al. 1995). 

Genetic Code [10]
   
Euplotid Nuclear
    
AAs  =   FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
  
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#  Version 3.1 - 1995
#     Addition of Eubacterial by J.Ostell at NCBI
#  Version 3.2 - 6/24/95
#     Eubacterial renamed to Bacterial as most alternate starts
#                 have been found in Archaea
#
# Differences from the Standard Code: 
# 
# None 
# 
# Alternative Initiation Codons: 
# 
# GUG, UUG, AUU, CUG 
# 
# Systematic Range and Comments: 
# 
# Table 11 is used for Bacteria, Archaea, prokaryotic viruses and
# chloroplast proteins.  As in the standard code, initiation is most
# efficient at AUG.  In addition, GUG and UUG starts are documented in
# Archaea and Bacteria (Kozak 1983, Fotheringham et al.  1986, Golderer et
# al.  1995, Nolling et al.  1995, Sazuka & Ohara 1996, Genser et al. 
# 1998, Wang et al.  2003).  In E.  coli, UUG is estimated to serve as
# initiator for about 3% of the bacterium's proteins (Blattner et al. 
# 1997).  CUG is known to function as an initiator for one plasmid-encoded
# protein (RepA) in Escherichia coli (Spiers and Bergquist, 1992).  In
# addition to the NUG initiations, in rare cases Bacteria can initiate
# translation from an AUU codon as e.g.  in the case of poly(A) polymerase
# PcnB and the InfC gene that codes for translation initiation factor IF3
# (Polard et al.  1991, Liveris et al.  1993, Sazuka & Ohara 1996, Binns &
# Masters 2002).  The internal assignments are the same as in the standard
# code though UGA codes at low efficiency for Trp in Bacillus subtilis
# and, presumably, in Escherichia coli (Hatfiled and Diamond, 1993). 

Genetic Code [11]

Bacterial and Plant Plastid
 
AAs  =   FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ---M---------------M------------MMMM---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#   Version 3.1 - 1995
#      Initial base data set from Andrzej Elzanowski (PIR/NCBI)
#      Addition of Alternative Yeast by J.Ostell at NCBI
#
# Differences from the Standard Code: 
# 
#            Code 12      Standard
# 
#  CUG       Ser          Leu
#        
# 
# Alternative Initiation Codons: 
# 
# CAG may be used in Candida albicans (Santos et al., 1993). 
# 
# Systematic Range: 
# 
# Endomycetales (yeasts): Candida albicans, Candida cylindracea, Candida
# melibiosica, Candida parapsilosis, and Candida rugosa (Ohama et al.,
# 1993). 
# 
# Comment: 
# 
# However, other yeast, including Saccharomyces cerevisiae, Candida azyma,
# Candida diversa, Candida magnoliae, Candida rugopelliculosa, Yarrowia
# lipolytica, and Zygoascus hellenicus, definitely use the standard
# (nuclear) code (Ohama et al., 1993). 


Genetic Code [12]
  
Alternative Yeast Nuclear

AAs  =   FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -------------------M---------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#   Version 3.1 - 1995
#      Initial base data set from Andrzej Elzanowski (PIR/NCBI)
#
# Differences from the Standard Code: 
# 
#           Code 13     Standard
# 
#  AGA      Gly  G        Arg  R
#  AGG      Gly  G        Arg  R
#  AUA      Met  M        Ile  I
#  UGA      Trp  W        Ter  *
# 
# 
# Systematic Range and Comments: 
# 
# There is evidence from a phylogenetically diverse sample of tunicates
# (Urochordata) that AGA and AGG code for glycine.  In other organisms,
# AGA/AGG code for either arginine or serine and in vertebrate
# mitochondria they code a STOP.  Evidence for glycine translation of
# AGA/AGG has been found in Pyura stolonifera (Durrheim et al.  1993),
# Halocynthia roretzi (Kondow et al.  1999, Yokobori et al., 1993,
# Yokobori et al.  1999) and Ciona savignyi (Yokobori et al.  2003).  In
# addition, the Halocynthia roretzi mitochondrial genome encodes an
# additional tRNA gene with the anticodon U*CU that is thought to enable
# the use of AGA or AGG codons for glycine and the gene has been shown to
# be transcribed in vivo (Kondow et al.  1999, Yokobori et al.  1999). 


Genetic Code [13]

Ascidian Mitochondrial
  
AAs  =   FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG


# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#   Version 3.1 - 1995
#      Initial base data set from Andrzej Elzanowski (PIR/NCBI)
#
# Differences from the Standard Code: 
# 
#           Code 14      Standard
# 
#  AAA      Asn  N       Lys  K
#  AGA      Ser  S       Arg  R
#  AGG      Ser  S       Arg  R
#  UAA      Tyr  Y       Ter  *
#  UGA      Trp  W       Ter  *
# 
# 
# Systematic Range: 
# 
# Platyhelminthes (flatworms) 
# 
# Comments:
#
# Code 14 differs from code 9 only by translating UAA to Tyr rather than
# STOP.  A recent study [PMID:11027335] has found no evidence that the
# codon UAA codes for Tyr in the flatworms but other opinions exist. 
# There are very few GenBank records that are translated with code 14 but
# a test translation shows that retranslating these records with code 9
# can cause premature terminations.  Therefore, GenBank will maintain code
# 14 until further information become available. 

Genetic Code [14]

Alternative Flatworm Mitochondrial
  
AAs  =   FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG


# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
# Version 3.2 - 6/24/95
#    Blepharisma Macronuclear code added
#
# Differences from the Standard Code: 
# 
#           Code 10       Standard
# 
# UAG       Gln  Q        Ter  *
# 
# 
# Systematic Range: 
# 
# Ciliata: Blepharisma (Liang and Heckman, 1993) 


Genetic Code [15]
  
Blepharisma Nuclear
  
AAs  =   FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
     
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#  Version 3.5
#     Added code 16, Chlorophycean Mitochondrial
#       (TAG can translated to Leucine instaed to STOP in chlorophyceans
#        and fungi)
#
# Systematic Range:
#
# Chlorophyceae: Hayashi-Ishiimaru, Y, T.  Ohama, Y.  Kawatsu, K. 
# Nakamura, S.  Osawa, 1996.  Current Genetics 30: 29-33

Genetic Code [16]

Chlorophycean Mitochondrial

AAs  =   FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#  Version 3.5
#     Added code 21, Trematode Mitochondrial
#       (as deduced from: Garey & Wolstenholme,1989; Ohama et al, 1990)
#
# Systematic Range:
#
# Trematoda: Ohama, T, S.  Osawa, K.  Watanabe, T.H.  Jukes, 1990.  J. 
# Molec Evol.  30 Garey, J.R.  and D.R.  Wolstenholme, 1989.  J.  Molec. 
# Evol.  28: 374-387 329-332. 
# 
# Other Alternative Initiation Codons
#
# GUG, UUG (and possibly CUG) in the Archaea (Noelling et al., unpublished) 
#
# AUA, GUG, UUG, and AUC or AAG may be used (at least in experimental
# systems) by the yeasts Saccharomyces cerevisiae (Olsen, 1987, and references
# therein). 
#
# ACG initiates translation of certain proteins in the adeno-associated virus
# type 2 (Becerra et al., 1985), the phage T7 mutant CR17 (Anderson and
# Buzash-Pollert, 1985), Sendai virus (Gupta and Patwardhan, 1988), and rice
# chloroplast (Hiratsuka et al., 1989). Also, it is the most effective non-AUG
# initiation codon in mammalin cells (Koepke and Leggatt, 1991). 
#
# CUG is the initiation codon for one of the two alternative products of the
# human c-myc gene (Hann et al., 1987). 
#

Genetic Code [21]

Trematode Mitochondrial

AAs  =   FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG
Starts = -----------------------------------M---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#  Version 3.6
#     Added code 22 TAG-Leu, TCA-stop
#        found in mitochondrial DNA of Scenedesmus obliquus
#        submitted by Gertraude Berger, Ph.D.
#        Organelle Genome Megasequencing Program, Univ Montreal
#
# Other Alternative Initiation Codons
#
# GUG, UUG (and possibly CUG) in the Archaea (Noelling et al., unpublished) 
# AUA, GUG, UUG, and AUC or AAG may be used (at least in experimental
# systems) by the yeasts Saccharomyces cerevisiae (Olsen, 1987, and references
# therein). 
#
# ACG initiates translation of certain proteins in the adeno-associated virus
# type 2 (Becerra et al., 1985), the phage T7 mutant CR17 (Anderson and
# Buzash-Pollert, 1985), Sendai virus (Gupta and Patwardhan, 1988), and rice
# chloroplast (Hiratsuka et al., 1989). Also, it is the most effective non-AUG
# initiation codon in mammalin cells (Koepke and Leggatt, 1991). 
#
# CUG is the initiation codon for one of the two alternative products of the
# human c-myc gene (Hann et al., 1987). 
#
# Systematic Range:
#
# Scenedesmus obliquus: Nedelcu A, Lee RW, Lemieux C, Gray MW and Burger G.
# "The complete mitochondrial DNA sequence of Scenedesmus obliquus reflects an
# intermediate stage in the evolution of the green algal mitochondrial genome."
# Genome Research (in press). 

Genetic Code [22]

Scenedesmus obliquus Mitochondrial

AAs  =   FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = -----------------------------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
# Genetic Code Table
#
# Obtained from: http://www3.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#  Version 3.7
#     Added code 23 Thraustochytrium mitochondrial code
#        formerly OGMP code 93
#        submitted by Gertraude Berger, Ph.D.
#
# This code has been created for the mitochondrial genome of the labyrinthulid
# Thraustochytrium aureum sequenced by the The Organelle Genome Megasequencing
# Program (OGMP).
#
# It is the similar to the bacterial code (trans_table 11) but it contains an
# additional stop codon (TTA) and also has a different set of start codons.

Genetic Code [23]

Thraustochytrium Mitochondrial

AAs  =   FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = --------------------------------M--M---------------M------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
