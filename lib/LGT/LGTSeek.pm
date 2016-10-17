
=head1 NAME

LGTSeek - Find Lateral Gene Transfer in sequencing data

=head1 SYNOPSIS

Need to put something useful here

=head1 DESCRIPTION

A module to run computes and process the output of data for purposes
of finding putative lateral gene transfer.

=head1 AUTHOR - Shaun Adkins, David R. Riley & Karsten B. Sieber

e-mail: sadkins@som.umaryland.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

=head2 LGTSeek.pm

 Title   : LGTSeek.pm
 Usage   : Suite of subroutines used to identify LGT
 Instance Routines:
        new                 : Create lgtseek object
        downloadSRA         : Download SRA file
        dumpFastq           : Convert sra file 2 fastq
        downloadCGHub       : Download bam from GC-Hub
        downloadEGA         : Download data from EGA
        prelim_filter       : Filter for potential LGT reads from human mapped bam
        runBWA              : Execute BWA aln mapping
        bwaPostProcess      : Filter LGT reads from host and donor
        prinseqFilterBam    : Make a new bam by removing low complexity and duplicate reads
        filter_bam_by_ids   : Make a new bam by filtering for read-ids
        sam2fasta           : Make fasta file from a bam
        splitBam            : Split 1 bam into multiple bams
        blast2lca           : Generate LCA using blast
        bestBlast2          : Run Blast keeping besthits only
        runLgtFinder        : Process blast results for LGT reads
        getGiTaxon          : Create gi2taxon object
        mpileup             : Make mpileup file for a bam
        empty_chk           : Check if an file is empty
        fail                : Gracefully error out
        _create_flag        : Create a sam flag
        _run_cmd            : Run a unix command

Class Routines:
		print_tab			: Print tab-delimited file of filtered counts by category
=cut

package LGT::LGTSeek;
our $VERSION = '1.3';	# SAdkins - Updated Version # since I rewrote a bit of it
use warnings;
use strict;
use version;

# Dependencies
use Carp;
$Carp::MaxArgLen = 0;
use File::Basename;
use Data::Dumper;
use GiTaxon;
use LGT::LGTBestBlast;
use LGT::LGTFinder;
use LGT::LGTbwa;
use LGT::Common;
use XML::Simple qw(:strict);
use Cwd;
$| = 1;

=head2 new

 Title   : new
 Usage   : my $lgtseek = LGTSeek->new({fastq => $fastq,...})
 Function: Creates a new LGTSeek object.
 Returns : An instance of LGTSeek
 Args    : A hash containing potentially several config options:

           fastq - Path fastq files.
           host_fasta - One or more fasta files to use a host references
           donor_fasta - One or more fasta files to use as donor references
           prinseq_bin - Path to prinseq perl script
           bwa_bin - Path to bwa binary
           aspera_path - Path to Aspera ascp binary
           taxon_dir - Directory where taxon information lives
           taxon_idx_dir - Directory where taxon indices can be created (or where they are)
           taxon_host - Hostname of taxon mongo database

=cut

sub new {
    my ( $class, $args ) = @_;

    my $self = {};

    bless $self, $class;
    $self->_init($args);
	return $self;
}

sub _init {
    my ( $self, $args ) = @_;

    # Useful suffix list for fileparse: @{$lgtseek->{'list'}}
    my $suffix_hash = {
        sam_suffix_list => [ '.sam.gz', '.sam' ],
        bam_suffix_list => [
            '_resorted.\d+.bam', '_resorted\.bam',
            '\.gpg\.bam',        '_prelim\.bam',
            '_name-sort\.bam',   '_pos-sort\.bam',
            '_psort\.bam',       '-psort\.bam',
            '\.psort\.bam',      '\.psrt\.bam',
            '_nsort\.bam',       '\.nsrt\.bam',
            '\.srt\.bam',        '\.sorted\.bam',
            '.bam'
        ],
        fastq_suffix_list => [
            qr/_[12]{1}\.f\w{0,3}q(.gz)?/,
            qr/_[12]{1}(\.\w+)?\.f\w*q(.gz)?/,
            qr/((_[12]{1})?\.\w+)?\.f\w*q(.gz)?/,
            '\.fastq\.gz',
            '\.f\w{0,3}q'
        ],
        fasta_suffix_list   => [ qr/.f\w{3}a(.gz)?/, '.fasta',        '.fa' ],
        mpileup_suffix_list => [ '.mpileup',         '_COVERAGE.txt', '.txt' ],
        suffix_regex        => qr/\.[^\.]+/
      };

    # add each key to current object
    foreach my $key ( keys %{$suffix_hash} ) {
        $self->{$key} = $suffix_hash->{$key};
    }

	# Do the same with the passed in args
	foreach my $key ( keys %{$args} ) {
		$self->{$key} = $args->{$key};
	}
}

=head2 getGiTaxon

 Title   : getGiTaxon
 Usage   : my $gi2tax = $lgtseek->getGiTaxon({'host' => 'foobar.com'...});
 Function: Retrieve a GiTaxon object ready to assign taxonomic information
 Returns : A GiTaxon object
 Args    : A hash options to pass to GiTaxon. These could have been
           passed in globally. They can also be overridden here.

           taxon_dir - Directory where taxon information lives
           taxon_idx_dir - Directory where taxon indices can be created (or where they are)
           taxon_host - Hostname of taxon mongo database

=cut

sub getGiTaxon {
    my ( $self, $config ) = @_;
    if ( $self->{verbose} ) {
        print STDERR "======== &getGiTaxon: Start ========\n";
    }

    # If we already have a gitaxon object we'll just return it.
    if ( !$self->{gitaxon} ) {

        # Apply any config options that came over.
        if ($config) {
            map { $self->{$_} = $config->{$_} } keys %$config;
        }

        # Create the object.
        $self->{gitaxon} = GiTaxon->new(
            {
                'taxon_dir'  => $self->{taxon_dir},
                'chunk_size' => 10000,
                'idx_dir'    => $self->{taxon_idx_dir},
                'host'       => $self->{taxon_host},
                'type'       => 'nucleotide',
                'verbose'    => $self->{verbose},
            }
        );
    }
    if ( $self->{verbose} ) {
        print STDERR "======== &getGiTaxon: Finished ========\n";
    }
    return $self->{gitaxon};
}

=head2 &prinseqFilterBam

 Title   : prinseqFilterBam
 Usage   : my $filteredBam = $LGTSeek->prinseqFilterBam({'input_bam' => '/path/to/file.bam'...})
 Function: Prinseq filter a bam file.
 Returns : A Hash $ref->{bam} = path to the filtered bam. $ref->{count} = # of reads passing filtering.
 Args    :
        input_bam    => /path/to/file.bam
        output_dir   => /path/for/output.bam
        overwrite    => <0|1> [0] 1= Overwrite output if it is found already.
        prinseq_bin  =>
        Picard_jar   =>
        java_opts    =>
        samtools_bin =>
=cut

sub prinseqFilterBam {
    my ( $self, $config ) = @_;
	if ( !$config->{input_bam} ) {
        $self->fail(
            "*** Error *** Must pass &prinseqFilterBam an input_bam =>\n");
    }
    if ( $self->empty_chk( { input => $config->{input_bam} } ) == 1 ) {
        $self->fail(
			"*** Error ***: &prinseqFilterBam input: $config->{input_bam} is empty.\n");
    }
    if ( $self->{verbose} ) {
        print STDERR "======== &prinseqFilterBam: Start ========\n";
    }

    # Override if it is provided
    $self->{prinseq_bin} =
      defined $config->{prinseq_bin}
      ? $config->{prinseq_bin}
      : $self->{prinseq_bin};
    $self->{Picard_jar} =
      defined $config->{Picard_jar}
      ? $config->{Picard_jar}
      : $self->{Picard_jar};
    $self->{java_opts} =
      defined $config->{java_opts} ? $config->{java_opts} : $self->{java_opts};
	unless ( defined($self->{java_bin})){
		$self->{java_bin} =
    	  defined $config->{java_bin} ? $config->{java_bin} : "/usr/bin/java";
  	}
	$self->{samtools_bin} = '/usr/bin/samtools' unless (defined $self->{samtools_bin});
    $self->{dedup} =
      defined $config->{dedup} ? $config->{dedup} : $self->{dedup};
    $self->{rm_low_cmplx} =
      defined $config->{rm_low_cmplx}
      ? $config->{rm_low_cmplx}
      : $self->{rm_low_cmplx};
    $self->{lc_method} =
      defined $config->{lc_method} ? $config->{lc_method} : "dust";
    $self->{lc_threshold} =
      defined $config->{lc_threshold} ? $config->{lc_threshold} : "7";

    my $overwrite = $config->{overwrite} ? $config->{overwrite} : 0;
    if ( $config->{output_dir} ) {
        $self->_run_cmd("mkdir -p $config->{output_dir}");
    }

    if ( !$self->{prinseq_bin} ) {
        $self->fail(
"*** Error *** Must provide a prinseq_bin parameter to run prinseq filtering\n"
        );
    }

    my $retval;
    if ( $self->{paired_end} ) {
        $retval = $self->_prinseqFilterPaired(
            $config->{input_bam}, $config->{output_dir},
            $config->{tmp_dir},   $overwrite
        );
    }
    else {
        $self->fail("*** Error *** Single end is currently not implemented\n");
    }
    if ( $self->{verbose} ) {
        print STDERR "======== &prinseqFilterBam: Finished ========\n";
    }
    return $retval;
}

=head2 &_prinseqFilterPaired

 Title   : _prinseqFilterPaired
 Usage   : *PRIVATE*
 Function: Prinseq filter a bam paired end file
 Returns : Path to the filtered bam file
 Args    : @_=($self,<bam_to_filter>,<output_dir>)

=cut

sub _prinseqFilterPaired {
    my ( $self, $original_bam, $output_dir, $temp_d, $overwrite ) = @_;
    my $bam_file = $original_bam;
    my ( $name, $path, $suff ) =
      fileparse( $bam_file, @{ $self->{bam_suffix_list} } );

    $output_dir = $output_dir ? $output_dir : $path;

    my $tmp_dir = defined $temp_d ? $temp_d : "$output_dir/tmp";

    # Create directories if they don't exist
    $self->_run_cmd("mkdir -p $output_dir") unless ( -d $output_dir );
    $self->_run_cmd("mkdir -p $tmp_dir")    unless ( -d $tmp_dir );

    my $prinseq_bin = $self->{prinseq_bin};
    my $samtools    = $self->{samtools_bin};
	my $java_opts_add = defined $self->{java_opts} ? "$self->{java_opts}" : "";
    my $Picard =
      "$self->{java_bin} $java_opts_add -jar $self->{Picard_jar}";
    my $dedup = defined $self->{dedup} ? $self->{dedup} : "1";
    my $rm_low_cmplx =
      defined $self->{rm_low_cmplx} ? $self->{rm_low_cmplx} : "1";
    my $lc_method    = $self->{lc_method};
    my $lc_threshold = $self->{lc_threshold};
    my $cmd;
    my $filtered;

    if ( -e "$output_dir/$name\_bad_ids.out" && $overwrite == 0 ) {
        if ( $self->{verbose} ) {
            print STDERR
"Already found the output for &prinseqFilter: $output_dir/$name\_bad_ids.out";
        }
        $filtered = $self->filter_bam_by_ids(
            {
                input_bam => $original_bam,
                header_comment =>
                  "\@CO\tID:PrinSeq-filtered\tPG:LGTseek\tVN:$VERSION",
                bad_list => "$output_dir/$name\_bad_ids.out",
            }
        );
        return $filtered;
    }
    else {

# Generate concatenated fastq files for prinseq derep filtering  ## Need to incorporate sam2fasta.pm KBS 01.07.14
        if ( $dedup == 1 ) {
            if ( $self->{verbose} ) {
                print STDERR "========= Deduplication Filtering =========\n";
            }
            if (
                $self->_run_cmd(
                    "$self->{samtools_bin} view -H $bam_file | head -n 1") !~
                /queryname/
              )
            {
                $self->_run_cmd(
                    "$self->{samtools_bin} sort -n -o $tmp_dir/$name\.bam $bam_file");
                $bam_file = "$tmp_dir/$name\.bam";
            }
            $self->_run_cmd(
"$Picard FixMateInformation I=$bam_file TMP_DIR=$tmp_dir SO=coordinate ASSUME_SORTED=1 VALIDATION_STRINGENCY=SILENT"
            );
            $self->_run_cmd(
"$Picard MarkDuplicates I=$bam_file TMP_DIR=$tmp_dir OUTPUT=$tmp_dir/$name\_dedup.bam METRICS_FILE=$tmp_dir/$name\_dedup-metrics.txt REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT"
            );
            open( my $BAM, "$self->{samtools_bin} view $tmp_dir/$name\_dedup.bam |" )
              or $self->fail(
"*** Error *** &prinseqFilterBam unable to open the deduped bam: $tmp_dir/$name\_dedup.bam"
              );
            open( my $BADS, "> $tmp_dir/$name\_derep_bad_ids.out" )
              or $self->fail(
"*** Error *** &prinseqFilterBam unable to open output file for dedup bad ids: $tmp_dir/$name\_derep_bad_ids.out"
              );
            while ( my $bam_line = $self->_read_bam_line($BAM) ) {
                if ( $bam_line->{flag}->{pcrdup} ) {
                    print $BADS "$bam_line->{id}\n";
                }
            }
            close $BAM;
            close $BADS;
        }

        if ( $rm_low_cmplx == 1 ) {

            # Generate single-read fastq for low complexity filtering
            if ( $self->{verbose} ) {
                print STDERR "========= Low Complexity Filter ========\n";
            }
            $self->_run_cmd(
"$Picard SamToFastq INPUT=$bam_file FASTQ=$tmp_dir/$name\_1.fastq SECOND_END_FASTQ=$tmp_dir/$name\_2.fastq VALIDATION_STRINGENCY=SILENT"
            );

            # Run prinseq for low complexity filtering
            $self->_run_cmd(
"/usr/bin/perl $prinseq_bin --fastq=$tmp_dir/$name\_1.fastq --out_good null --out_bad=$tmp_dir/$name\_lc_1_bad -lc_method $lc_method -lc_threshold $lc_threshold"
            );

            if ( -e "$tmp_dir/$name\_lc_1_bad.fastq" ) {

                # Pull out bad ids
                if ( $self->{verbose} ) {
                    print STDERR
                      "========= Pull Low-Cmplx-1 Bad ID's ========\n";
                }
                $self->_run_cmd(
"/usr/bin/perl -e 'while(<>){s/\\@//;s/\\/\\d//;print;<>;<>;<>;}' $tmp_dir/$name\_lc_1_bad.fastq > $tmp_dir/$name\_lc_1_bad_ids.out"
                );
            }
            else {
                if ( $self->{verbose} ) {
                    print STDERR
                      "Didn't find any low complexity sequences in read 1\n";
                }
                $self->_run_cmd("touch $tmp_dir/$name\_lc_1_bad_ids.out");
            }

            # Run prinseq for low complexity filtering
            $self->_run_cmd(
"/usr/bin/perl $prinseq_bin --fastq=$tmp_dir/$name\_2.fastq --out_good null --out_bad=$tmp_dir/$name\_lc_2_bad -lc_method dust -lc_threshold 7"
            );

            # Pull out bad ids
            if ( -e "$tmp_dir/$name\_lc_2_bad.fastq" ) {
                if ( $self->{verbose} ) {
                    print STDERR
                      "========= Pull Low-Cmplx-2 Bad ID's ========\n";
                }
                $self->_run_cmd(
"/usr/bin/perl -e 'while(<>){s/\\@//;s/\\/\\d//;print;<>;<>;<>;}' $tmp_dir/$name\_lc_2_bad.fastq > $tmp_dir/$name\_lc_2_bad_ids.out"
                );
            }
            else {
                if ( $self->{verbose} ) {
                    print STDERR
                      "Didn't find any low complexity sequences in read 2\n";
                }
                $self->_run_cmd("touch $tmp_dir/$name\_lc_2_bad_ids.out");
            }
        }

        # Merge bad ids from derep and lc filtering
        $cmd = "cat";
        if ( $dedup == 1 ) {
            $cmd = $cmd . " $tmp_dir/$name\_derep_bad_ids.out";
        }
        if ( $rm_low_cmplx == 1 ) {
            $cmd = $cmd
              . " $tmp_dir/$name\_lc_1_bad_ids.out $tmp_dir/$name\_lc_2_bad_ids.out";
        }
        $cmd = $cmd . " | sort -u -o $tmp_dir/$name\_prinseq-bad-ids.out";
        $self->_run_cmd($cmd);

        # Finally, filter based on dedup & lc bad ids
        $filtered = $self->filter_bam_by_ids(
            {
                input_bam  => $original_bam,
                output_dir => $output_dir,
                header_comment =>
                  "\@CO\tID:PrinSeq-filtered\tPG:LGTseek\tVN:$VERSION",
                bad_list => "$tmp_dir/$name\_prinseq-bad-ids.out",
            }
        );

        return $filtered;
    }
}

=head2 &sam2Fasta
## This needs to be updated to incorporate the sam2fasta.pm into LGTSeek.pm OR use sam2fastaconverter.pm
 Title   : sam2Fasta
 Usage   : my $fastas = $LGTSeek->sam2Fasta({'input' => '/path/to/file.bam'...})
 Function: Convert a bam/sam file to a fasta file
 Returns : a list of fasta/fastq files
 Args    : input => sam or bam file to convert to fasta
           output_dir => directory for output

=cut

sub sam2Fasta {
    my ( $self, $config ) = @_;
    my $bin = $self->{ergatis_bin};

    my $output_dir =
      $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    if ( $self->{verbose} ) {
        print STDERR "======== &sam2Fasta: Start ========\n";
    }

    # Make sure the output directory is present
    $self->_run_cmd("mkdir -p $output_dir");

    my $outfile;
    my ( $name, $path, $suff ) = fileparse( $config->{input},
        (qr/_filtered.bam$||.bam$||.sam$||.sam.gz$/) );
    my $cmd =
"/usr/bin/perl $bin/sam2fasta.pl --samtools_bin=$self->{samtools_bin} --input=$config->{input}";
    if ( $config->{fastq} ) {
        $outfile = "$output_dir/$name.fastq";
        $cmd .= " --fastq=1 --output_file=$outfile";
    }
    else {
        $outfile = "$output_dir/$name.fasta";
        $cmd .= " --fastq=0 --output_file=$outfile";
    }
    if ( $config->{combine_mates} ) {
        $cmd .= " --combine_mates=0";
    }
    if ( $config->{paired} || $self->{paired_end} ) {
        $cmd .= " --paired=1";
    }

    $self->_run_cmd($cmd);
    if ( $self->{verbose} ) {
        print STDERR "======== &sam2Fasta: Finished ========\n";
    }
    return $outfile;
}

=head2 &_run_cmd

 Title   : _run_cmd
 Usage   : *PRIVATE*
 Function: Run a unix command and fail if something goes wrong
 Returns : void
 Args    : Command to run

=cut

sub _run_cmd {

    my ( $self, $cmd ) = @_;

    if ( $self->{verbose} ) { print STDERR "CMD: $cmd\n"; }
    my $res = `$cmd`;
    if ($?) {
        print STDERR "FAIL_CMD: $cmd died with message: $res\n";
        print STDERR "Pausing 1 min and trying to run the cmd again.\n";
        sleep 60;
        chomp( $res = `$cmd` );
        if ($?) {
            $self->fail("*** Error *** $cmd died with message:\n$res\n\n");
        }
    }
    return $res;
}

=head2 &downloadSRA

 Title   : downloadSRA
 Usage   : $lgtseek->downloadSRA(({'experiment_id'} = 'SRX01234'})
 Function: Download sra files from the sequence read archive
 Returns : A list of the downloaded file paths
 Args    :

=cut

sub downloadSRA {
    my ( $self, $config ) = @_;
    if ( $self->{verbose} ) {
        print STDERR "======== &downloadSRA: Start ========\n";
    }

    # Check for all the aspera related options
    $self->{aspera_host} =
      $config->{aspera_host} ? $config->{aspera_host} : $self->{aspera_host};
    $self->{aspera_user} =
      $config->{aspera_user} ? $config->{aspera_user} : $self->{aspera_user};
    $self->{aspera_params} =
        $config->{aspera_params}
      ? $config->{aspera_params}
      : $self->{aspera_params};
    $self->{aspera_path} =
      $config->{aspera_path} ? $config->{aspera_path} : $self->{aspera_path};
    $self->{aspera_rate} =
      $config->{aspera_rate} ? $config->{aspera_rate} : $self->{aspera_rate};
    $self->{output_dir} =
      $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};

    my $retry_attempts =
      $config->{retry_attempts} ? $config->{retry_attempts} : 10;
    if ( !$self->{aspera_host} ) {
        $self->{aspera_host} = 'ftp-private.ncbi.nlm.nih.gov';
    }
    if ( !$self->{aspera_user} ) {
        $self->{aspera_user} = 'anonftp';
    }
    if ( !$self->{aspera_rate} ) {
        $self->{aspera_rate} = '200M';
    }

    if ( !$self->{aspera_path} ) {
        $self->fail(
"*** Error *** Need to specify an aspera_path (where is aspera installed) in order to download from the sra\n"
        );
    }
    if ( !$self->{output_dir} ) {
        $self->fail(
"*** Error *** Need to specify an output_dir in order to download from the sra\n"
        );
    }

    my $prefix;
    my $thousand;
    my $output_dir =
      $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    my $path_to_file;

    # We can pass an experiment_id, run_id or a full path to download
    if ( $config->{experiment_id} ) {
        my $exp_id = $config->{experiment_id};
        $exp_id =~ /^((.{3}).{3}).*/;
        $prefix   = $2;
        $thousand = $1;
        $path_to_file =
          "/sra/sra-instant/reads/ByExp/sra/$prefix/$thousand/$exp_id";
        $output_dir = "$output_dir/$prefix/$thousand/";
    }
    if ( $config->{run_id} ) {
        $config->{run_id} =~ /^((.{3}).{3}).*/;
        $prefix   = $2;
        $thousand = $1;
        $path_to_file =
"/sra/sra-instant/reads/ByRun/sra/$prefix/$thousand/$config->{run_id}";
        $output_dir = "$output_dir/$prefix/$thousand/";
    }
    elsif ( $config->{path} ) {
        $path_to_file = $config->{path};
    }

    # Make sure the output directory is present
    $self->_run_cmd("mkdir -p $output_dir");

    my $cmd_string =
"$self->{aspera_path}/bin/ascp -QTd -l$self->{aspera_rate} -i $self->{aspera_path}/etc/asperaweb_id_dsa.putty $self->{aspera_user}\@$self->{aspera_host}:$path_to_file $output_dir -L $output_dir -o Overwrite=diff 2>&1";

    #Retry the download several times just incase.
    my $retry = 1;

    my $retries = 0;
    while ($retry) {

        # Doing this echo y to ensure we accept any certs.
        my $out = $self->_run_cmd("echo y | $cmd_string");

        # We can actually exit non-0 and still succeed if the
        if ( $out =~ /Error/ ) {
            print STDERR
"Had a problem downloading $self->{aspera_host}:$path_to_file to $output_dir\n";
            print STDERR "$cmd_string";
            if ( $retries < $retry_attempts ) {
                $retries++;
                sleep $retries * 2;    # Sleep for 2 seconds per retry.
            }
            else {
                print STDERR
                  "Retries exhausted. Tried $retries times.\n$cmd_string";
                exit(1);
            }
        }
        else {
            $retry = 0;
            print STDERR
"$? $out Successfully downloaded $self->{aspera_host}:$path_to_file to $output_dir\n";
        }
    }

    my @files = `find $output_dir -name *.sra`;
    if ( $self->{verbose} ) {
        print STDERR "======== &downloadSRA: Finished ========\n";
    }
    return \@files;
}

=head2 &dumpFastq

 Title   : dumpFastq
 Usage   : $lgtseek->dumpFastq({'sra_file' => 'SRX01234'})
 Function: Run the sratoolkit program dump-fastq
 Returns : An object with a path and basename of the output as well as a list of output files
 Args    : An object with element 'sra_file' and optionally the path to the sratoolkit install

=cut

sub dumpFastq {
    my ( $self, $config ) = @_;

    $self->{sratoolkit_path} =
        $config->{sratoolkit_path}
      ? $config->{sratoolkit_path}
      : $self->{sratoolkit_path};
    if ( $self->{verbose} ) {
        print STDERR "======== &dumpFastq: Start ========\n";
    }

    # If we don't have a path provided we'll hope it's in our path.
    my $fastqdump_bin =
      $self->{sratoolkit_path}
      ? "$self->{sratoolkit_path}/fastq-dump"
      : "fastq-dump";

    $config->{sra_file} =~ s/\/\//\//g;

# Need to pull the version of the sratoolkit to determine if we need the --split-3 parameter.
    my $ret = `$fastqdump_bin -V`;
    my $version;
    my $cutoff_version;
    if ( $ret =~ /fastq-dump : ([\d.]+)/ ) {
        $version        = version->parse($1);
        $cutoff_version = version->parse('2.1.0');
    }
    else {
        $self->fail("*** Error *** $? $ret $fastqdump_bin\n");
    }
    if ( $version > $cutoff_version ) {

        $fastqdump_bin .= " --split-3 ";
    }
    $self->{output_dir} =
      $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    if ( !$self->{output_dir} ) {
        $self->fail(
"*** Error *** Need to specify an output_dir in order to download from the sra\n"
        );
    }

    my ( $name, $path, $suff ) = fileparse( $config->{sra_file}, ".sra" );
    chomp $name;
    my $res = $self->_run_cmd(
        "find $self->{output_dir} -maxdepth 1 -name '$name*.fastq'");
    my @files = split( /\n/, $res );

    if ( !@files && !$config->{overwrite} ) {
        my $cmd = "$fastqdump_bin -O $self->{output_dir} $config->{sra_file}";
        $self->_run_cmd($cmd);
        my $res = $self->_run_cmd(
            "find $self->{output_dir} -maxdepth 1 -name '$name*.fastq'");
        @files = split( /\n/, $res );
    }

    if ( $self->{verbose} ) { print STDERR "@files\n"; }

    my $retval = {
        'files'      => \@files,
        'path'       => $self->{output_dir},
        'basename'   => $name,
        'paired_end' => 0
    };
    if ( $files[0] =~ /_\d.fastq/ ) {
        $retval->{paired_end} = 1;
    }
    if ( $self->{verbose} ) {
        print STDERR "======== &dumpFastq: Finished ========\n";
    }
    return $retval;
}

=head2 downloadCGHub

 Title   : downloadCGHub
 Usage   : $lgtseek->downloadCGHub({
        analysis_id     => '00007994-abeb-4b16-a6ad-7230300a29e9',
        xml             => '/file/path/files_to_download.xml'
        output_dir      => '/file/path/for/output/dir/',
        threads         =>
        rate_limit      =>
        cghub_key       => '/path/to/key.file'

        });
 Function: Download TCGA bam files from the CGHub
 Returns : A list of the downloaded file paths
 Args    :

=cut

sub downloadCGHub {
    my ( $self, $config ) = @_;
    if ( $self->{verbose} ) {
        print STDERR "======== &downloadCGHub: Start ========\n";
    }

    # Check for all the genetorrent related options
    my $output_dir =
      defined $config->{output_dir}
      ? $config->{output_dir}
      : $self->{output_dir};
    $self->_run_cmd("mkdir -p $output_dir");
    my $genetorrent_path =
      defined $config->{genetorrent_path}
      ? $config->{genetorrent_path}
      : $self->{genetorrent_path};    ## Defaults to lgtseek.conf
    my $python =
      defined $self->{'python_2_7_path'}
      ? $self->{'python_2_7_path'}
      : "python";
    my $cghub_key =
      defined $config->{cghub_key} ? $config->{cghub_key} : $self->{cghub_key};
    my $max_retry_attempts =
      defined $config->{retry_attempts}
      ? $config->{retry_attempts}
      : $self->{retry_attempts};      ## Defaults to lgtseek.conf
    my $max_children =
      defined $config->{threads}
      ? $config->{threads}
      : $self->{threads};             ## Defaults to lgtseek.conf
    my $rate_limit =
      defined $config->{rate_limit}
      ? " -r $config->{rate_limit}"
      : "-r $self->{rate_limit}";     ## Defaults to lgtseek.conf

    if ( !$self->{cghub_key} && !$config->{cghub_key} ) {
        die "Need to specify the path to the cghub_key\n";
    }
    if ( !$self->{genetorrent_path} && !$config->{genetorrent_path} ) {
        die
"Need to specify an genetorrent_path (where is genetorrent/gtdownload installed) in order to download from CGHub\n";
    }
    if ( !$self->{output_dir} && !$config->{output_dir} ) {
        die "Need to specify an output_dir in order to download for CGHub\n";
    }

# We can pass an analysis_id (UUID), URI, .xml .gto to download. They can all be called analysis_id and will work the same way.
    my $download =
      $config->{analysis_id} ? $config->{analysis_id} : $config->{xml};

    # Retry the download several times just incase.
    my @bams_downloaded_list;
    my @bais_downloaded_list;
    my $retry               = 1;
    my $retry_attempts_made = 0;

  DOWNLOAD_BAM:
    while ( $retry == 1 && ( $retry_attempts_made <= $max_retry_attempts ) ) {
        ## Download the files
        if ( $self->{verbose} ) {
            print STDERR
"========= &downloadCGHub: Downloading: analysis_id\=$download . . . ========\n";
        }
        if ( $self->{verbose} ) {
            print STDERR
"CMD: $genetorrent_path\/gtdownload -l stdout -v -t -k 300 --max-children $max_children $rate_limit -p $output_dir -c $cghub_key -d $download \&>$output_dir/gtdownload.log\n";
        }
`$genetorrent_path\/gtdownload -l stdout -v -t -k 300 --max-children $max_children $rate_limit -p $output_dir -c $cghub_key -d $download \&>$output_dir/gtdownload.log`;
        if ($?) {
            $retry_attempts_made++;
            ## Every 3rd attempt delete and start over.
            if ( $retry_attempts_made % 3 ) {
                if ( $self->{verbose} ) {
                    print STDERR
                      "*** WARNING **** Download Failed!  ========\n";
                }
                if ( $self->{verbose} ) {
                    print STDERR
"========= &downloadCGHub: Retrying the download again.\n";
                }
                sleep 30;
                goto DOWNLOAD_BAM;
            }
            ## If the $retry_attempts_made is an even number we delete the original download and start it over cleanly
            else {
                if ( $self->{verbose} ) {
                    print STDERR
                      "*** WARNING **** Download Failed!  ========\n";
                }
                if ( $self->{verbose} ) {
                    print STDERR
"========= &downloadCGHub: Delete download and restart.\n";
                }
                my $bam_to_delete_string = $self->_run_cmd(
                    "find $output_dir -mindepth 1 -name \'*.bam\'");
                my @bam_to_delete_array = split( /\n/, $bam_to_delete_string );
                $self->_run_cmd("rm -rf $bam_to_delete_array[0]");
                $self->_run_cmd("rm -rf $bam_to_delete_array[0]\.bai")
                  if ( -e "$bam_to_delete_array[0]\.bai" );
                $self->_run_cmd("rm $output_dir/gtdownload.log");
                goto DOWNLOAD_BAM;
            }
        }
        else {
            $retry = 0;
        }
    }

    #  Making  a hash of the bam filename = md5 sum.
    ## This will allow retrying to download if the bam/bai are not finished or correct.
    if ( $config->{analysis_id} ) {
        if ( $self->{verbose} ) {
            print STDERR
"========= &downloadCGHub: Downloading cgquery.xml file to get md5 numbers.\n";
        }
        my $cmd =
"$python $genetorrent_path\/cgquery \"analysis_id\=$config->{analysis_id}\" -o $output_dir/cgquery.xml";
        if ( $self->{verbose} ) { print STDERR "$cmd\n"; }

        my $cgquery_exec_return = `$cmd`;
        if ($?) {
            print STDERR
"***Error*** :: $cmd :: failed with message: $cgquery_exec_return :: $?. Will now retry cgquery.\n";
            sleep 360;
            $cgquery_exec_return = `$cmd`;
        }
    }
    if ( !-e "$output_dir/cgquery.xml" ) {
        $self->fail("*** Error *** No $output_dir/cgquery.xml");
    }
    my $xml = $config->{xml} ? $config->{xml} : "$output_dir/cgquery.xml";
    my %cghub_md5_hash;
    my $ref = XMLin(
        $xml,
        KeyAttr => { Result => 'id' },
        ForceArray => [ 'Result', 'file' ],
        GroupTags => { files => 'file' }
      )
      or confess
"======== &downloadCGHub: &XMLin is unable to read in the cgquery.xml file: $xml because: $!\n";
    foreach my $sample ( sort keys %{ $ref->{'Result'} } ) {
        foreach my $file ( @{ $ref->{'Result'}->{$sample}->{files} } ) {
            if ( $file->{filename} =~ /\.bam$|\.bai$/ ) {
                $cghub_md5_hash{ $file->{filename} } =
                  $file->{checksum}->{content};
            }
        }
    }

    if ( $self->{verbose} ) {
        print STDERR
"========= &downloadCGHub: Finished downloading: $download. Checking md5's ========\n";
    }
    ## Find a list of the files downloaded
    my @files;
    foreach my $files_to_download ( keys %cghub_md5_hash ) {
        my $found_file =
          $self->_run_cmd("find $output_dir -name \'*$files_to_download\'");
        chomp($found_file);
        push( @files, $found_file );
    }

    ## Calculate the md5sum
    my $md5_fail_status = 0;
    foreach my $full_file_path (@files) {
        my ( $file_name, $path_to_file ) = fileparse($full_file_path);
        my $md5 =
          $self->_run_cmd("md5sum -b $full_file_path | cut -f1 -d \" \"");
        chomp($md5);
        if ( $md5 eq $cghub_md5_hash{$file_name} ) {
            if ( $self->{verbose} ) {
                print STDERR
"========= &downloadCGHub: md5 is correct for: $full_file_path ========\n";
            }
            if ( $file_name =~ /\.bam$/ ) {
                push( @bams_downloaded_list, $full_file_path );
            }
            elsif ( $file_name =~ /\.bai$/ ) {
                push( @bais_downloaded_list, $full_file_path );
            }
        }
        elsif ( $md5 ne $cghub_md5_hash{$file_name} ) {
            if ( $file_name =~ /\.bam$/ ) { $md5_fail_status++; }
            if ( $self->{verbose} ) {
                print STDERR
"========= &downloadCGHub: md5 Inconsistency for file: $full_file_path Calculated_md5: $md5 Expected-md5: $cghub_md5_hash{$file_name}.  ========\n";
            }
        }
    }

    if ( $md5_fail_status >= 1 ) {
        $self->fail(
"***Error*** &downloadCGHub: BAM download failed; the calculated md5 is inconsistent with the expected md5.\n"
        );
    }

    $self->_run_cmd("rm $output_dir/*.xml");
    $self->_run_cmd("rm $output_dir/*.gto");

    if ( $self->{verbose} ) {
        print STDERR "======== &downloadCGHub: Finished ========\n";
    }
    return \@bams_downloaded_list;
}

=head2 runBWA

 Title   : runBWA
 Usage   : $lgtseek->runBWA(({'base' => 'SRX01234','path' => '/path/to/files'})
 Function: Run bwa using the lgt_bwa wrapper
 Returns : The path to the bam file
 Args    : The input fasq/bam files and references which can be done a few different ways:

           # For files like /path/to/files/SRR01234_1.fastq and /path/to/files/SRR01234_2.fastq
           {
              'input_dir' => '/path/to/files/',
              'input_base' => 'SRR01234',
              'reference' => '/path/to/references/hg19.fa'
           }

           # For bam files and a list of references
           {
              'input_bam' => '/path/to/files/SRR01234.bam',
              'reference_list' => '/path/to/references/all_refs.list'
           }

=cut

sub runBWA {
    my ( $self, $config ) = @_;
    if ( $self->{verbose} ) {
        print STDERR "======== &runBWA: Start ========\n";
    }
    $self->{ergatis_bin} =
      $config->{ergatis_bin} ? $config->{ergatis_bin} : $self->{ergatis_bin};
    my $output_dir =
      $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};

# Check for a bwa path. If we don't have one we'll just hope it's in our global path.
    $self->{bwa_path} =
      $config->{bwa_path} ? $config->{bwa_path} : $self->{bwa_path};
    $self->{bwa_path} = $self->{bwa_path} ? $self->{bwa_path} : 'bwa';
    my $num_threads =
      defined $config->{threads} ? $config->{threads} : $self->{threads};
    $self->_run_cmd("mkdir -p $output_dir");

    my $conf = {
        num_aligns => 3,
        bwa_path   => $self->{bwa_path},
        output_dir => $output_dir,
        threads    => $num_threads,
    };

# Build the command string;
# my @cmd = ("$self->{ergatis_dir}/lgt_bwa --num_aligns=3 --bwa_path=$self->{bwa_path}");

    my $suff = '.sam';
    if ( $config->{output_bam} ) {
        $suff = '.bam';
        $conf->{output_bam} = 1;
    }

    my $basename = $config->{input_base};

    # Handle making up the lgt_bwa command with a bam file
    if ( $config->{input_bam} ) {
        if ( $self->empty_chk( { input => $config->{input_bam} } ) == 1 ) {
            print STDERR
              "*** Error *** &runBWA input: $config->{input_bam} is empty.\n";
            return [ $config->{input_bam} ];
        }
        my ( $name, $path, $suff ) =
          fileparse( $config->{input_bam}, ( "_prelim.bam", ".bam" ) );
        $basename           = $name;
        $conf->{input_base} = $basename;
        $conf->{input_bam}  = $config->{input_bam};
    }
    elsif ( $config->{input_dir} && $config->{input_base} ) {
        $conf->{input_dir}  = $config->{input_dir};
        $conf->{input_base} = $config->{input_base};
    }
    else {
        $self->fail(
"Must provide either a value to either input_bam or to input_base and input_dir\n"
        );
    }

    my $pre = '';
    if ( $config->{reference} ) {
        my ( $name, $dir, $suff ) =
          fileparse( $config->{reference}, qr/\.[^\.]+/ );
        $pre = "$name\_";
        $conf->{ref_file} = $config->{reference};
    }
    elsif ( $config->{reference_list} ) {
        $conf->{ref_file_list} = $config->{reference_list};
    }
    else {
        die
"Must provide a value for either reference or reference_list to run bwa\n";
    }

    $conf->{overwrite} = $config->{overwrite};
    map { $conf->{$_} = $config->{other_opts}->{$_}; }
      keys %{ $config->{other_opts} };
    $conf->{run_lca}     = $config->{run_lca};
    $conf->{lgtseek}     = $self;
    $conf->{cleanup_sai} = $config->{cleanup_sai};
    $conf->{out_file}    = $config->{out_file};
    LGT::LGTbwa::runBWA($conf);

    my @files = split( /\n/,
        $self->_run_cmd("find $output_dir -name '*$pre$basename$suff'") );
    map { chomp $_; } @files;
    if ( $self->{verbose} ) {
        print STDERR join( "\n", @files );
        print STDERR "\n";
    }
    if ( $self->{verbose} ) {
        print STDERR "======== &runBWA: Finished ========\n";
    }
    return \@files;
}

=head2 bwaPostProcess

 Title   : bwaPostProcess
 Usage   : $lgtseek->bwaPostProcess(({'donor_bams' => \@donors,'host_bams' => \@hosts})
 Function: Classify the results of a short read mapping (UM, UU, MM, UM_UM etc.)
 Returns : An object with counts of the different classes as well as the path to bam files
           containing these reads.
           $object->{files}->{           }
                            'lgt_donor'   => "$self->{output_dir}/".$prefix."lgt_donor.bam",
                            'lgt_host'    => "$self->{output_dir}/".$prefix."lgt_host.bam",
                            'integration_site_donor_donor' => "$self->{output_dir}/".$prefix."integration_site_donor_donor.bam",
                            'integration_site_donor_host' => "$self->{output_dir}/".$prefix."integration_site_donor_host.bam",
                            'microbiome_donor' => "$self->{output_dir}/".$prefix."microbiome.bam",
                            'output_dir' => '/dir/for/output.bams',
            $object->{counts}->{lgt}
            $object->{counts}->{microbiome}

 Args    : An object with donor and optionally host bam files.

=cut

sub bwaPostProcess {
    my ( $self, $config ) = @_;
    if ( $self->{verbose} ) {
        print STDERR "======== &bwaPostProcess: Start ========\n";
    }
    my $retval;

    # Do we have both donor and host bams?
    if ( $config->{donor_bam} && $config->{host_bam} ) {
		print STDERR "--- Detected donor and host paired-end files.\n";
        $retval = $self->_bwaPostProcessDonorHostPaired($config);
    } elsif ($config->{donor_bam} || $config->{host_bam}) {
		print STDERR "--- Detected just a single paired-end file.\n";
		$retval = $self->_bwaPostProcessSingle($config);
	} else {
		confess "***ERROR*** Detected no BAM files for either the donor or host!\n"
	}
    if ( $self->{verbose} ) {
        print STDERR "======== &bwaPostProcess: Finished ========\n";
    }
    return $retval;
}

sub _getPairedClass {
    my ( $self, $config ) = @_;

    my $fh         = $config->{fh};
    my $more_lines = 1;

    my $r1_class;
    my $r1_line;
    my $r2_class;
    my $r2_line;

    # Next establish the class of the donor read
    my $r1 = <$fh>;
    my $r2 = <$fh>;

    # Should check if these ended at the same time?
	# This means we are at the end of the file
    if ( !($r1 && $r2) ) {
        $more_lines = 0;
		print STDERR "Parsed to end of file\n";
        last;
    }

	chomp $r1;
	chomp $r2;

    if ( $config->{strip_xa} ) {
        $r1 =~ s/\tXA:Z\S+$//;
        $r2 =~ s/\tXA:Z\S+$//;
    }

    my $r1_flag = parse_flag( ( split( /\t/, $r1 ) )[1] );

    if ( !$r1_flag->{'qunmapped'} ) {
        $r1_line  = $r1;
        $r1_class = 'M';
    }
    elsif ( !$r1_class ) {
        $r1_line  = $r1;
        $r1_class = 'U';
    }
    if ( !$r1_flag->{'munmapped'} ) {
        $r2_line  = $r2;
        $r2_class = 'M';
    }
    elsif ( !$r2_class ) {
        $r2_line  = $r2;
        $r2_class = 'U';
    }

    my $class = "$r1_class$r2_class";

    return {
        class      => $class,
        r1_line    => $r1_line,
        r2_line    => $r2_line,
        more_lines => $more_lines
    };
}

sub _bwaPostProcessSingle {
     my ( $self, $config ) = @_;

     $self->{samtools_bin} =
       $self->{samtools_bin} ? $self->{samtools_bin} : '/usr/local/bin/samtools';
     my $samtools = $self->{samtools_bin};
     my $output_dir =
       $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};

     my $classes_each = {
         'MM' => 'paired',
         'UM' => 'single',
         'MU' => 'single',
         'UU' => 'none'
     };

     my $prefix = $config->{output_prefix} ? $config->{output_prefix} : '';

     # Use case requires only the single mappings and "none" mappings
     my $class_to_file_name = {
		 'single_map' => "$output_dir/" . $prefix . ".single_map.bam",
		 'no_map'	=> "$output_dir/" . $prefix . ".no_map.bam",
		 'single_paired_map' => "$output_dir/" . $prefix . ".single_paired_map.bam",
	 };

     my $class_counts = {
         'paired'  => undef,
         'single'  => undef,
         'none' => undef,
     };

     # Here are a bunch of file handles we'll use later.
     if ( $self->{verbose} ) {
         print STDERR "$output_dir/" . $prefix . ".single_map.bam\n";
         print STDERR "$output_dir/" . $prefix . ".no_map.bam\n";
         print STDERR "$output_dir/" . $prefix . ".single_paired_map.bam\n";
     }
     open(
         my $single_map,
         "| $samtools view -S -b -o $output_dir/" . $prefix . ".single_map.bam -"
     ) or die "Unable to open LGT single map file for writing\n";

     open(
         my $no_map,
         "| $samtools view -S -b -o $output_dir/" . $prefix . ".no_map.bam -"
     ) or die "Unable to open LGT 'no' map file for writing\n";

     open(
         my $single_paired_map,
         "| $samtools view -S -b -o $output_dir/" . $prefix . ".single_paired_map.bam -"
     ) or die "Unable to open LGT single/paired map file for writing\n";

	 # Perhaps in the future I can change these file names to rely on extensions like the Donor/Host subroutine relies on _donor and _host for assigning to the right file - SAdkins
     my $class_to_file = {
         'single_map'  => $single_map,
		 'no_map'	=> $no_map,
		 'single_paired_map' => $single_paired_map
     };

     my $bam = defined $config->{donor_bam} ? $config->{donor_bam} : $config->{host_bam};
	 if (! $bam) {
		confess "***ERROR*** Passed config neither has a defined donor nor host BAM";
	 }
     my $fh;
     my $head;

     # Open the BAM file for reading
     if ( $self->{verbose} ) { print STDERR "Opening $bam\n"; }
     if ( $bam =~ /.bam$/ ) {
         $head = `$samtools view -H $bam`;
         open( $fh, "-|", "$samtools view $bam" );
     }
     elsif ( $bam =~ /.sam.gz$/ ) {
         $head = `zcat $bam | $samtools view -H -S -`;
         open( $fh, "-|", "zcat $bam | $samtools view -S -" );
     }

     # Prime the files with headers.
     map {
        my @headers = split( /\n/, $head );
        print STDERR "Printing header to $_ file\n";
        print { $class_to_file->{$_} }
        join( "\n", grep( /^\@SQ/, @headers ) );
        print { $class_to_file->{$_} } "\n";
        my @pg_headers = grep( /^\@PG|^\@CO/, @headers );
        my %pg_hash;
        foreach my $pgs (@pg_headers) {
            chomp($pgs);
            $pg_hash{$pgs}++;
        }
		my $last_pg_id;
        if ( $headers[-1] =~ /^\@PG|^\@CO/ ) {
            my $last_pg = $headers[-1];
            $last_pg_id = ( split /\t/, $last_pg )[1];
            $last_pg_id =~ s/ID\:/PP\:/;
            $last_pg_id = "\t" . $last_pg_id;
        }
        foreach my $pgs ( keys %pg_hash ) {
            print { $class_to_file->{$_} } "$pgs\n";
        }
        print { $class_to_file->{$_} } "\@CO\tID:PostProcess\tPN:LGTseek\tVN:$VERSION";
		# Sometimes there are no previous PG-ID's so check for that.
		if (defined $last_pg_id) {
			print {$class_to_file->{$_}} "$last_pg_id\n";
		} else {
			print {$class_to_file->{$_}} "\n";
		}
     } keys %$class_to_file;

     my $more_lines = 1;
     my $line_num = 0;

     while ($more_lines) {
         # Get the class of the donor mappings
         my $obj = $self->_getPairedClass( { fh => $fh } );
         my $class = $obj->{class};
         $more_lines = $obj->{more_lines};
         my $r1_line = $obj->{r1_line};
         my $r2_line = $obj->{r2_line};

         if ($more_lines) {
             my $paired_class = $class;

  			 # print the single lines to the single_map file (if we are keeping this output file)
             if ( $classes_each->{$paired_class} eq "single" ) {
				 print { $class_to_file->{"single_map"} } "$r1_line\n$r2_line\n";
				 print { $class_to_file->{"single_paired_map"} } "$r1_line\n$r2_line\n";
             }

             if ( $classes_each->{$paired_class} eq "none" ) {
                 print { $class_to_file->{"no_map"} } "$r1_line\n$r2_line\n";
             }

             if ( $classes_each->{$paired_class} eq "paired" ) {
                 print { $class_to_file->{"single_paired_map"} } "$r1_line\n$r2_line\n";
             }

             # Increment the count for this class
             if ( $classes_each->{$paired_class} ) {
                 $class_counts->{ $classes_each->{$paired_class} }++;
             }
             # Increment the total count
             $line_num++;
         }
     }


     # Close up the file handles
     map {
         if ( $_ =~ /_map$/ ) {
             if ( $self->{verbose} ) { print STDERR "closing $_ file\n"; }
             close $class_to_file->{$_};
         }
     } keys %$class_to_file;

     # Set the total
     $class_counts->{total} = $line_num;

     # Return the files and the counts
     return {
         counts => $class_counts,
         files  => $class_to_file_name
    };
}

sub _bwaPostProcessDonorHostPaired {
     my ( $self, $config ) = @_;

     $self->{samtools_bin} =
       $self->{samtools_bin} ? $self->{samtools_bin} : '/usr/local/bin/samtools';
     my $samtools = $self->{samtools_bin};
     my $output_dir =
       $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};

     # Classes have the donor on the left and the host on the right
     my $classes_both = {
         'UM_MU' => 'lgt',
         'MU_UM' => 'lgt',
         'UM_MM' => 'integration_site_host',
         'MU_MM' => 'integration_site_host',
         'UU_UM' => 'integration_site_host',
         'UU_MU' => 'integration_site_host',
         'MM_UM' => 'integration_site_donor',
         'MM_MU' => 'integration_site_donor',
         'MU_UU' => 'integration_site_donor',
         'UM_UU' => 'integration_site_donor',
         'MM_UU' => 'microbiome',
         'UU_MM' => 'host',
         'UU_UU' => 'no_map',
         'MM_MM' => 'all_map',
         'UM_UM' => 'single_map',
         'MU_MU' => 'single_map'
     };

     my $prefix = $config->{output_prefix} ? $config->{output_prefix} : '';
     my $class_to_file_name = {
         'lgt_donor' => "$output_dir/" . $prefix . ".lgt_donor.bam",
         'lgt_host'  => "$output_dir/" . $prefix . ".lgt_host.bam",
         'integration_site_donor_donor' => "$output_dir/"
           . $prefix
           . ".integration_site_donor_donor.bam",
         'integration_site_donor_host' => "$output_dir/"
           . $prefix
           . ".integration_site_donor_host.bam",
         'microbiome_donor' => "$output_dir/" . $prefix . ".microbiome.bam",
     };

     my $class_counts = {
         'lgt'                    => undef,
         'integration_site_host'  => undef,
         'integration_site_donor' => undef,
         'microbiome'             => undef,
         'host'                   => undef,
         'no_map'                 => undef,
         'all_map'                => undef,
         'single_map'             => undef
     };

     # Here are a bunch of file handles we'll use later.
     if ( $self->{verbose} ) {
         print STDERR "$output_dir/" . $prefix . ".lgt_donor.bam\n";
     }
     open(
         my $lgtd,
         "| $samtools view -S -b -o $output_dir/" . $prefix . ".lgt_donor.bam -"
     ) or die "Unable to open LGT donor file for writing\n";
     open( my $lgth,
         "| $samtools view -S -b -o $output_dir/" . $prefix . ".lgt_host.bam -" )
       or die "Unable to open LGT recipient file for writing\n";
     open(
         my $int_site_donor_d,
         "| $samtools view -S -b -o $output_dir/"
           . $prefix
           . ".integration_site_donor_donor.bam -"
     ) or die "Unable to open donor integration site file for writing\n";
     open(
         my $int_site_donor_h,
         "| $samtools view -S -b -o $output_dir/"
           . $prefix
           . ".integration_site_donor_host.bam -"
     ) or die "Unable to open recipient integration site file for writing\n";
     open(
         my $microbiome_donor,
         "| $samtools view -S -b -o $output_dir/" . $prefix . ".microbiome.bam -"
     ) or die "Unable to open donor microbiome file for writing\n";

     my $class_to_file = {
         'lgt_donor'                    => $lgtd,
         'lgt_host'                     => $lgth,
         'integration_site_donor_donor' => $int_site_donor_d,
         'integration_site_donor_host'  => $int_site_donor_h,
         'microbiome_donor'             => $microbiome_donor
     };

     my $donor_bam = $config->{donor_bam};
     my $host_bam  = $config->{host_bam};
     my $donor_fh;
     my $host_fh;
     my $donor_head;
     my $host_head;

     # Open all the donor files
     if ( $self->{verbose} ) { print STDERR "Opening $donor_bam\n"; }
     if ( $donor_bam =~ /.bam$/ ) {
         $donor_head = `$samtools view -H $donor_bam`;
         open( $donor_fh, "-|", "$samtools view $donor_bam" );
     }
     elsif ( $donor_bam =~ /.sam.gz$/ ) {
         $donor_head = `zcat $donor_bam | $samtools view -H -S -`;
         open( $donor_fh, "-|", "zcat $donor_bam | $samtools view -S -" );
     }

    # Open all the host files
    if ( $self->{verbose} ) { print STDERR "Opening $host_bam\n"; }
    if ( $host_bam =~ /.bam$/ ) {
     	$host_head = `$samtools view -H $host_bam`;
     	open( $host_fh, "-|", "$samtools view $host_bam" );
    } elsif ( $host_bam =~ /.sam.gz$/ ) {
     	$host_head = `zcat $host_bam | $samtools view -H -S -`;
     	open( $host_fh, "-|", "zcat $host_bam | $samtools view -S -" );
    }

     # Prime the files with headers.
     map {
         if ( $_ =~ /_donor$/ ) {
             my @donor_headers = split( /\n/, $donor_head );
             print STDERR "Printing header to $_ donor file\n";
             print { $class_to_file->{$_} }
               join( "\n", grep( /^\@SQ/, @donor_headers ) );
             print { $class_to_file->{$_} } "\n";
             my @pg_headers = grep( /^\@PG|^\@CO/, @donor_headers );
             my %pg_hash;
             foreach my $pgs (@pg_headers) {
                 chomp($pgs);
                 $pg_hash{$pgs}++;
             }
			 my $last_pg_id;
             if ( $donor_headers[-1] =~ /^\@PG|^\@CO/ ) {
                 my $last_pg = $donor_headers[-1];
                 $last_pg_id = ( split /\t/, $last_pg )[1];
                 $last_pg_id =~ s/ID\:/PP\:/;
                 $last_pg_id = "\t" . $last_pg_id;
             }
             foreach my $pgs ( keys %pg_hash ) {
                 print { $class_to_file->{$_} } "$pgs\n";
             }
             print { $class_to_file->{$_} }
               "\@CO\tID:PostProcess\tPN:LGTseek\tVN:$VERSION";
			# Sometimes there are no previous PG-ID's so check for that.
			if (defined $last_pg_id) {
				print {$class_to_file->{$_}} "$last_pg_id\n";
			} else {
				print {$class_to_file->{$_}} "\n";
  			}
         }
         elsif ( $_ =~ /_host$/) {
     		# Host file does not exist when BWA is only aligned to donor
             my @host_headers = split( /\n/, $host_head );
             print STDERR "Printing header to $_ host file\n";
             print { $class_to_file->{$_} }
               join( "\n", grep( /^\@SQ/, @host_headers ) );
             print { $class_to_file->{$_} } "\n";
             my @pg_headers = grep( /^\@PG|^\@CO/, @host_headers );
             my %pg_hash;
             foreach my $pgs (@pg_headers) {
                 chomp($pgs);
                 $pg_hash{$pgs}++;
             }
             my $last_pg_id;
             if ( $host_headers[-1] =~ /^\@PG|^\@CO/ ) {
                 my $last_pg = $host_headers[-1];
                 $last_pg_id = ( split /\t/, $last_pg )[1];
                 $last_pg_id =~ s/ID\:/PP\:/;
                 $last_pg_id = "\t" . $last_pg_id;
             }
             foreach my $pgs ( keys %pg_hash ) {
                 print { $class_to_file->{$_} } "$pgs\n";
             }
             print { $class_to_file->{$_} }
               "\@CO\tID:PostProcess\tPN:LGTseek\tVN:$VERSION";
			# Sometimes there are no previous PG-ID's so check for that.
			if (defined $last_pg_id) {
				print {$class_to_file->{$_}} "$last_pg_id\n";
			} else {
				print {$class_to_file->{$_}} "\n";
         	}
		}
     } keys %$class_to_file;

     #   exit;
     my $more_lines = 1;

     my $line_num = 0;

     while ($more_lines) {

         # Get the class of the host mappings
         my $obj = $self->_getPairedClass( { fh => $host_fh } );
         my $hclass = $obj->{class};
         $more_lines = $obj->{more_lines};
         my $hr1_line = $obj->{r1_line};
         my $hr2_line = $obj->{r2_line};

         # Get the class of the donor mappings
         $obj = $self->_getPairedClass( { fh => $donor_fh } );
         my $dclass = $obj->{class};
         $more_lines = $obj->{more_lines};
         my $dr1_line = $obj->{r1_line};
         my $dr2_line = $obj->{r2_line};

         if ($more_lines) {
             my $paired_class = "$dclass\_$hclass";

   # print the donor lines to the donor file (if we are keeping this output file)
             if ( $class_to_file->{ $classes_both->{$paired_class} . "_donor" } )
             {
                 print {
                     $class_to_file->{ $classes_both->{$paired_class}
                           . "_donor" } } "$dr1_line\n$dr2_line\n";
             }

     # print the host lines to the host file (if we are keeping this output file)
             if ( $class_to_file->{ $classes_both->{$paired_class} . "_host" } )
             {
                 print {
                     $class_to_file->{ $classes_both->{$paired_class} . "_host" }
                 } "$hr1_line\n$hr2_line\n";
             }

             # Increment the count for this class
             if ( $classes_both->{$paired_class} ) {
                 $class_counts->{ $classes_both->{$paired_class} }++;
             }

             # Increment the total count
             $line_num++;
         }
     }

     # Close up the file handles
     map {
         if ( $_ =~ /_donor$/ ) {
             if ( $self->{verbose} ) { print STDERR "closing $_ donor file\n"; }
             close $class_to_file->{$_};
         }
         elsif ( $_ =~ /_host$/ ) {
             if ( $self->{verbose} ) { print STDERR "closing $_ host file\n"; }
             close $class_to_file->{$_};
         }
     } keys %$class_to_file;

     # Set the total
     $class_counts->{total} = $line_num;

     # Return the files and the counts
     return {
         counts => $class_counts,
         files  => $class_to_file_name
    };
}

=head2 blast2lca

 Title   : blast2lca
 Usage   : my $lca_file = $lgtseek->blast2lca({'blast' => 'path_to_blast-m8.txt','output_dir' => 'path_to_out_dir'})
 Function: Take a blast -m8 report and calculate the LCA'
 Returns : A object with the path to the output file(s) with LCA's.
            $lca_file->{independent} = Path to the file with the LCA's for each unique ID.
            $lca_file->{PE_lca}      = Path to the file with the conservative and liberal LCA for combining mate information
 Args    :
            blast               =>
            output_dir          =>
            evalue_cutoff       =>  ###  [1] Max evalue allowed for a hit. Example : 1e-5.
            best_hits_only      => <0|1> [0] 1= Parse the Blast file for only best hits.
            combine_PE_lca dd     => <0|1> [0] 1= Merge the PE reads with the conservative and liberal methods. ## NOT IMPLEMENTED YET ##

=cut

sub blast2lca {
    my ( $self, $config ) = @_;
    my $gi2tax = $config->{gi2tax};
    if ( $self->{verbose} ) {
        print STDERR "======== &blast2lca: Start ========\n";
    }
    open IN, $config->{blast} or confess "Unable to open " . $config->{blast} . "\n";
    my $output_dir =
      $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    $self->_run_cmd("mkdir -p $output_dir");

    my ( $name, $directories, $suffix ) =
      fileparse( $config->{blast}, qr/\.[^.]*/ );
    $name = $config->{output_prefix} ? $config->{output_prefix} : $name;
    open OUT, ">$output_dir/$name\_lca-independent.out"
      or confess
"*** Error *** &blast2lca unable to open $output_dir/$name\_lca-independent.out";
    if ( $config->{best_hits_only} == 1 ) {
        print OUT join( "\t", ( 'read_id', 'lca', 'best_evalue' ) );
    }
    else {
        print OUT
          join( "\t", ( 'read_id', 'lca', 'best_evalue', 'highest_evalue' ) );
    }    ## KBS
    print OUT "\n";
    my $hits_by_readname = {};
    my $id = '';
    my $evalue_cutoff = 1;
    my $best_evalue;
    my $worst_evalue;
    my $lca;

    my @vals;
    while (<IN>) {
        my @fields = split(/\t/);
        my $new_id = $fields[0];

        if ( $config->{evalue_cutoff} ) {
            $evalue_cutoff = $config->{evalue_cutoff};
        }
        next if ( $new_id eq $id && $fields[10] > $evalue_cutoff );
        my $taxon = $gi2tax->getTaxon( $fields[1] );
        $taxon->{scientific_name} =~ /^(\w+) /;
        my $genera = $1;
        if ( !defined($id) && $fields[10] < $evalue_cutoff ) {
            $id = $new_id;
            if ( $config->{best_hits_only} == 1 ) {
                $evalue_cutoff = $fields[10];
            }
            $best_evalue  = $fields[10];
            $worst_evalue = $fields[10];
            $lca          = $taxon->{lineage};
            $hits_by_readname->{ $fields[0] }->{ $taxon->{taxon_id} } =
              $taxon->{lineage};
        }
        elsif ( $new_id eq $id && $fields[10] <= $evalue_cutoff ) {
            $hits_by_readname->{ $fields[0] }->{ $taxon->{taxon_id} } =
              $taxon->{lineage};
            my $newlca = &find_lca( [ $lca, $taxon->{lineage} ] );
            $lca = $newlca;
            if ( $fields[10] > $worst_evalue ) { $worst_evalue = $fields[10]; }
        }
        elsif ( $new_id ne $id ) {
            if ( $config->{best_hits_only} == 1 ) {
                print OUT "$id\t$lca\t$best_evalue\n";
            }
            else {
                print OUT "$id\t$lca\t$best_evalue\t$worst_evalue\n";
            }
            $id = $new_id;
            if ( $config->{best_hits_only} == 1 ) {
                $evalue_cutoff = $fields[10];
            }
            $best_evalue  = $fields[10];
            $worst_evalue = $fields[10];
            $lca          = $taxon->{lineage};
        }
    }

    # Print out last line
    if ( $config->{best_hits_only} == 1 ) {
        print OUT "$id\t$lca\t$best_evalue\n";
    }
    else {
        print OUT "$id\t$lca\t$best_evalue\t$worst_evalue\n";
    }
    close OUT;

    if ( $config->{combine_PE_lca} == 1 ) {
        open( IN, "<", "$output_dir/$name\_lca-independent.out" )
          or confess
"*** Error *** &blast2lca unable to open input: $output_dir/$name\_lca-independent.out because: $!\n";
        open( OUT1, ">", "$output_dir/$name\_lca-conservative.out" )
          or confess
"*** Error *** &blast2lca unable to open output: $output_dir/$name\_lca-conservative.out because: $!\n";
        open( OUT2, ">", "$output_dir/$name\_lca-liberal.out" )
          or confess
"*** Error *** &blast2lca unable to open output: $output_dir/$name\_lca-liberal.out because: $!\n";
        my $continue = 1;
        while ( $continue == 1 ) {
            my $line1 = <IN>;
            my $line2 = <IN>;
            if ( !$line1 || !$line2 ) { $continue = 0; }
          NEXT:
            my @f1 = split( /\t/, $line1 );
            my @f2 = split( /\t/, $line2 );
            my $id1 = $f1[0];
            my $id2 = $f2[0];
            $id1 =~ /([A-Za-z0-9-.|:]+)(\_?\/?[1,2]?)/;
            my $id1_short = $1;
            $id2 =~ /([A-Za-z0-9-.|:]+)(\_?\/?[1,2]?)/;
            my $id2_short = $1;

            if ( $id1_short ne $id2_short ) {
                $line1 = $line2;
                $line2 = <IN>;
                goto NEXT;
            }
            my $lca1             = $f1[1];
            my $lca2             = $f2[1];
            my $conservative_lca = &find_lca( [ $lca1, $lca2 ] );

            # Calculate Liberal LCA
            my $liberal_lca;
            if ( $lca1 =~ $lca2 || $lca2 =~ $lca1 ) {
                if ( length($lca1) >= length($lca2) ) {
                    $liberal_lca = $lca1;
                }
                else {
                    $liberal_lca = $lca2;
                }
            }
            print OUT1 "$id1_short\t$conservative_lca\n";
            print OUT2 "$id1_short\t$liberal_lca\n";
        }
        close OUT1;
        close OUT2;
    }
    if ( $self->{verbose} ) {
        print STDERR "======== &blast2lca: Finished ========\n";
    }
    return "$output_dir/$name\_lca-independent.out";
}


=head2 &bestBlast2

 Title   : bestBlast2
 Usage   : $lgtseek->bestBlast2(({'fasta' => 'fasta_file', 'ref' => '/path/to/ref/db'}))
 Function: Run blast (or megablast) against a reference database to find best blast hits.
 Returns : Hash with:
            overall_blast => $overallfile,
            out1file_file => $out1file,
            out2file_file => $out2file
 Args    : An object with the input bame files and the reference database.
            fasta =>
            bam =>
            threads =>
            lineage1 =>
            lineage2 =>
            db =>
            output_dir =>

=cut

sub bestBlast2 {
    my ( $self, $config ) = @_;
    if ( $self->{verbose} ) {
        print STDERR "======== &bestBlast2: Start ========\n";
    }
    my $output_dir =
      defined $config->{output_dir}
      ? "$config->{output_dir}"
      : "$self->{output_dir}";
    $self->_run_cmd("mkdir -p $output_dir");

    # Convert bams to fasta
    my $fasta;
    if ( $config->{'bam'} ) {
        my $newfasta = $self->sam2Fasta(
            {
                input      => $config->{'bam'},
                output_dir => $config->{output_dir},
                paired     => 1
            }
        );
        $fasta = $newfasta;
    }
    elsif ( $config->{'fasta'} ) {
        $fasta = $config->{fasta};
    }

    # Blast fasta @ database and filter for best hits.
    my $files = LGT::LGTBestBlast->filterBlast(
        {
            blast_bin  => $config->{blast_bin},
            fasta      => $fasta,
            db         => $config->{db},
            threads    => $self->{threads},
            gitaxon    => $self->getGiTaxon( {} ),
            lineage1   => $config->{lineage1},
            lineage2   => $config->{lineage2},
            output_dir => $output_dir,
        }
    );

    if ( $self->{verbose} ) {
        print STDERR "======== &bestBlast2: Finished ========\n";
    }
    return $files;
}

=head2 &runLgtFinder

 Title   : runLgtFinder
 Usage   : $lgtseek->runLgtFinder(({'inputs' => \@files})
 Function: Take output from LGTBestBlast and look for putative LGT's both within reads and across read pairs
 Returns : A file with all the LGT's hit information
 Args    : An object from BestBlast2
            max_overlap => 20,
            min_length  => 0,
            ref_lineage => 'Homo',
            output_dir  => '/somewhere/'

=cut

sub runLgtFinder {
    my ( $self, $config ) = @_;
    if ( $config->{output_dir} ) {
        $self->_run_cmd("mkdir -p $config->{output_dir}");
    }
    my $ret = LGT::LGTFinder->findLGT($config);
    return $ret;
}

=head2 splitBam

 Title   : splitBam
 Usage   : $lgtseek->splitBam(({'input' => $file})
 Function: Split a bam file into smaller chunks.
 Returns : A list of the bam files split up
 Args    :
            input = An object with the input bam files
            seqs_per_file = # of seqs per file for each split file
            output_dir = Directory for output
            samtools_bin = bin directory with samtools

=cut

sub splitBam {
    my ( $self, $config ) = @_;
    if ( !$config->{input} ) { die "Must give &splitBam an input =>\n"; }
    if ( $self->empty_chk( { input => $config->{input} } ) == 1 ) {
        print STDERR
          "*** Warning *** &split_bam input: $config->{input} is empty.\n";
    }
    my $output_dir =
      $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    if ( $self->{verbose} ) {
        print STDERR "======== &splitBam: Start ========\n";
    }
    $self->_run_cmd("mkdir -p $output_dir");
    my $seqs_per_file =
      $config->{seqs_per_file} ? $config->{seqs_per_file} : 50000000;
    my $samtools = $self->{samtools_bin};

    # Parse out the pieces of the filename
    my ( $fn, $path, $suffix ) = fileparse( $config->{input}, '.bam' );

    # Pull out the header
    my $header = $self->_run_cmd("$samtools view -H $config->{input}");

    # Open the input file
    open( IN, "-|", "$samtools view $config->{input}" );

    # Open the first output file
    my $count = 0;
    my $ofile = "$output_dir/$fn\_$count.bam";
    open( my $ofh, "| samtools view -S -b -o $ofile -" );
    my @outfiles = ($ofile);

    # Print the header
    print $ofh $header;

    my $i = 0;
    while ( my $line = <IN> ) {
        my @fields = split( /\t/, $line );
        my $flag = parse_flag( $fields[1] );

        # Increment the counter
        $i++;

        # Strip out the XA tag
        $line =~ s/\s+XA:Z:\S+//;

        # Make sure we don't accidentally split read pairs.
        if ( $flag->{'last'} && $i >= $seqs_per_file ) {

            # Print out the line
            print $ofh $line;

            # Close the old file
            close $ofh;

            # Open the new file
            $count++;
            $ofile = "$output_dir/$fn\_$count.bam";
            open( $ofh, "| samtools view -S -b -o $ofile -" );
            push( @outfiles, $ofile );

            # Print the header to the new file
            print $ofh $header;

            # Reset the counter
            $i = 0;
            next;
        }

        # If we haven't filled the current file yet just keep printing.
        print $ofh $line;

    }
    close $ofh;

    ## Check to make sure the last file isn't empty.
    if ( $self->empty_chk( { input => $ofile } ) == 1 ) {
        $self->_run_cmd("rm $ofile");
        my $remove_empty_bam = pop(@outfiles);
        $count--;
    }

    if ( $self->{verbose} ) {
        print STDERR
"Split $config->{input} into $count bams, each with $seqs_per_file sequences per bam.\n";
    }
    my $files = \@outfiles;
    if ( $self->{verbose} ) {
        print STDERR "======== &splitBam: Finished ========\n";
    }
    return $files;
}

sub _create_flag {
    my ( $self, $data ) = @_;

    my $flag = 0;
    $flag += $data->{paired} *        ( 2**0 );
    $flag += $data->{proper} *        ( 2**1 );
    $flag += $data->{qunmapped} *     ( 2**2 );
    $flag += $data->{munmapped} *     ( 2**3 );
    $flag += $data->{qrev} *          ( 2**4 );
    $flag += $data->{mrev} *          ( 2**5 );
    $flag += $data->{first} *         ( 2**6 );
    $flag += $data->{last} *          ( 2**7 );
    $flag += $data->{secondary} *     ( 2**8 );
    $flag += $data->{failqual} *      ( 2**9 );
    $flag += $data->{pcrdup} *        ( 2**10 );
    $flag += $data->{supplementary} * ( 2**11 );

    return $flag;
}

=head2 decrypt

 Title   : decrypt
 Usage   : my $decrypted_bam=$lgtseek->decrypt(({'input' => <file.bam.gpg>, 'key' => <Path to key>})
 Function: Decrypt a .bam.gpg input file.
 Returns : An un-encrypted bam
 Args    :
    input       => encrypted bam for decryption
    url         => url to download the key from
    key         => path to key file
    output_dir  => directory for output

=cut

sub decrypt {
    my ( $self, $config ) = @_;
    if ( !$config->{input} || $config->{input} !~ /\.bam\.gpg$/ ) {
        die
          "*** Error *** Must pass an encrypted .bam.gpg to this subroutine.\n";
    }
    if ( !$config->{url} && !$config->{key} ) {
        die
"*** Error *** Must give an url to download the key from or the path to the key.\n";
    }
    my $output_dir =
      $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    if ( $self->{verbose} ) {
        print STDERR "======== &Decrypt: Start ========\n";
    }
    $self->_run_cmd("mkdir -p $output_dir");
    my ( $bam, $path, $suffix ) = fileparse( $config->{input}, '.gpg' );
    my $key;
    if ( $config->{url} ) {
        $self->_run_cmd("wget $config->{url}");
        $key = `find . -name '*.key' | cut -f2 -d '/'`;
    }
    else {
        $key = $config->{key};
    }
    $self->_run_cmd("gpg --import $key");
    $self->_run_cmd("gpg -o $output_dir/$bam -d $config->{input}");
    my $outfile = "$output_dir/$bam";
    $self->_run_cmd("rm $config->{input}");
    $self->_run_cmd("rm $key");
    if ( $self->{verbose} ) {
        print STDERR "======== &Decrypt: Finished ========\n";
    }
    return $outfile;
}

=head2 mpileup

 Title   : mpileup
 Usage   : my $mpileup_file=$lgtseek->mpileup({'input' => <bam>})
 Function: Calculate samtools mpileup
 Returns : File path to mpileup output.
 Args    :
    input       => unsorted bam
    srtd_bam     => position sorted bam input
    ref         =>Reference (Optional)
    max         => max per-BAM depth to avoid excessive memory usage [250] (Optional)

=cut

sub mpileup {
    my ( $self, $config ) = @_;
    if ( !$config->{input} && !$config->{srtd_bam} ) {
        die "*** Error *** Must give an input bam to calculat coverage on.\n";
    }
    my $input = $config->{input} ? $config->{input} : $config->{srtd_bam};
    if ( $self->empty_chk( { input => $input } ) == 1 ) {
        print STDERR "*** Warning *** &mpileup input: $input is empty.\n";
    }
    my ( $fn, $path, $suffix ) = fileparse( $config->{input}, '.bam' );
    my $output_dir =
      $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    my $output = "$output_dir/$fn\.mpileup";
    if ( $self->{verbose} ) {
        print STDERR "======== &mpileup: Start ========\n";
    }

    my $samtools = $self->{samtools_bin};
    my $srtd_bam;
    my $max_opt = $config->{max}     ? "-d $config->{max}" : "";
    my $ref_opt = $config->{ref}     ? "-f $config->{ref}" : "";

    my $cmd = "$samtools";

    ## If BAM isn't sorted, do that first
    if ( $config->{input} ) {
        $cmd .= " sort -m 5000000000 -o - $config->{input} | $samtools mpileup";
        $srtd_bam = "-";
    } else {
        $cmd      .= " mpileup";
        $srtd_bam = $config->{srtd_bam};
    }
    $cmd .= " -A -o $output $max_opt $ref_opt $srtd_bam";
    $self->_run_cmd($cmd);
    if ( $self->{verbose} ) {
        print STDERR "======== &mpileup: Finished ========\n";
    }
    return $output;
}

=head2 &prelim_filter

 Title   : prelim_filter
 Usage   : my $potential_LGT_bam = $lgtseek->prelim_filter({input => <bam>})
 Function: Removes M_M reads
 Returns : A bam will all other reads
 Args    : A hash containing potentially several config options:
        input_bam       => bam
        output_dir      => directory for output
        keep_softclip   => <0|1> [0] 1= Keep the soft clipped M_M reads
        overwrite       => <0|1> [1] 1= Overwrite the output if it already exists
        split_bam       => <0|1> [0] 1= Split by by #seqs_per_file
        seqs_per_file   => [50000000]
        name_sort_input => <0|1> [0] 1= Resort input bam on names

    ## TODO: I should be able to implement the following samtools command to create a list of reads that pass prelim_filter criteria, then re-parse the input file for the reads and split the bams acocordingly. This may be significantly faster than the current implementation: "samtools view -F 3844 $input -u | view -f0x8 -" --> picard SamToFastQ
=cut

sub prelim_filter {

    ## General code scheme: (1) Resort based on names. (2) Filter out M_M reads. (3) Then finally split into chunks.

    my ( $self, $config ) = @_;
    my $input =
      $config->{input_bam} ? $config->{input_bam} : $self->{input_bam};
    if ( $self->empty_chk( { input => $input } ) == 1 ) {
        print STDERR "*** Warning *** &prelim_filter input: $input is empty.\n";
    }
    if ( $input !~ /\.bam$/ ) {
        $self->fail(
"*** Error *** Input file: $input doesn't look like a bam without a .bam suffix."
        );
    }
    my $output_dir =
      $config->{output_dir} ? $config->{output_dir} : $self->{output_dir};
    $self->_run_cmd("mkdir -p $output_dir");

    my $keep_softclip =
      defined $config->{keep_softclip}
      ? $config->{keep_softclip}
      : $self->{keep_softclip};
    my $softclip_min =
      defined $config->{softclip_min}
      ? $config->{softclip_min}
      : $self->{softclip_min};
    my $overwrite =
      defined $config->{overwrite} ? $config->{overwrite} : $self->{overwrite};
    my $split_bam =
      defined $config->{split_bam} ? $config->{split_bam} : $self->{split_bam};
    my $prelim_filter =
      defined $config->{prelim_filter}
      ? $config->{prelim_filter}
      : $self->{prelim_filter};
    my $seqs_per_file =
        $config->{seqs_per_file}
      ? $config->{seqs_per_file}
      : $self->{seqs_per_file};
    my $name_sort_input =
      defined $config->{name_sort_input}
      ? $config->{name_sort_input}
      : $self->{name_sort_input};
    my $name_sort_check =
      defined $config->{name_sort_check}
      ? $config->{name_sort_check}
      : $name_sort_input;

    if ( $self->{verbose} ) {
        print STDERR "======== &prelim_filter: Start ========\n";
    }
    my ( $fn, $path, $suffix ) = fileparse( $input, ( '.srt.bam', '.bam' ) );
    my $header       = $self->_run_cmd("samtools view -H $input");
    my @split_header = split( /\n/, $header );
    my $last_pg      = $split_header[-1];
    my $last_pg_id   = ( split /\t/, $last_pg )[1];
    $last_pg_id =~ s/ID\:/PP\:/;
    $header = $header
      . "\@CO\tID:Prelim-filter\tPN:LGTseek\t$last_pg_id\tVN:$VERSION\n";
    my @output_list;    ## Array of output bams to return

    ## Check if the output exists already
    ## $output isn't set right at this point in the scrip to check if the final output already exist. FIX later.
    my $files_exist = 0;
    if ( -e "$output_dir$fn\_0_prelim.bam" ) { $files_exist = 1; }
    if ( $files_exist == 1 && $overwrite == 0 ) {
        if ( $self->{verbose} ) {
            print STDERR
"*** Warning *** Already found the output for &prelim_filter: $output_dir$fn\_0_prelim.bam\n";
        }
        chomp( my @output_list = `find $output_dir -name '$fn\*_prelim.bam'` );
        return \@output_list;
    }

    ## Check we dont' have the name-sorted.bam already
    my $sorted_bam = "$output_dir$fn\_name-sorted.bam";
    if ( -e "$sorted_bam" ) { $files_exist = 1; }
    if ( $files_exist == 1 && $overwrite == 0 ) {
        if ( $self->{verbose} ) {
            print STDERR
"*** Warning *** Already found &prelim_filter \$sorted_bam, starting to filter: $sorted_bam.\n";
        }
        $name_sort_input = 0;
    }

    ## (1). Sort the bam by name instead of position
    if ( $name_sort_input == 1 ) {
        if ( $self->{verbose} ) {
            print STDERR "======== &prelim_filter: Sort Start ========\n";
        }
        my $cmd = "samtools sort -n  -o $output_dir$fn\_name-sorted\.bam $input";
        if ( $self->{verbose} ) {
            print STDERR "======== &prelim_filter Sort: $cmd ===\n";
        }
        $self->_run_cmd("$cmd");
        if ( $self->{verbose} ) {
            print STDERR "======== &prelim_filter: Sort Finished ========\n";
        }
    }
    else {
        $sorted_bam = $input;
    }

    ## (2). QUICK/dirty Check to make sure the reads are PE and name sorted properly.
    if ( defined $self->{name_sort_check} and $self->{name_sort_check} == 1 ) {
        if ( $self->{verbose} ) {
            print STDERR
"========= &prelim_filter 2x Check PE & name sort: Start ========\n";
        }
        open( my $chk, "samtools view $sorted_bam | head |" )
          or $self->fail(
"*** Error *** &prelim_filter can't open the name-sorted-input.bam: $sorted_bam\n"
          );
        my %dbl_chk_hash;
        while (<$chk>) {
            chomp;
            my @line = split;
            my $read = $line[0];
            $read =~ s/(\/\d{1})$//;
            $dbl_chk_hash{$read}++;
        }
        close $chk;

        my $die_count = 0;
        foreach my $reads ( keys %dbl_chk_hash ) {
            if ( $dbl_chk_hash{$reads} != 2 ) {
                $die_count++;
                print STDERR
                  "*** Warning *** This read is missing the mate: $reads \n";
            }
        }

        if ( $die_count >= 2 ) {
            die
"*** Error *** This bam is missing atleast 2 proper-pairs in the first 10 reads. It is highly likely this bam isn't PE or name sorted properly.
            Please manually inspect this bam: $sorted_bam . If the bam is correct, relaunch prelim_filter with --name_sort_check=0.\n";
        }
    }
    if ( $self->{verbose} ) {
        print STDERR
"========= &prelim_filter 2x Check PE & name sort:: Finished ========\n";
    }

    ## (3). Prelim filtering.
    my $more_lines        = 1;
    my $num_pass          = 0;
    my $num_null          = 0;
    my $num_singletons    = 0;
    my $num_secondary     = 0;
    my $num_supplementary = 0;
    my $MM                = 0;
    my $MU                = 0;
    my $UU                = 0;
    my $SC                = 0;

    if ( $prelim_filter == 1 ) {
        if ( $self->{verbose} ) {
            print STDERR "========= &prelim_filter: Filter Start ========\n";
        }
        ## Setup and open output
        my $i = 0;   ## Used to count number of lines per bam as being processed
        my $count = 0;    ## Used to count number of split output bams
        my $output = "$output_dir/$fn\_$count\_prelim.bam";
        open( my $out, "| samtools view -S - -bo $output" )
          || $self->fail("*** Error *** Can't open: $output because: $!\n");
        print $out "$header";

        # print $out "\@CO\tID:Prelim-filter\tPG:LGTseek\tVN:$VERSION\n";
        push( @output_list, $output );

        ## Open input bam
        open( my $infh, "samtools view $sorted_bam |" )
          || $self->fail(
"*** Error *** &prelim_filter can't open the \$sorted_bam: $sorted_bam because: $1\n"
          );
        ## (2). Actual Prelim Filtering starting:
        ##      Read through and filter bam reads, splitting output as we go

        while ( $more_lines == 1 ) {
            my $read1 = $self->_read_bam_line($infh);
            my $read2 = $self->_read_bam_line($infh);

            my $print = 0;
            ## (3). Split output: close the current output and open a new output
            if ( $i >= $seqs_per_file && $split_bam == 1 ) {
                close $out
                  or $self->fail(
"*** Error *** Unable to close &prelim_filter output: $output\n"
                  );
                $count += 1;
                $output = "$output_dir/$fn\_$count\_prelim.bam";
                open( $out, "| samtools view -S - -bo $output" )
                  || $self->fail(
                    "*** Error *** Can't open: $output because: $!\n");
                push( @output_list, $output );
                print $out "$header";
                $i = 0;
            }

            ## Skip "null" reads
            while ( $read2->{id} =~ /null/ ) {
                $num_null++;
                $read2 = $self->_read_bam_line($infh);
            }
            while ( $read1->{id} =~ /null/ ) {
                $num_null++;
                $read1 = $read2;
                $read2 = $self->_read_bam_line($infh);
            }
            ## Skip Secondary & Supplementary reads
            while ($read1->{flag}->{secondary}
                or $read1->{flag}->{supplementary} )
            {
                if ( $read1->{flag}->{secondary} )     { $num_secondary++; }
                if ( $read1->{flag}->{supplementary} ) { $num_supplementary++; }
                $read1 = $read2;
                $read2 = $self->_read_bam_line($infh);
            }
            while ($read2->{flag}->{secondary}
                or $read2->{flag}->{supplementary} )
            {
                if ( $read2->{flag}->{secondary} )     { $num_secondary++; }
                if ( $read2->{flag}->{supplementary} ) { $num_supplementary++; }
                $read2 = $self->_read_bam_line($infh);
            }

            ## Make sure we have the same read ID's
            while ( $read1->{id} ne $read2->{id} ) {
                $num_singletons++;
                $read1 = $read2;
                $read2 = $self->_read_bam_line($infh);
            }

            ## Count number of read types & filter for potential LGT reads (MU, UM, UU)
            if ( !$read1->{flag}->{qunmapped} && !$read2->{flag}->{qunmapped} )
            {
                $MM += 2;
            }    ## Count M_M
            if ( $read1->{flag}->{qunmapped} or $read2->{flag}->{qunmapped} ) {
                $print = 1;
                if ( $read1->{flag}->{qunmapped}
                    && !$read2->{flag}->{qunmapped} )
                {
                    $MU += 2;
                }    ## U_M
                if (  !$read1->{flag}->{qunmapped}
                    && $read2->{flag}->{qunmapped} )
                {
                    $MU += 2;
                }    ## M_U
                if (   $read1->{flag}->{qunmapped}
                    && $read2->{flag}->{qunmapped} )
                {
                    $UU += 2;
                }    ## U_U
            }

            ## FILTER FOR soft clipped reads
            if ( $keep_softclip == 1 ) {
                map {
                    if ( $_ =~ /(\d+)S/ && $1 >= $softclip_min ) {
                        $print = 1;
                        $SC += 2;
                    }    ## 10.21.14 Testing more liberal softclip parsing
                } ( $read1->{cigar}, $read2->{cigar} );
            }

            if ( $print == 1 ) {
                ## Fix flags if one of the reads is both first & second
                ( $read1, $read2 ) = $self->_fix_flags( $read1, $read2 );

                ## Create filtered/fixed bam lines
                map {
                    chomp( my $original_line = $_->{line} );
                    my @split_original_line = split( /\t/, $original_line );
                    my $fixed_flag = $self->_create_flag( $_->{flag} );

# Create the polished output (removed \1 & \2 from id, removed flag=First & Second, removed tags)
                    my $new_line = join( "\t",
                        $_->{id}, $fixed_flag,
                        @split_original_line[ 2 .. 10 ] );
                    print $out "$new_line\n";
                } ( $read1, $read2 );

                $i        += 2;
                $num_pass += 2;
            }
        }
        close $out;

        ## Check to make sure the last file isn't empty.
        if ( $self->empty_chk( { input => $output } ) == 1 ) {
            $self->_run_cmd("rm $output");
            my $remove_empty_bam = pop(@output_list);
        }
        if ( $name_sort_input == 1 && $sorted_bam =~ /name-sorted.bam$/ ) {
            $self->_run_cmd("rm $sorted_bam");
        }
        if ( $self->{verbose} ) {
            print STDERR "========= &prelim_filter: Filter Finished ========\n";
        }
    }
    else {
        push( @output_list, $sorted_bam );
    }

    my @sort_out_list = sort @output_list;
    if ( $self->{verbose} ) {
        print STDERR
"========= &prelim_filter: MM:$MM | Bad Singletons:$num_singletons | Null:$num_null | Secondary:$num_secondary | Supplementary:$num_supplementary ==\n";
        print STDERR
"========= &prelim_filter: Pass:$num_pass | MU:$MU | UU:$UU | SC:$SC ==\n";
        print STDERR "======== &prelim_filter: Finished ========\n";
    }
    return \@sort_out_list;
}

=head2 &_read_bam_line

 Title   : _read_bam_line
 Usage   : my $read = $self->_read_bam_line($filehandle)
 Function: Read in a bam line, split important data into hash for each read, convert flag
 Returns : hash->{id}           = id w/ \1 or \2 substituted off
           hash->{flag}         = converted flag
           hash->{sequence}
           hash->{cigar}
 Args    :
=cut

sub _read_bam_line {
    my ( $self, $fh ) = @_;
    my $line = <$fh>;
    if ( !$line ) { last; }
    chomp($line);
    my ( $id, $flag, $cigar, $sequence ) = ( split /\t/, $line )[ 0, 1, 5, 9 ];
    $id =~ s/(.+)\/\d+/$1/;
    my $converted_flag = parse_flag($flag);
    return {
        id       => $id,
        flag     => $converted_flag,
        sequence => $sequence,
        cigar    => $cigar,
        line     => $line,
    };
}

=head2 &_fix_flags

 Title   : _fix_flags
 Usage   : ($read1,$read2)=$self->_fix_flags($read1,$read2);
 Function: Fix reads with first & second flags
 Returns :
 Args    :
=cut

sub _fix_flags {
    my ( $self, $read1, $read2 ) = @_;

    # If read1 is First and Second ...
    if ( $read1->{flag}->{first} && $read1->{flag}->{last} ) {
        print STDERR "TEST1\n";

        # If read2 is last, make read1 first
        if ( $read2->{flag}->{last} ) {
            $read1->{flag}->{last} = 0;
        }
        elsif ( $read2->{flag}->{first} ) {
            $read1->{flag}->{first} = 0;
        }
    }

    # If read2 is First and Second ...
    if ( $read2->{flag}->{first} && $read2->{flag}->{last} ) {
        print STDERR "TEST2\n";

        # If read1 is first, make read2 last
        if ( $read1->{flag}->{first} ) {
            $read2->{flag}->{first} = 0;
        }
        elsif ( $read1->{flag}->{last} ) {
            $read2->{flag}->{last} = 0;
        }
    }
    return ( $read1, $read2 );

}

=head2 &filter_bam_by_ids

 Title   : filter_bam_by_ids
 Usage   : my $filtered_bam = $lgtseek->filter_bam_by_ids({input => <bam>, good_list=> <list of desired ids>})
           output = output_dir/prefix_suffix.bam
 Function: Filter bam by ids
 Returns : A hash{bam} = new filtered bam
           A hash{count} = # ids found
 Args    : A hash containing potentially several config options:
        input_bam           => bam for filtering
        output_dir          => directory for output
        output_prefix       => Prefix for output
        output_suffix       => Suffix for output
        header_comment      => "\@PG\tID:lgtseeq\tPG:filter-criteria\tVN:$VERSION"
        output              => Full Path, prefix, name, and suffix for output.
        good_list           => File path with list of desired reads
        bad_list            => File path with a list of reads to discard
=cut

sub filter_bam_by_ids {
    my ( $self, $config ) = @_;
    my $samtools = $self->{samtools_bin};
    if ( $self->{verbose} ) {
        print STDERR "======== &filter_bam_by_ids: Start ========\n";
    }
    if ( !$config->{input_bam} ) {
        $self->fail(
"*** Error *** Must pass &filter_bam_by_ids an input_bam => <BAM_TO_FILTER>\n"
        );
    }
    if ( $self->empty_chk( { input => "$config->{input_bam}" } ) == 1 ) {
        print STDERR
"*** Warning ***: Can not filter ids from an empty input bam: $config->{input_bam}\n";
        return { count => 0, file => $config->{input_bam} };
    }
    if ( !$config->{good_list} && !$config->{bad_list} ) {
        $self->fail(
"*** Error *** Must pass &filter_bam_by_ids a file with a list of reads to filter on. Use good_list => || bad_list => \n"
        );
    }
    ## Setup hash of ids
    my $good_ids = {};
    my $bad_ids  = {};
    my %found_ids;
    if ( $config->{good_list} ) {
        $good_ids = $self->_read_ids( { list => $config->{good_list} } );
    }
    elsif ( $config->{bad_list} ) {
        $bad_ids = $self->_read_ids( { list => $config->{bad_list} } );
    }
    ## Setup input and output bams
    my $input = $config->{input_bam};
    my ( $fn, $path, $suf ) = fileparse( $input, ".bam" );
    my $out_dir = defined $config->{output_dir} ? $config->{output_dir} : $path;
    my $prefix =
      defined $config->{output_prefix} ? $config->{output_prefix} : $fn;
    my $suffix =
      defined $config->{output_suffix} ? $config->{output_suffix} : "filtered";
    my $out =
      defined $config->{output}
      ? "$config->{output}"
      : "$out_dir/$prefix\_$suffix\.bam";
    my $header = $self->_run_cmd("$samtools view -H $input");
    open( my $in, "-|", "$samtools view $input" )
      or $self->fail(
"*** Error *** &filter_bam_by_ids can't open input bam: $input because: $!\n"
      );
    open( my $fh, "| $samtools view -S - -bo $out" )
      or $self->fail(
"*** Error *** &filter_bam_by_ids can't open  output bam: $out because: $!\n"
      );
    print $fh "$header";

    if ( defined $config->{header_comment} ) {
        my @split_header = split( /\n/, $header );
        my $last_pg      = $split_header[-1];
        my $last_pg_id   = ( split /\t/, $last_pg )[1];
        $last_pg_id =~ s/ID\:/PP\:/;
        chomp( $config->{header_comment} );
        print $fh "$config->{header_comment}\t$last_pg_id\n";
    }

    while (<$in>) {
        chomp;
        my @fields = split(/\t/);
        if ( $config->{good_list} && $good_ids->{ $fields[0] } ) {
            print $fh "$_\n";
            $found_ids{ $fields[0] }++;
        }
        if ( $config->{bad_list} && !$bad_ids->{ $fields[0] } ) {
            print $fh "$_\n";
            $found_ids{ $fields[0] }++;
        }
    }
    close $in;
    close $fh;
    my $count = 0;
    foreach my $keys ( keys %found_ids ) { $count++; }
    if ( $self->{verbose} ) {
        print STDERR "======== &filter_bam_by_ids: Finished ========\n";
    }
    return {
        count => $count,
        bam   => $out
    };
}

=head2 &_read_ids

 Title   : _read_ids
 Usage   : my $id_hash = $self->_read_ids({list => <ID_LIST_FILE>})
 Function: Create a hash of id's from a list file
 Returns : A hash of ids
 Args    : A hash containing potentially several config options:
        list     => File with ids
=cut

sub _read_ids {
    my ( $self, $config ) = @_;
    if ( !$config->{list} ) {
        $self->fail("*** Error *** Must pass &_read_ids a list =>\n");
    }
    if ( $self->empty_chk( { input => $config->{list} } ) == 1 ) {
        carp "*** Warning ***: &_read_ids input: $config->{list} is empty\n";
    }
    my $file = $config->{list};
    my %hash;
    open IN, "<$file"
      or $self->fail(
        "*** Error *** &_read_ids unable to open: $file because: $!\n");
    while (<IN>) {
        chomp;
        $hash{$_} = 1;
    }
    close IN;
    return \%hash;
}

=head2 &empty_chk

 Title   : empty_chk
 Usage   : last if($self->empty_chk({input => <FILE>}) == 1);
 Function: Check to make sure a file is not empty
 Returns : 1 = Empty; 0 = input is not empty
 Args    : A hash containing potentially several config options:
        input     => Check if input is empty
=cut

sub empty_chk {
    my ( $self, $config ) = @_;
    if ( !$config->{input} ) {
        $self->fail("*** Error *** Must pass &empty_chk a file.");
    }
    my $file  = $config->{input};
    my $empty = 0;                  ## 0 = False, 1=True, file is empty.
    my $count;
	$self->{samtools_bin} = '/usr/bin/samtools' unless (defined $self->{samtools_bin});
    if ( $file =~ /\.bam$/ ) {
        $count = `$self->{samtools_bin} view $file | head | wc -l`;
    }
    elsif ( $file =~ /\.gz/ ) {
        $count = `zcat $file | head | wc -l`;
    }
    else {
        $count = `head $file | wc -l`;
    }
    if ($?) {
        $self->fail(
            "*** Error *** &empty_chk could not check: $config->{input}\n");
    }
    if ( $count == 0 ) { $empty = 1; }
    return $empty;
}

=head2 &fail

 Title   : fail
 Usage   : if(1+1!=2){$self->fail("Because you can't add!");}
 Function: Die with either die or confess depending on verbosity.
 Args    : Die message.
=cut

sub fail {
    my ( $self, $message ) = @_;
    chomp($message);
    if ( $self->{verbose} )      { confess "$message\n"; }
    if ( $self->{verbose} != 1 ) { die "$message\n"; }
}

=head2 &validated_bam

 Title   : validated_bam
 Usage   : my $validated_bam = lgtseek->validated_bam({by_clone=> <input> });
 Function: Creates a bam based on LGTFinder Results for reads with Human/Euk & Bacteria LGT
 Returns:  Hash->{count} && Hash->{file} bam based on LGTFinder output.
 Args    : A hash containing potentially several config options:
        input           =>  Input bam to parse LGT reads from
        by_clone        =>  LGTFinder Output
        by_trace        =>  Not implemented yet.
        output_dir      =>  /path/for/output/
        output_prefix   =>  {$prefix}_$suffix.bam  [ {$input_fn}.bam ]
        output_suffix   =>  $prefix_{$suffix}.bam  [ "filtered" ]
        output          =>  /path/and/name/for/output.bam
=cut

sub validated_bam {
    my ( $self, $config ) = @_;
    if ( $self->{verbose} ) {
        print STDERR "======== &validated_bam: starting ========\n";
    }
    if ( !$config->{input} ) {
        $self->fail(
"*** Error *** Must pass &validated_bam an input bam to parse reads from with input =>.\n"
        );
    }
    if ( !$config->{by_clone} && !$config->{by_trace} ) {
        $self->fail(
"*** Error *** Must pass &validated_bam an LGTFinder output with by_clone => or by_trace=>.\n"
        );
    }
    if ( $self->empty_chk( { input => $config->{input} } ) == 1 ) {
        if ( $self->{verbose} ) {
            print STDERR
"========= &validated_bam: Skipping this input because it is empty! ========\n";
        }
        return {
            count => "0",
            file  => undef,
        };
    }
    my $input = "$config->{input}";
    my ( $fn, $path, $suf ) = fileparse( $input, ".bam" );
    my $out_dir = defined $config->{output_dir} ? $config->{output_dir} : $path;
    my $prefix =
      defined $config->{output_prefix} ? $config->{output_prefix} : $fn;
    my $suffix =
      defined $config->{output_suffix} ? $config->{output_suffix} : "filtered";
    my $out =
      defined $config->{output}
      ? "$config->{output}"
      : "$out_dir/$prefix\_$suffix\.bam";
    if ( $config->{output} ) {
        my ( $foo, $bar, $salad ) =
          fileparse( $config->{output}, qr/\.[^\.]+/ );
        $out_dir = $bar;
    }

	#my $host_blast_otu =
	#	defined $config->{host_blast_otu} ? $config->{host_blast_otu} : $self->{host_blast_otu};    ## Example: Homo
	#my $host_lin       =
	#	defined $config->{host_lineage}   ? $config->{host_lineage}   : $self->{host_lineage};    ## Example: Eukaryota
	#my $donor_lin      =
	#	defined $config->{donor_lineage}  ? $config->{donor_lineage}  : $self->{donor_lineage};     ## Example: Bacteria

    open( IN, "<", "$config->{by_clone}" )
      || $self->fail(
        "*** Error *** &validated_bam can not open input: $config->{by_clone}\n"
      );
    open( OUT, ">", "$out_dir/$prefix\_valid_blast_read.list" )
      || $self->fail(
"*** Error *** &validated_bam can not open output: $out_dir/$prefix\_valid_blast.list\n"
      );
    ##  Make a list of reads from lgt_finder that have a "valid blast."
    while (<IN>) {
        chomp;
        my ( $read, $otu1, $otu2, $lca1, $lca2 ) = ( split /\t/, $_ )[ 0, 1, 2, 6, 11 ];
			next if ( $lca1 =~ /Eukaryota|Homo/ && $lca2 =~ /Eukaryota|Homo/ );
		    next if ( $lca1 =~ /Bacteria/  && $lca2 =~ /Bacteria/ );
		#next if ( $otu1 =~ /$host_blast_otu/ && $otu2 =~ /$host_blast_otu/ );    ## Example: /Homo/
		#next if ( $otu1 !~ /$host_blast_otu/ && $otu2 !~ /$host_blast_otu/ );
		#next if ( $lca1 =~ /$host_lin/       && $lca2 =~ /$host_lin/ );          ## Example: /Eukaryota/
		#next if ( $lca1 =~ /$donor_lin/      && $lca2 =~ /$donor_lin/ );        ## Example: /Bacteria/
        print OUT "$read\n";
    }
    close IN
      || $self->fail(
"*** Error *** &validated_bam can't close input by_clone: $config->{by_clone} because: $!\n"
      );
    close OUT
      || $self->fail(
"*** Error *** &validated_bam can't close output valid_blast.list: $out_dir/$prefix\_valid_blast_read.list because: $!\n"
      );

    ## Use the list of blast validated LGT's create a new bam.
    ## $valid_bam is a hash->{count} & ->{file}
    my $valid_bam = $self->filter_bam_by_ids(
        {
            input_bam => $input,
            header_comment =>
              "\@CO\tID:Blast-validated\tPG:LGTseek\tVN:$VERSION",
            good_list => "$out_dir/$prefix\_valid_blast_read.list",
            output    => $out
        }
    );
    return $valid_bam;
}

=head2 downloadEGA

 Title   : downloadEGA
 Usage   : my $bam_file = $lgtseek->downloadEGA(
            {
                input           =>      'ID#_to_search_for' or 'ID#_1_1.rnaseq.fastq.gz.gpg,ID#_1_2.rnaseq.fastq.gz.gpg'
                output_dir      =>      '/file/path/for/output/dir/',
                aln_human       =>      <0|1> [0] 1= Aln fastq downloads to hg19
                ega_DL_client   =>      '/path/to/EGAs/webin-data-streamer-EGA-Download-Client.jar'
                overwrite       =>      < 0|1 > [ .lgtseek.conf ]
            });
 Function:  Download fastq files from EGA, decrypt, and map to human.
            ## This subroutine does NOT support BAM downloads yet##

 Returns :  A bam file mapped at the human genome.
 Args    :

=cut

sub downloadEGA {
    my ( $self, $config ) = @_;
    if ( $self->{verbose} ) {
        print STDERR "======= &downloadEGA: Start ========\n";
    }

    if ( !$config->{input} ) {
        $self->fail("*** Error *** &downloadEGA must have an input.\n");
    }

    # Check for all the genetorrent related options
    my $output_dir =
      defined $config->{output_dir}
      ? $config->{output_dir}
      : $self->{output_dir};
    $self->_run_cmd("mkdir -p $output_dir");
    my $ega_DL_client =
      defined $config->{ega_DL_client}
      ? $config->{ega_DL_client}
      : $self->{ega_DL_client};
    my $java = defined $self->{java_bin} ? $self->{java_bin} : "java";
    my $overwrite =
      defined $config->{overwrite} ? $config->{overwrite} : $self->{overwrite};
    my $max_retry_attempts =
      defined $config->{retry_attempts}
      ? $config->{retry_attempts}
      : $self->{retry_attempts};    ## Defaults to lgtseq.conf

# First, determine if the input is an ID to search for or a specific file(s) to download.
    my @files_to_download;
    my @fastqs;
    ## If the input looks like a filename:
    if ( $config->{input} =~ /(_[12]{1}){0,1}(\.\w+)?\.f\w{0,3}q(.gz)?(.gpg)?/ )
    {
        my @split_filenames = split( /,/, $config->{input} );
        if ( scalar(@split_filenames) > 2 ) {
            $self->fail(
"*** Error *** The input looks like more than 2 fastq files to download. Only pass 2 fastq per input for download and mapping.\n"
            );
        }
        foreach my $files (@split_filenames) {
            $files =~ /(.+)\.gpg/;
            push( @fastqs, $1 );
            if (
                ( !-e "$output_dir/$1.cip" and !-e "$output_dir/$1" )
                or ( ( -e "$output_dir/$1.cip" or -e "$output_dir/$1" )
                    and $overwrite == 1 )
              )
            {
                push( @files_to_download, $files );
            }
        }
    }
    else {
        if ( $self->{verbose} ) {
            print STDERR
              "========= &downloadEGA: Query FASTQs to download ========\n";
        }
        my $retry_attempts_made = 0;
        my $retry               = 1;
      QUERY_EGA:
        while ( $retry == 1 && ( $retry_attempts_made <= $max_retry_attempts ) )
        {
            my @ega_query = split(
                /\n/,
                $self->_run_cmd(
"$java -jar $ega_DL_client -p -list -user $self->{ega_user} -pass $self->{ega_pass} | grep $config->{input} | cut -f1"
                )
            );
            ## Check we got something from the query. If not, retry!
            if ( $ega_query[0] !~ /\w+/ ) {
                $retry_attempts_made++;
                sleep 90;
                goto QUERY_EGA;
            }
            else { $retry = 0; }

            # Add each file from the ega_query to the list of files to download
            map {
                my $fq_fn = fileparse( $_, ".gpg" );
                push( @fastqs, $fq_fn );
                if (
                    (
                            !-e "$output_dir/$fq_fn.cip"
                        and !-e "$output_dir/$fq_fn"
                    )
                    or (
                        (
                               -e "$output_dir/$fq_fn.cip"
                            or -e "$output_dir/$fq_fn"
                        )
                        and $overwrite == 1
                    )
                  )
                {
                    push( @files_to_download, $_ );
                }
                elsif ( $self->{verbose} ) {
                    print STDERR
"========= &downloadEGA: Skipping download, $fq_fn exists ========\n";
                }
            } @ega_query;
        }
    }

    # Second, Download the files
    ## Change to the appropriate directory for downloading
    my $current_directory =
      getcwd;  ## So we can change back after download and decrypt are complete.
    chdir($output_dir);
    ## Build the -files string for the download client.
    my $ega_DL_client_files_to_DL_string = " -f ";
    $ega_DL_client_files_to_DL_string = join( ",", @files_to_download );
    $ega_DL_client_files_to_DL_string =
      "-f " . $ega_DL_client_files_to_DL_string;
    ## Now, download the file(s)
    if ( scalar( @files_to_download >= 1 ) ) {
        if ( $self->{verbose} ) {
            print STDERR
              "========= &downloadEGA: Downloading FASTQs ========\n";
        }
        my $retry               = 1;
        my $retry_attempts_made = 0;
      DOWNLOAD_EGA:
        while ( $retry == 1 && ( $retry_attempts_made <= $max_retry_attempts ) )
        {
            my $download =
"$java -jar $ega_DL_client -p -user $self->{ega_user} -pass $self->{ega_pass} -d -re LGTSeq -destination $output_dir $ega_DL_client_files_to_DL_string > $output_dir/EGA_download_log.txt";
            if ( $self->{verbose} ) { print STDERR "CMD: $download\n"; }
            `$download`;
            ## If we can detect there was a problem, either retry or exit.
            if (   $?
                or $self->_run_cmd("tail -n 1 $output_dir/EGA_download_log.txt")
                =~ /RETURN CODE 1/ )
            {
                if ( $retry_attempts_made == $max_retry_attempts ) {
                    $self->fail(
                        "======== &downloadEGA: Multiple download errors.\n");
                }
                else {
                    print STDERR
                      "======== &downloadEGA: Download failed, trying again.";
                    foreach my $files (@files_to_download) {
                        $self->_run_cmd("rm $output_dir/$files*");
                    }
                    $self->_run_cmd("rm $output_dir/EGA_download_log.txt");
                    $retry_attempts_made++;
                    goto DOWNLOAD_EGA;
                }
            }
            ## Otherwise report a successful download and move on
            else {
                if ( $self->{verbose} ) {
                    print STDERR
                      "========= &downloadEGA: FASTQs downloaded ========\n";
                }
                $retry = 0;
            }
        }
    }

    # Third, Decrypt the downloaded files
    my @tmp_return_list;
    foreach my $fastq (@fastqs) {
        if ( -e "$fastq.cip" and ( !-e $fastq or $overwrite == 1 ) ) {
            if ( $self->{verbose} ) {
                print STDERR
                  "========= &downloadEGA: Decrypting FASTQ: $fastq ========\n";
            }
            $self->_run_cmd(
"$java -jar $ega_DL_client -p -dc -re LGTSeq -destination $output_dir -f $fastq.cip"
            );
        }
        elsif ( $self->{verbose} ) {
            print STDERR
"========= &downloadEGA: Skipping decrypt, $fastq already decrypted ========\n";
        }
        my $tmp_fastq = "$output_dir/$fastq";
        $tmp_fastq =~ s/\/{2,}/\//g;
        push( @tmp_return_list, $tmp_fastq );
    }
    ## Move back to the original directory.
    chdir($current_directory);
    if ( -e "$output_dir/log.txt" ) {
        $self->_run_cmd("rm $output_dir/log.txt");
    }

    if ( $self->{verbose} ) {
        print STDERR "======= &downloadEGA: Finished ========\n";
    }
    return \@tmp_return_list;
}

=head2	print_tab

Title	:	print_tab
Usage	:	LGT:LGTSeek->print_tab($filename, $header_array_ref, $counts_array_ref);
Function	:	Prints out counts based on various means of filtering the BAM reads
Returns	:	File with tab-delimited counts

=cut

sub print_tab {
    my ( $class, $file, $header, $vals ) = @_;
    open OUT, ">$file" or confess("Couldn't open $file for writing\n");
    print OUT join( "\t", @$header );
    print OUT "\n";
    print OUT join( "\t", @$vals );
    print OUT "\n";
    close OUT;
    return;
}

1;

__END__
