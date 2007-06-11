#!/usr/local/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Ergatis::Logger;
use File::Find;
use XML::Twig;
use Data::Dumper;

######################### GLOBALS AND CONSTANTS ####################################
my $input_file;
my @btab_files = ();
my @aragorn_files = ();
my @trna_files = ();
my $info_file;
my $tmp_dir = "./";
my $other_opts;
my $output_file;
my $polypeptide_multi_fsa;
my @hmms;
my @hmm_frags;

my $PHAGE_FINDER_DIR;
######################### OPTIONS AND LOGGER########################################
my %options = ();
my $results = GetOptions (\%options, 
                          'input_bsml_file|i=s',
                          'input_polypeptide_multi_fsa|p=s',
                          'btab_list|b=s',
                          'aragorn_list|a=s',
                          'trna_scan_list|t=s',
                          'info_file|f=s',
                          'other_opts|o=s',
                          'tmp_dir|m=s',
                          'output_file|u=s',
                          'phage_finder_dir|g=s',
                          'log|l=s',
                          'debug=s',
                          'help') || &pod;

#Setup the logger
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

&check_parameters( \%options );

my $PHAGE_FINDER_BIN = "$PHAGE_FINDER_DIR/bin/phage_finder_ergatis.pl";
my $HMM_SEARCH_EXE = "/usr/local/bin/hmmsearch";
my $PHAGE_HMM = "$PHAGE_FINDER_DIR/PHAGE_HMMs_dir";
my $PHAGE_HMM_FRAGS = "$PHAGE_FINDER_DIR/PHAGE_FRAG_HMMS_dir";

##################################### MAIN ##############################################

#Get the identifier from the defline of the input fsa files
my $identifier;
my $assembly_fsa;

my $twig = new XML::Twig( 'twig_handlers' => {
    'Sequence[@class="assembly"]' => sub {
        $identifier = $_[1]->att('id');
        my $sdi = $_[1]->first_child('Seq-data-import');
        $assembly_fsa = $sdi->att('source');
    }
});

$twig->parsefile( $input_file );
$logger->logdie("Unable to parse assembly identifier out of $input_file") unless( $identifier );
$logger->logdie("Unable to find assembly fasta from $input_file.  Check the assembly Seq-data-import element")
    unless( $assembly_fsa && -e $assembly_fsa );

#Find the files related to this sequence.
my $aragorn_file = &get_related_files( $identifier, @aragorn_files ) if( @aragorn_files > 0 );
my $trna_file = &get_related_files( $identifier, @trna_files ) if( @trna_files > 0 );

#Get a list of the hmm files to search with 
find( \&phage_hmm, $PHAGE_HMM);
find( \&phage_hmm_frags, $PHAGE_HMM_FRAGS );

my $combined_hmm = "$tmp_dir/$identifier.GLOCAL.combined.htab";
if( -e $combined_hmm ) {
    $logger->warn("removing $combined_hmm to make room for new");
    system( "rm $combined_hmm" );
}

my $hmm_count = scalar( @hmms );

#do all the hmms
print "\n\nRunning HMMs\n";
my $count = 0;
foreach my $hmm ( @hmms ) {
    $count++;
    my $hmm_cmd = "$HMM_SEARCH_EXE $hmm $polypeptide_multi_fsa >> $combined_hmm";
    system( $hmm_cmd );
    print "\r$count/$hmm_count";
}
print "\nDone.\n";

my $frag_count = scalar( @hmm_frags );
$count = 0;
my $combined_frag_hmm = "$tmp_dir/$identifier.FRAG.combined.htab";
if( -e $combined_frag_hmm ) {
    $logger->warn( "Removing $combined_frag_hmm" );
    system( "rm $combined_frag_hmm" );
}

foreach my $hmm_frag ( @hmm_frags ) {
    $count++;
    my $frag_hmm_cmd = "$HMM_SEARCH_EXE $hmm_frag $polypeptide_multi_fsa >> $combined_frag_hmm";
    system( $frag_hmm_cmd );
    print "\r$count/$frag_count";
}
print "\nDone with hmms.\n";
print "\n\nStarting to combine blast hits\n";

#Concat all the blast outputs into one file
my $blast_out = "$tmp_dir/$identifier.blast.combined.btab";
system( "rm -f $blast_out") if( -e $blast_out );

my $blast_count = scalar( @btab_files );
$count = 0;
foreach my $blast_file ( @btab_files ) {
    $count++;
    system("cat $blast_file >> $blast_out");
    print "\r$count/$blast_count";
}
print "\nDone combining blast outputs\n";


 #Generate the command line string
my $cmd = "$PHAGE_FINDER_BIN -i $info_file -t $blast_out -g $combined_hmm -f $combined_frag_hmm -b $tmp_dir ";
$cmd .= "-r $trna_file " if( $trna_file );
$cmd .= "-n $aragorn_file " if( $aragorn_file );
$cmd .= "-A $assembly_fsa " if( $assembly_fsa );
$cmd .= $other_opts if( $other_opts );
$cmd .= " > $output_file ";
print "\n$cmd\n";
system( $cmd );

if ($? == -1) {
	print "failed to execute: $!\n";
}
elsif ($? & 127) {
	printf "child died with signal %d, %s coredump\n",
    ($? & 127),  ($? & 128) ? 'with' : 'without';
}
else {
	printf "child exited with value %d\n", $? >> 8;
}


####################################### SUB ROUTINES #####################################
sub phage_hmm {
    return unless( $_ =~ /\.HMM$/ );
    push( @hmms, $File::Find::name );
}

sub phage_hmm_frags {
    return unless( $_ =~ /\.FRAG$/ );
    push( @hmm_frags, $File::Find::name );

}

sub get_related_files {
    my ($ident, @files) = @_;
    my @retval;

    foreach my $file ( @files ) {
        push( @retval, $file ) if( $file =~ /$ident/ );
    }

    return $retval[0];   

}

sub check_parameters {
    my $opts = shift;

    #input_file option is required
    if( $opts->{'input_bsml_file'} ) {
        $logger->logdie("Option input_bsml_file [$opts->{'input_bsml_file'}] does not exist")
            unless( -e $opts->{'input_bsml_file'} );
        $input_file = $opts->{'input_bsml_file'};
    } else {
        $logger->logdie("Option input_bsml_file is required");
    }

    #btab_list is required
    if( $opts->{'btab_list'} ) {
        &get_files_from_list( $opts->{'btab_list'}, \@btab_files );
    } else {
        $logger->logdie("Option btab_list is required");
    }

    #aragorn list is optional
    if( $opts->{'aragorn_list'} ) {
        &get_files_from_list( $opts->{'aragorn_list'}, \@aragorn_files );
    }

    #tRNAscan list is optional
    if( $opts->{'trna_scan_list'} ) {
        &get_files_from_list( $opts->{'trna_scan_list'}, \@trna_files );
    }

    #info_file is required
    if( $opts->{'info_file'} ) {
        $logger->logdie("info file [$opts->{'info_file'}] does not exist") 
            unless( -e $opts->{'info_file'} );
        $info_file = $opts->{'info_file'};
    } else {
        $logger->logdie("Option info_file is required");
    }

    #other opts is optional
    if( $opts->{'other_opts'} ) {
        $other_opts = $opts->{'other_opts'};
    }

    #tmp_dir is required
    if( $opts->{'tmp_dir'} ) {
        $tmp_dir = $opts->{'tmp_dir'};
    }

    #output_file is required
    unless( $opts->{'output_file'} ) {
        $logger->logdie("Option output_file is required");
    } else {
        $output_file = $opts->{'output_file'};
    }

    #input_polyeptide_fsa_list is required
    unless( $opts->{'input_polypeptide_multi_fsa'} ) {
        $logger->logdie("Option input_polypeptide_multi_fsa is required");
    } else {
        $logger->logdie("$opts->{'input_polypeptide_multi_fsa'} does not exist") 
            unless( -e $opts->{'input_polypeptide_multi_fsa'} );
        $polypeptide_multi_fsa = $opts->{'input_polypeptide_multi_fsa'};
    }

    #Where phage_finder is installed
    unless( $opts->{'phage_finder_dir'} ) {
        $PHAGE_FINDER_DIR = "/usr/local/devel/ANNOTATION/phage_finder";
    } else {
        $logger->logdie("Input option phage_finder_dir is not a directory or does not exist") 
            unless( -d $opts->{'phage_finder_dir'} );
        $PHAGE_FINDER_DIR = $opts->{'phage_finder_dir'};
    }

}

sub get_files_from_list {
    my ($list_file, $array_ref) = @_;
    open( IN, "< $list_file") or
        $logger->logdie("Unable to open $list_file ($!)");
    chomp( @{$array_ref} = <IN> );
    close(IN);    
}

sub pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

