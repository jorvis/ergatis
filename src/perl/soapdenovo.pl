#!/usr/bin/env perl

eval 'exec /usr/bin/env perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/env perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

soapdenovo.pl - Wrapper to facilitate the execution of SOAPdenovo through ergatis.

=head1 SYNOPSIS

USAGE ./soapdenovo.pl --soapdenovo_exe=/path/to/soapdenovo --process=soapdenovo_process --config_file=/path/to/config --args=<SOAPDENOVO ARGS>

=head1 OPTIONS

B<--soapdenovo_exe, -s>
        Path to the soapdenovo executable.

B<--process, -p> 
        Choose soapdenovo process to run.
	Options include all, pregraph, contig, map, scaff
        
B<--config_file, -c>
        Config file with info location of input files, library info, etc

B<--args, -a> 
        Any command-line arguments that should be passed to soapdenovo

B<--log, -l>
        OPTIONAL. Log file
        
=head1 DESCRIPTION

The script executes the SOAPdenovo script from the SOAPdenovo software package (version 1.04).

=head1 INPUT

The main input for this script is the config file which has the following format:
	#maximal read length
	max_rd_len=75
	[LIB]
	avg_ins = <LIBSZ> # average insert size
	reverse_seq = 0 # forward/reverse library
	asm_flags = 3 # reads used for contigging and scaffolding
	pair_num_cutoff = 2 # num of mates needed to scaffold across a gap
	map_len = 60 # minimum length of read mapping to a contig
	q1 = <prefix>.1.fastq # read1 in fastq
	q2 = <prefix>.2.fastq # read2 in fastq
	q=<prefix>.singleton.fastq

more details: http://soap.genomics.org.cn/soapdenovo.html

=head1 CONTACT

        Kemi Abolude
        kabolude@som.umaryland.edu

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Ergatis::Logger;
use File::Copy;
use File::Basename;

my %options = ();
my $results = GetOptions (\%options,
			'soapdenovo_exec|s=s',
			'process|p=s',
                        'config_file|c=s',
                        'args|a=s',
                        'log|l=s',
                        'help|h') || pod2usage();
                                                   
if ( $options{'help'} ) {
        pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}                                                  
                                                                                        
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger( 'LOG_FILE'    =>      $logfile,
                                  'LOG_LEVEL'   =>      $options{'debug'} );                    
                                          
$logger = Ergatis::Logger::get_logger();
               

my $process;
my $soapdenovo_exec;
my $config_file;
my $args;


## make sure everything passed was peachy
check_parameters(\%options);


## Check args being passed into mothur
#my $cmd = "$options{'soapdenovo_exe'} $options{'process'} -s $options{'config_file'} \"$options{'args'}\""; 
my $cmd = "$soapdenovo_exec $process -s $config_file $args";

##run command
run_system_cmd($cmd);



##Subroutines

sub run_system_cmd {
        my $cmd = shift;
        my $res = system($cmd);
        $res = $res >> 8;
        
    unless ($res == 0) {
#        unlink($temp_input);
                $logger->logdie("Could not run $cmd");
        }
        
        return $res;
}

sub check_parameters {
   my $opts = shift;
   my $logfh;
   my $ofh;

   if( $opts->{'help'} ) {
       &_pod;
   }

   if( $opts->{'log'} ) {
       open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
   }

   foreach my $req ( qw(soapdenovo_exec process config_file args) ) {
       die("ERROR: Option $req is required") unless( $opts->{$req} );
   }
   $soapdenovo_exec = $opts->{'soapdenovo_exec'};
   $process = $opts->{'process'};
   $config_file = $opts->{'config_file'};
   $args = $opts->{'args'} ;

}


                    
