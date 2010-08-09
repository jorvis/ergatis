#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

mother.pl - Wrapper to facilitate the execution of mothur through ergatis.

=head1 SYNOPSIS

USAGE ./mothur.pl --mothur_exe=/path/to/mothur --input_file=/path/to/input --output_file=/path/to/output --args=<MOTHUR ARGS>

=head1 OPTIONS

B<--mothur_exe, -m>
	Path to the mothur executable.
	
B<--input_file, -i> 
	Input file that to be fed into mothur
	
B<--output_dir, -o>
	Desired directory where mothur output should go

B<--args, -a> 
	Any command-line arguments that should be passed to mothur

B<--log, -l>
	OPTIONAL. Log file
	
=head1 DESCRIPTION

This wrapper is written to get around the funky way that mothur handles where it writes any output to. The output directory
is not user-controlled in the normal sense with output being written to the directory where the primary input file is house.
The primary input file is different depending on the mothur command invoked.

=head1 INPUT

Input file for the mothur command to be invoked; this can be anything from a dist file to a FASTA file.

=head1 CONTACT

	Cesar Arze
	carze@som.umaryland.edu

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Ergatis::Logger;
use File::Copy;
use File::Basename;

my %options = ();
my $results = GetOptions (\%options,
						   'mothur_exe|m=s',
						   'input_file|i=s',
						   'output_dir|o=s',
						   'args|a=s',
						   'log|l=s',
						   'help|h') || pod2usage();
						   
if ( $options{'help'} ) {
	pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}						   
						   					
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger( 'LOG_FILE' 	=>	$logfile,
							   	  'LOG_LEVEL'	=>	$options{'debug'} );							   	  
$logger = Ergatis::Logger::get_logger();

## Because mothur handles where output is generated to a bit wonky
## we not only need to copy our input file to the desired output
## directory but also change the working directory to the desired output
## directory
my ($name) = basename($options{'input_file'});
chdir($options{'output_dir'});
my $temp_input = $options{'output_dir'} . "/" . $name;
copy_input_files($options{'input_file'}, $temp_input);

## Check args being passed into mothur
$options{'args'} =~ s/\,\s*\)$/\)/;
my $args = &check_args($options{'args'});
my $cmd = "$options{'mothur_exe'} \"$args\""; 
run_system_cmd($cmd);

unlink($temp_input);

## mothur operates in an interactive shell mode so 
## we need to scan the mothur log file looking for 
## an "Error - ..." line to see if our analysis 
## completed successfully.
&check_mothur_logfile($options{'output_dir'});

#####################################################
#													#
#				    SUBROUTINES						#
#													#
#####################################################

sub check_mothur_logfile {
    my $out_dir = shift;

    # Need to find the logfile
    opendir(DIR, $out_dir) or $logger->logdie("Could not open output folder $out_dir");
    my @files = grep { /\.(logFile|logfile)/ && -f "$out_dir/$_" } readdir(DIR);
    @files = map { $out_dir . "/" . $_ } @files;
    close (DIR);
    
    ## Horrible assumption to make here, but only one log file should exist per 
    ## output directory so we should be ok pulling the first element out of our 
    ## files array
    open (LOGFILE, $files[0]) or $logger->logdie("Could not open mothur logFile $files[0]: $!");
    while (my $line = <LOGFILE>) {
        chomp ($line);
        if ($line =~ /Error/) {
            my $err_msg = ( split(":", $line) )[1];
            $logger->logdie("mothur execution failed -- $line");
        }
    }

    close (LOGFILE);
}

sub check_args {
    my $raw_cmd = shift;
    my @tokens = split(";", $raw_cmd);
    my $formatted_cmd;

    foreach my $token (@tokens) {
        $token =~ s/(^\s*|\s*$)//;
        $token =~ /(.*)\((.*)\)/;
        my $mothur_cmd = $1;
        my $raw_args = $2;
    
        my @args = split(",", $raw_args);
        my $formatted_args;
        foreach my $arg (@args) {
            my ($key, $val) = split("=", $arg);
            $key =~ s/(^\s*|\s*$)//;
            $val =~ s/(^\s*|\s*$)//;

            ## If the key exists but the value is null exclude the arg
            next unless ($val);
            $formatted_args .= "$key=$val, ";    
        }
    
        $formatted_args =~ s/(\s*)?\,(\s*)$//;
        $formatted_args = "(" . $formatted_args . ")";
        $formatted_cmd .= $mothur_cmd . $formatted_args . "; ";
    }        
    
    $formatted_cmd =~ s/\;\s*$//;
    return $formatted_cmd;
}

sub copy_input_files {
	my ($file, $output_file) = @_;

	eval {
		copy($file, $output_file);
	};
	if ($@) {
		$logger->logdie("Could not copy input file to output directory: $@");
	}
}

sub run_system_cmd {
	my $cmd = shift;
	my $res = system($cmd);
	$res = $res >> 8;
	
    unless ($res == 0) {
        unlink($temp_input);
		$logger->logdie("Could not run $cmd");
	}
	
	return $res;
}
