#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{
	foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}
	$ENV{'GRID_CONFIG'} = "/var/www/.grid_request.conf";
};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

run_qsub_cmd.pl - Wrapper to facilitate the execution of a qsub command through ergatis.

=head1 SYNOPSIS

USAGE: ./run_qsub_cmd.pl --cmd=<QSUB CMD GOES HERE> --output_dir=/path/to/output/directory

=head1 OPTIONS

B<--cmd, -c>
	The qsub command to be executed.

B<--output_dir, -o>
    Output directory.
	
B<--log, -l>
	Log file
	
B<--debug, -d>
	Debug level
	
=head1 DESCRIPTION

This wrapper is written primarily to support the use of the celera-assembler in sge grid mode. The runCA executable
issues qsub commands to be issued which are fed into this wrapper script to be executed under ergatis.

=head1 INPUT

Either a string command or a file containing a command.

=head1 CONTACT

	Cesar Arze
	carze@som.umaryland.edu
	
=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Ergatis::Logger;
use Grid::Request;
use Tie::File;

## TODO: Remove this hard-coded path
my $HARVEST_EXE = "/opt/vappio-scripts/sge/submit_dir_for_harvesting.sh";

my %options = ();
my $results = GetOptions (\%options,
						   'cmd|c=s',
                           'output_dir|o=s',
						   'log|l=s',
						   'debug=s',
						   'help|h') || pod2usage();
						   
if ( $options{'help'} ) {
	pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
						  
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger( 'LOG_FILE' 	=>	$logfile,
							   	  'LOG_LEVEL'	=>	$options{'debug'} );							   	  
$logger = Ergatis::Logger::get_logger();

my $raw_cmd = $options{'cmd'};
my $cmd_string;

# Check if the cmd is a file and if so extract the qsub command.
if (-f $raw_cmd && -e $raw_cmd && -r $raw_cmd) {
	$cmd_string = &parse_cmd_from_file($raw_cmd);				
} else {
	&verify_qsub_cmd($raw_cmd);
	$cmd_string = $raw_cmd;
}

run_grid_job($cmd_string, $options{'output_dir'});

#####################################################
#						    #
#	            SUBROUTINES                     #
#						    #
#####################################################

sub run_grid_job {
	my ($cmd, $harvest_dir) = @_;
    my $work_dir = "/mnt/scratch/";
    my $log_dir = "/mnt/scratch/";

    # We need to extract all the information we need from the qsub command.
    $cmd =~ s/\\//g;
    $cmd =~ /.*\-N\s+(\w+)\s+\-t\s+(.*)\s+\-j.*\-o\s+(.*)\s+(.*)/;
    my $job_name = $1;
    my $array_num = $2;
    my $output_dir = $3;
    my $CA_cmd = $4;

    $output_dir =~ s/\s*$//;
    my $times = ( split("-", $array_num) )[1];
    chomp($times);
    $times =~ s/\s*$//;

    ## We need to append our custom harvesting step to the CA shell script
    ## TODO: Think of a non-hackish way to do this.
    modify_CA_shell_script($CA_cmd, $harvest_dir, $times);
	
    my $request = Grid::Request->new(
                                        opsys       =>  'Linux',
                                        initialdir  =>  "$work_dir",
                                        project     =>  "global",
                                        command     =>  "$CA_cmd",
                                        output      =>  "$output_dir",
                                        error       =>  "$log_dir/" . 'CA.err',
                                        times       =>  "$times"
                                    );
    

    ## Check if times is equal to 1, if so Grid::Request will not submit 
    ## the job as a an array job, so we must export  the environmental variable 
    ## $TASK_ID as one.     
    $request->set_env_list( "CLOVR_TASK_ID=1" ) if ($times eq "1"); 

    $request->pass_through("-S /bin/sh -q exec.q -cwd -N $job_name -j y -b n");
    $request->submit_and_wait();
 
    my $id_refs = $request->ids();
    my $success = &get_overall_status($id_refs, $request);
    
    unless ($success) {
        $logger->logdie("Error executing grid job");
    }
}

sub get_overall_status {
    my ($id_refs, $request) = @_;
    my $success;

    foreach my $id (@$id_refs) {
        if ($request->get_status($id) eq "FINISHED") {
            $success = 1;
        } else {
            $success = 0;
            last;
        }
    }

    return $success;
}

sub modify_CA_shell_script {
    my ($file, $dir, $times) = @_;
    my @bash_script = ();

    tie @bash_script, 'Tie::File', $file or $logger->logdie("Could not open file $file: $!");
    my $last_line = pop(@bash_script);
    
    ## if the last line is exit we need to add the harvesting script before exit
    ## otherwise the last line might be something important (i.e. post processing to output)
    ## and should occur before we call the harvesting script.
    if ($last_line =~ /exit/i) {
        push (@bash_script, "$HARVEST_EXE $dir");
	push (@bash_script, $last_line);
    } else {
       push (@bash_script, $last_line);
       push (@bash_script, "$HARVEST_EXE $dir");
    }

    ## If the number of exec nodes in use is equal to 1 we must manually pass the 
    ## SGE_TASK_ID env variable to whatever shell script is being executed. Because
    ## SGE does not allow the user to touch the SGE_TASK_ID env var we must use 
    ## A custom env varibale (CLOVR_TASK_ID) and modify the shell script to pull the jobid
    ## from there.
    if ($times eq "1") {
	foreach my $line (@bash_script) {
	    if ($line =~ /jobid=\$SGE_TASK_ID/) {
		$line = "jobid=\$CLOVR_TASK_ID";
	    }
	} 
    }

    untie @bash_script;
}

# TODO: Need to check if we get two qsub commands; should error out.
sub parse_cmd_from_file {
	my $file = shift;
	my (@cmd, $qsub_cmd);
	
	open (CMDFILE, $file) or $logger->logdie("Could not open file $file: $!");
	while (my $line = <CMDFILE>) {
		chomp($line);

        if ($line =~ /qsub/) {
            $qsub_cmd = $line;
            
            @cmd = <CMDFILE>;
            my $cmd_string = join("", @cmd);
            $cmd_string =~ s/\n/ /g;
            $qsub_cmd .= " " . $cmd_string;
		    $qsub_cmd =~ s/(.*.sh)\s(.*)/$1/;

            &verify_qsub_cmd($line);		
            last;
        } 
	}

    return $qsub_cmd;
}	

# TODO: Flesh this out a little more.
sub verify_qsub_cmd {
	my $cmd = shift;
	
	if ($cmd !~ /qsub/) {
		$logger->logdie("qsub command not well formed.");
	}
}

