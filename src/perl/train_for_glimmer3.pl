#!/usr/local/bin/perl

use lib(@INC, "/usr/local/devel/ANNOTATION/ard/current/lib/5.8.8");


=head1 NAME

train_for_glimmer.pl - you should put a one-line description of your script 
    here

=head1 SYNOPSIS

USAGE: train_for_glimmer.pl 
            --training_seqs=/path/to/somefile.fsa
            --icm_file=/path/to/glimmmer/training.icm
            --input_file=/path/to/input.fsa
            --output_file=/path/to/out.icm
          [ --log=/path/to/log.file
            --help
          ]

=head1 OPTIONS

B<--some_argument,-i>
    here you'll put a longer description of this argument that can span multiple lines. in 
    this example, the script has a long option and a short '-i' alternative usage.

B<--another_argument,-o>
    here's another named required argument.  please try to make your argument names all
    lower-case with underscores separating multi-word arguments.

B<--optional_argument,-s>
    optional.  you should preface any long description orf optional arguments with the
    optional word as I did in this description.  you shouldn't use 'optional' in the
    actual argument name like I did in this example, it was just to indicate that
    optional arguments in the SYNOPSIS should go between the [ ] symbols.

B<--optional_argument2,-f>
    optional.  if your optional argument has a default value, you should indicate it
    somewhere in this description.   ( default = foo )

B<--debug> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

put a longer overview of your script here.

=head1  INPUT

the input expectations of your script should be here.  pasting in examples of file format
expected is encouraged.

=head1  OUTPUT

the output format of your script should be here.  if your script manipulates a database,
you should document which tables and columns are affected.

=head1  CONTACT

    Kevin Galens
    (kgalens@tigr.org)

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Workflow::Logger;

########GLOBALS#############
use constant BUILD_ICM => '/usr/local/devel/ANNOTATION/glimmer/glimmer3.01/bin/build-icm';
my $training_seqs;
my $icm_file;
my $input_list;
my $output;
my $dontCreate = 0;
my $tmp_dir;
############################

my %options = ();
my $results = GetOptions (\%options, 
                          'training_seqs|t=s',
                          'input_list|i=s',
                          'output_file|o=s',
                          'tmp_dir|t=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

## make sure everything passed was peachy
&check_parameters(\%options);

unless($training_seqs || $icm_file) {
    
    #make the input file for build-icm
    open(IN, "< $input_list") or &_die("Cannot open input list");
    my @inputFiles = <IN>;
    close(IN);

    $training_seqs = &makeTrainingFile($tmp_dir, \@inputFiles);
                                      

}

unless($dontCreate) {

    &makeIcmFile($output, $training_seqs);

}

exit(0);

############################ SUB ROUTINES #######################################

sub makeTrainingFile {
    my ($tmp_dir, $inputFiles) = @_;
    my $outFile = "$tmp_dir/tmpFile.seqs";
    
    foreach my $file(@{$inputFiles}) {
        chomp $file;
        system("long-orfs $file > $tmp_dir/out.coords");
        system("extract $file $tmp_dir/out.coords >> $outFile")
    }

    system("rm -f $tmp_dir/out.coords");

    return $outFile;

}

sub makeIcmFile {
    my ($output) = @_;
    
    my $cmd = BUILD_ICM." $output < $training_seqs";

    system($cmd);

}

sub check_parameters {
    my $options = shift;
    
    if($options{'training_seqs'}) {
        &_die("option training_seqs does not exist") unless(-e $options{'training_seqs'});
        $training_seqs = $options{'training_seqs'};
    } elsif($options{'input_list'}) {
        &_die("option input_list does not exist") unless(-e $options{'input_list'});
        $input_list = $options{'input_list'};
    } else {
        &_die("one of input_list, training_seqs must exist");
    }

    unless($options{'output_file'}) {
        &_die("Options output_file is required");
    }
    $output = $options{'output_file'};

    unless($options{'tmp_dir'}) {
        &_die("Option tmp_dir is required");
    }
    $tmp_dir = $options{'tmp_dir'};
    $tmp_dir =~ s/\/$//;

    if(-e $output) {
        $dontCreate == 1;
    }

}

sub _die {
    my $msg = shift;
    $logger->logdie($msg);
}

########EOF
