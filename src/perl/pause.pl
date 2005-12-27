#!/usr/local/bin/perl

=head1  NAME 

pause.pl - Spins until a condition is met

=head1 SYNOPSIS

USAGE:  pause.pl -f file_name -e email_nofication

=head1 OPTIONS

=cut 

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

use Pod::Usage;
BEGIN {
    require '/usr/local/devel/ANNOTATION/ard/chado-v1r5b1/lib/site_perl/5.8.5/Workflow/Logger.pm';
    import Workflow::Logger;
}

my %options = ();
my $results = GetOptions (\%options, 
			  'file|f=s',
			  'email|e',
			  'debug=s',
			  'log|l=s',
			  'help|h') || pod2usage();

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

if(! (-s $options{'file'})){
    &create_file($options{'file'});
    &email_notification($options{'file'},"Awaiting user intervention") if($options{'email'});
}

while(&check_file($options{'file'})){
    sleep 3000;
}
&email_notification($options{'file'},"User intervention complete") if($options{'email'});


sub create_file{
    my ($file) = @_;
    open FILE, "+>$file" or $logger->get_logger()->logdie("Can't create file $file: $?");
    print FILE "interrupted";
    close FILE;
}

sub check_file{
    my ($file) = @_;
    open FILE, "$file" or $logger->get_logger()->logdie("Can't read file $file: $?");
    while (my $line=<FILE>){
	if($line =~ /^complete/){
	    return 0;
	}
    }
    close FILE;
    return 1;
}

sub email_notification{
    my ($file,$msg) = @_;
    my $user = `whoami`;
    chomp $user;
    print STDERR "Sending notification email $msg to user $user\n";
    print `mail $user\@tigr.org -s "$msg ($file)" < /dev/null`;
}
 
