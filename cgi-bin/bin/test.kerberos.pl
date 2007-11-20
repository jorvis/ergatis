#!/usr/bin/perl

use Authen::Simple::Kerberos;
use Term::ReadKey;
use Log::Cabin;
$|++;

my $realm = 'SOM.UMARYLAND.EDU';

print "\ntesting kerberos realm $realm\n\n";

print "enter user name: ";
my $username = <STDIN>;
chomp $username;

print "password: ";

ReadMode('noecho');
my $password = ReadLine(0);
chomp $password;
ReadMode('restore');

print "\n\n";

my $logsimple = new Log::Cabin;
   $logsimple->level( 8 );
   $logsimple->set_output(*STDOUT);
   
my $logger = $logsimple->get_logger('kerberos');

my $kerberos = Authen::Simple::Kerberos->new(
    realm => $realm,
    log => $logger,
);

if ( $kerberos->authenticate( $username, $password ) ) {
    print "\n\nin like flynn\n\n";

} else {
    print "\n\nno soup for you!\n\n";
}
