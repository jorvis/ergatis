#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

order_paired_list.pl - ensures list file orders files correctly: member 1 then member 2

=head1 SYNOPSIS

USAGE: order_paired_list.pl 
            --input_list=/path/to/somefile.list 
            --output_list=/path/to/somefile.list 

=head1 OPTIONS

B<--input_list,-i>
    Required. The input list file containing the paths of two paired files.

B<--output_list,-o>
    Optional. List file containing the paths of two paired files in the correct order: member 1 then 
    member 2. Output will be named ordered_files.list if user does not provide one.

=head1  DESCRIPTION

This script is used to reorder a list file containing two sequences. The script ensures that the list
file contains exactly two paths, and that they are in the correct order. 

=head1  INPUT

The input is defined with --input_list and should be a single list file. The file should contain paths 
to exactly two fasta or fasq files, each representing a member of the paired end reads. The header for 
on of the fasta/fastq files must end with "/1", while the second "/2".
 

=head1  CONTACT

    Kemi Abolude
    kabolude@som.umaryland.edu

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use POSIX;
BEGIN {
    use Ergatis::Logger;
}

my %options = ();
my $results = GetOptions (\%options, 
                          'input_list|i=s',
                          'output_list|o=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || pod2usage();

my $input_list;
my $output_list;
my $file1;
my $file2;

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				 'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

open(my $fh, "<$input_list") || $logger->logdie("couldn't open $input_list");

chomp ($file1 = <$fh>);
chomp ($file2 = <$fh>);

close($fh);

open(my $fh1, "$file1") || $logger->logdie("couldn't open $file1 for reading");
open(my $fh2, "$file2") || $logger->logdie("couldn't open $file2 for reading");

my $header1;
my $header2;
chomp ($header1 = <$fh1>);
chomp ($header2 = <$fh2>);

print "$header1\n$header2\n";

open(my $outfh, ">$output_list") || $logger->logdie("couldn't create $output_list list file");

my $sub1 = substr $header1, -2;
my $sub2 = substr $header2, -2;

if (($sub1 eq "/1") && ($sub2 eq "/2")) {
    print $outfh "$file1\n$file2";
}elsif (($sub1 eq "/2") && ($sub2 eq "/1")) {
    print $outfh "$file2\n$file1";
}else{
    $logger->logdie("The format for the fasta or fastq files are incorrect");
}

close($fh1);
close($fh2);
close($outfh);

exit(0);


sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( input_list);
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }

    ## make sure the input files exists
    if (! -e $options{'input_list'}) {
       $logger->logdie("input_list $options{'input_list'} does not exist")
       }
    $input_list = $options{'input_list'};

    ##
    if ($options{'output_list'}) {
	$output_list = $options{output_list};
    }else{
	$output_list = "ordered_files.list";
    }

}





