#!/usr/bin/perl

=head1 NAME

run_humann.pl - will create directory structure for the hummann pipeline and run it.

=head1 SYNOPSIS

USEAGE: run_humann.pl
       --meta_file=/path/to/metadata.dat
       --input_file=/path/to/input/file.txt
       --input_list=/path/to/input/file.list
       --input_format=txt.gz
       --max_id=1.0
       --hits=20
       --humann_dir=/path/to/humann/dir
       --output_dir=/path/to/output/dir
     [ --log
       --help
     ]

=head1 OPTIONS

B<--meta_file, -m>    
    Optional. Filename from which metadata annotations are read.
    
B<--input_file, -i>
    Input file. Format should be tab-delimited output file, from blastx (outfmt 6), mapx or mblastx. input_file or input_list_file is required.

B<--input_list, -i>
    Input list file. Format should be tab-delimited output file, from blastx (outfmt 6), mapx or mblastx. input_file or input_list_file is required.  

B<--input_format, -f>
    Required. Examples: txt, txt.gz, mapx.bz2, mapx.gz, mblastx.gz.

B<--max_id, -x>
    Optional. Keeps nothing above this identity. Default is 1.0

B<--hits, -t>
    Optional. Keep top t hits.

B<--humann_dir, -h>
    Required. Directory where human is installed.

B<--output_dir, -o>
    Required. Directory for HUMAnN output.

B<--log,-l>
    Log file

B<--help, h>
    This help message

=head1 DESCRIPTION

Creates directory structure and run the HUMAnN pipeline.

=head1 CONTACT

    Kemi Abolude
    kabolude@som.umaryland.edu

=cut


use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

my $meta_file;
my $input_list;
my $format;
my $id;
my $hits;
my $humann_dir;
my $output_dir;
my $file;
my @dirs = ();
my @files =();

my %options = ();
my $results = GetOptions (\%options, 
			  'meta_file|m:s',
			  'input_file|i=s',
			  'input_list|l=s',
                          'input_format|f=s',
                          'max_id|x:f',
                          'hits|t:i',
			  'humann_dir|d=s',
			  'output_dir|o=s',
                          'log|l=s',
                          'help|h') || pod2usage();

&check_options(\%options);

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "can't create log file: $!";
}


    #load files
    system ("cp -r $humann_dir/* $output_dir/");
    
    foreach my $file(@files) {
	system ("cp -p $file $output_dir/input/") ;
    }	
	
    my $meta_file_name;
    if (defined($meta_file)){
	my @fname = split(/\//, $meta_file);
	$meta_file_name =$fname[-1]; 
	system ("cp -p $meta_file $output_dir/input/") ;	
    }
    my $temp;
    #Create SConstruct file 
    open (SCONSTRUCT,">$output_dir/SConstruct") or die $!;
    open (FILE,"<$output_dir/temp.txt") or die $!;
    while (<FILE>){
        if ((/hmp_metadata.dat/) && defined($meta_file)){
		$_ =~ s/hmp_metadata.dat/$meta_file_name/; 
                print SCONSTRUCT "$_";
	}elsif(/txt.gz/){
		$_ =~ s/txt.gz/$format/;
		print SCONSTRUCT "$_";
        }elsif(/placeholder/){
		my @val = split("\\.", $format);
		if ($val[0] eq "txt"){
			$_ =~ s/placeholder//;
		}else{
		     if (defined($id)){
			if (defined($hits)){
				$_ =~ s/placeholder/\"$val[0]\", \"$id\", \"$hits\"/;
			}else{
				$_ =~ s/placeholder/\"$val[0]\", \"$id\"/;
			}
		      }else{	
                	$_ =~ s/placeholder/\"$val[0]\"/;
		      }
		    }
		print SCONSTRUCT "$_";
	}else {
	   print SCONSTRUCT "$_";
	}
    }
    close FILE;
    close SCONSTRUCT;

    #set python 2.7

    #run scons
    chdir ("$output_dir") or die "$!";
    system ("scons"); 


exit(0);


sub _log {
    my $msg = shift;
    print $logfh "$msg\n" if $logfh;
}

sub check_options {
    my ($opts) = @_;

    if( $opts->{'help'} ) {
        pod2usage();
        exit(0);
    }

    if($opts->{'input_list'} && $opts->{'input_list'} ne "") {
        &_die("input_list [$opts->{'input_list'}] does not exist") unless( -e $opts->{'input_list'});
        open(IN, "< $opts->{'input_list'}") or &_die("Unable to open $opts->{'input_list'}");
        while(<IN>) {
            push(@files, $_);
        }
        close(IN);
    } elsif($opts->{'input_file'} && $opts->{'input_file'} ne "") {
        &_die("input_file [$opts->{'input_file'}] does not exist") unless( -e $opts->{'input_file'});
        push(@files,$opts->{'input_file'});
    } else {
        &_die("Either input_list or input_file must be provided");
    }

    my @reqs = qw( input_format humann_dir output_dir);
    foreach my $req ( @reqs ) {
        die("Option $req is required") unless( exists( $opts->{$req} ) );
    }
    my @format_a = qw(txt txt.gz mapx.bz2 mapx.gz mblastx mblastx.gz);
    my %format_h = map { $_ => 1 } @format_a;
    unless(exists($format_h{$opts->{'input_format'}})) {
       die("Input format is not valid!"); 
    }

    $meta_file = $opts-> {'meta_file'};
    $format = $opts->{'input_format'};
    $id = $opts->{'max_id'};
    $hits = $opts->{'hits'};
    $humann_dir = $opts->{'humann_dir'};
    $output_dir = $opts->{'output_dir'};
}
