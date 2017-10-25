#!/usr/bin/env perl

use strict;
#use warnings;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my %opts = ();
GetOptions(\%opts,
            'tag-name|i=s',
            'pipeline_id|p=s',
            'url|u=s',
            'metadata|m:s',
            'log|l:s',
            'help') || pod2usage();

&pod2usage( {-exitval => 1, -verbose => 1, -output => \*STDOUT}) if ($opts{'help'} );

if($opts{url} =~ /clovr-/) {
    $opts{url} =~ s/clovr-//;
    $opts{url} =~ s/-/\./g;
}

my $metadata = $opts{metadata};

my @metadata_fields = split(/,/,$metadata);
push(@metadata_fields,"website=http://$opts{url}");
$metadata = join(" -m ",@metadata_fields);
my $cmd = "vp-add-dataset --overwrite --tag-name=".$opts{'pipeline_id'}."\_".$opts{'tag-name'}." -m $metadata";

print STDERR "$cmd\n";
print `$cmd`;
