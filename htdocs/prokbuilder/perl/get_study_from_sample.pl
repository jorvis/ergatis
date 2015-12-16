#!/usr/bin/perl

use strict;
use warnings;
use IPD::Client;
use IPD::IPDObject::Sample;

IPD::Client->set_client('production');
my $sample = defined($ARGV[0]) ? $ARGV[0] : die ("missing sample ID\n");

my @id;
my $study_id;
my $study;
my $category;

# make this either reference sample objects or project objects depending on if a project ID was specified in the initial args
$sample = IPD::IPDObject::Sample->get_sample($sample);
$study_id = defined( $sample->get_study_id() ) ? $sample->get_study_id() : 0;
print $study_id;

exit (0)
