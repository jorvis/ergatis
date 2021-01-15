#!/usr/local/bin/perl

use strict;
use FindBin qw( $RealBin );
use lib $RealBin;
use warnings;
use IPD::Client;
use IPD::IPDObject::StudyStage;

my $client = defined($ARGV[0]) ? $ARGV[0] : die ("No client specified.  Please contact sadkins\@igs.umaryland.edu\n");
IPD::Client->set_client($client);
my $study = defined($ARGV[1]) ? $ARGV[1] : die ("missing study ID\n");

my @id;
my $ss;
my $category;

# make this either reference sample objects or project objects depending on if a project ID was specified in the initial args
$ss = IPD::IPDObject::StudyStage->new (
		'study_id' => $study,
		'type' => 'automatic annotation',
		'status' => 'active'
	);

$ss = $ss->create_study_stage();

print $ss->get_id();

exit (0);

