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
my $ss_id = -1;
my $ss;
my $category;

# make this either reference sample objects or project objects depending on if a project ID was specified in the initial args
$ss = IPD::IPDObject::StudyStage->get_study_stages_by_study($study);
foreach my $stage (@{$ss}) {
    if ($stage->get_type() eq 'automatic annotation') {
    	$ss_id = $stage->get_id();
	last;
    }
}

#if ($ss_id  == -1) {

#	$ss = IPD::IPDObject::StudyStage->new(
#		'type' => 'automatic annotation',
#		'study_id' => $study, 
#		'status' => 'active',
#	);
	#$ss = $client->create_study_stage($ss);
	#$ss_id = $ss->get_id();
#}

print $ss_id;

exit (0);

