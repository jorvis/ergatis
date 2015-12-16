#!/usr/bin/perl

use strict;
use warnings;
use IPD::Client;
use IPD::IPDObject::Project;
use IPD::IPDObject::Sample;
use IPD::IPDObject::Study;

my $client = defined($ARGV[0]) ? $ARGV[0] : die ("No client specified.  Please contact sadkins\@igs.umaryland.edu\n");
IPD::Client->set_client($client);
my $project = $ARGV[1] if (defined ($ARGV[1]));

my @id;
my $obj;
my $category;

# make this either reference sample objects or project objects depending on if a project ID was specified in the initial args
$obj = defined ($project) ? IPD::IPDObject::Sample->get_samples_by_project($project) : IPD::IPDObject::Project->get_all_projects();

if (defined $project) {
    foreach my $o(@{$obj}) {
	my $s = $o->get_study_id;
	my $study = IPD::IPDObject::Study->get_study($s);
	my $status = $study->get_status;
	if ($status eq 'active') {
            print $o->get_organism_name . " " . $o->get_strain . "\t(ID: " . $s . ")\n" if defined($o->get_organism_name);
	    print $o->get_name . "\t(ID: " . $s . ")\n" unless defined($o->get_organism_name);
	}
    }    
} else {
    foreach my $o(@{$obj}) {
        print $o->get_name . "\t(" . $o->get_lims_project_id . ")" . "\t---" . $o->get_id . "---\n" if defined($o->get_lims_project_id);
        print $o->get_name . "\t(12345)" . "\t---" . $o->get_id . "---\n" unless defined($o->get_lims_project_id);
    }     
}


exit (0)
