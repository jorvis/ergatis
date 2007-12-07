#!/usr/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

new_id.pl - generate a single id for a given project and feature type, using Ergatis::IdGenerator.

=head1 SYNOPSIS

USAGE: new_id.pl 
        --repository_root=/usr/local/annotation/EHA3
        --feat_type=exon
        --project=eha3

=head1 OPTIONS

B<--repository_root>
    The directory in which the files of --project are stored.

B<--feat_type>
    SO type of feature for which to generate an id/uniquename.

B<--project>
    Project/database name for new id.

B<--help>
    Display usage information.

B<--man>
    Display detailed help information.

=head1 DESCRIPTION

Generate a single id for a given project and feature type, using Ergatis::IdGenerator.

=cut

use strict;
use File::Spec;
use Getopt::Long;
use Pod::Usage;

use Ergatis::ConfigFile;
use Ergatis::IdGenerator;

# ------------------------------------------------------------
# Globals
# ------------------------------------------------------------
my $ID_REPOS_PROP = '$;PROJECT_ID_REPOSITORY$;';

# ------------------------------------------------------------
# Input
# ------------------------------------------------------------
my($repository_root, $feat_type, $project, $help, $man);

&GetOptions("repository_root=s" => \$repository_root,
	    "feat_type=s" => \$feat_type,
	    "project=s" => \$project,
 	    "help" => \$help,
	    "man" => \$man,
	    );

&pod2usage(1) if $help;
&pod2usage({-verbose => 2}) if $man;

# check that all required arguments are present
&pod2usage({-message => "Error: --feat_type not specified\n", -exitstatus => 1, -verbose => 0}) if (!defined($feat_type));
&pod2usage({-message => "Error: --repository_root not specified\n", -exitstatus => 1, -verbose => 0}) if (!defined($repository_root));
&pod2usage({-message => "Error: --project not specified\n", -exitstatus => 1, -verbose => 0}) if (!defined($project));

# ------------------------------------------------------------
# Main program
# ------------------------------------------------------------

# figure out location of id repository for this project
my $id_repository = undef;

# ergatis v1 uses workflow_config_files/sharedconf.ini, property name = [init]/$;PROJECT_ID_REPOSITORY$;
my $ergatisV1Config = File::Spec->catfile($repository_root, 'workflow_config_files', 'sharedconf.ini');

# ergatis v2 uses workflow/project.config, property name = [project]/$PROJECT_ID_REPOSITORY$;
my $ergatisV2Config = File::Spec->catfile($repository_root, 'workflow', 'project.config');

# NOTE - this code relies on ergatis v2 being able to parse v1 config files
# (since *this* script could be run from v2 on a v1 installation)

# check for v2 first in case a project has config. files for both;
if (-e $ergatisV2Config) {
    my $shared_cfg = new Ergatis::ConfigFile( -file => $ergatisV2Config );
    $id_repository = $shared_cfg->val('project', $ID_REPOS_PROP);
    die "could not find [project] property $ID_REPOS_PROP in $shared_cfg (ergatis v2)" if (!defined($id_repository));
} elsif (-e $ergatisV1Config) {
    my $shared_cfg = new Ergatis::ConfigFile( -file => $ergatisV1Config );
    $id_repository = $shared_cfg->val('init', $ID_REPOS_PROP);
    die "could not find [init] property $ID_REPOS_PROP in $shared_cfg (ergatis v1)" if (!defined($id_repository));
} else {
    die "couldn't find project config file at either $ergatisV1Config (ergatis v1) or $ergatisV2Config (ergatis v2)";
}

die "id_repository $id_repository does not appear to exist" if (!-e $id_repository);

# generate and print the new id
my $idgen = Ergatis::IdGenerator->new( id_repository => $id_repository, logging => 0 );
$idgen->set_pool_size( $feat_type => 1 );
my $id = $idgen->next_id( 'type' => $feat_type, 'project' => lc($project) );
print "$id\n";

exit(0);

