#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

id_test.pl - generate IDs for testing

=head1 SYNOPSIS

USAGE: id_test.pl 
        --feat_type=exon
        --id_count=10000
        --id_repository=/path/to/someids
        --project=aa1 

=head1 OPTIONS

B<--feat_type> 
    SO type of features to pull.

B<--id_count> 
    Number of IDs to pull.

B<--id_repository>
    Required for creating feature identifiers.  Each project should have
    its own id_repository directory - use the full path to it here.  This
    is used by the IdGenerator.pm module.

B<--project> 
    Project name to use in the IDs pulled.

=head1   DESCRIPTION


=head1 INPUT


=head1 OUTPUT


=head1 CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Workflow::IdGenerator;
use Workflow::Logger;

my %options = ();
my $results = GetOptions (\%options, 
			  'feat_type=s',
              'id_count=i',
              'id_repository=s',
              'project=s',
			  'help|h') || pod2usage();

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

## make sure all passed options are peachy
&check_parameters(\%options);

## we're going to generate ids
my $idcreator = new Workflow::IdGenerator( id_repository => $options{id_repository} );

for (my $i=0; $i<$options{id_count}; $i++) {
    $idcreator->next_id( type => $options{feat_type}, project => $options{project} );
}


exit;


sub check_parameters {
    
    ## required params
    my @required = qw( feat_type id_count id_repository project );
    for ( @required ) {
        if (! defined $options{$_}) {
            $logger->logdie( "$_ is a required option" );
        }
    }
    
    return 1;
}

