package IPD::Client;

=head1 NAME 

IPD::Client - Module responsible for storing the URL path of the client.

=head1 DESCRIPTION

 IPD::Client->set_client('production') - sets client to the 'production' designated URl.
 IPD::Client->get_client() - retrieves the URL that has been set

 There are currently 2 valid client keywords:
 'production' => 'http://projects.igs.umaryland.edu/', 
 'devel' => 'http://projects-devel.igs.umaryland.edu/',

=head1 CONTACT

 Shaun Adkins
 sadkins@som.umaryland.edu

=over 4
 
=cut

my $client = "",

my %valid_ipd_clients = 
    ('production' => 'http://projects.igs.umaryland.edu/', 
     'devel' => 'http://projects-devel.igs.umaryland.edu/',
     'test' => 'http://projects-test.igs.umaryland.edu/'	#does not currently exist
     );

### Sets the client path
sub set_client {
    my $class = shift;
    $client = shift;
    if (exists $valid_ipd_clients{ $client }) {
	$client = $valid_ipd_clients{ $client };
    } else {
	die " \'$client\' is not a valid client name.  Please choose between \'production\', \'devel\', or \'test\': $!\n"; 
    }
}

### Retrieves the client path set in the class
sub get_client {
    my $class = shift;
    return $client;
}

1;
