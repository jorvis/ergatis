package IPD::Net;
use LWP::UserAgent;
use HTTP::Request;
use HTTP::Request::Common qw(PUT);

=head1 NAME 

IPD::Net - Module that takes parameters from IPD::Client and sends an XML request to the requested server through HTTP.  It then relays the content back to IPD::Client or it dies, printing the HTTP status code and the XML-embedded error message where applicable.

=cut

sub send_request {
    my $class = shift;
    my ($loc, $ext, $method, $xml) = @_;
    $url = $loc . $ext;
    $url =~ s/\s/+/g;	# substituting spaces with percent-encoded equivalent
    my $ua = LWP::UserAgent->new;
    $ua->timeout(999);	# some commands (get all projects and get all studies) may take a while
    $ua->env_proxy;
#print $url, "\n";
#print $xml, "\n\n";
    my $req = HTTP::Request->new($method => $url);
    $req->content_type('text/xml') if (defined $xml);
    $req->content($xml) if (defined $xml);
    my $res = $ua->request($req);

    if (! $res->is_success) {
	print $method, "\t", $res->status_line . "\n";    
	die $res->content, "\n";
    }
    return $res->content;




    #handle_300($1) if $res->status_line =~ /(3\d\d)/;
    #handle_400($1) if $res->status_line =~ /(4\d\d)/;
    #handle_500($1, $2) if $res->status_line =~ /(5\d\d)\s(.+^)/;

}



#### functions I have not implemented but are here in case we want to handle specific error codes
sub handle_300 {
    my $error = shift;
}

sub handle_400 {
    my $error = shift;

    print "Check the syntax of your URL\n" if $error == 400;
    print "Server doesn't understand header information\n" if $error == 406;
    print "Can't process entry\n" if $error == 422;
}

sub handle_500 {
    my $error = shift;
    my $msg = shift;
    chomp $msg;

    print "Unknown error occurred\n" if ($error == 500 && $msg =~ /Internal Server Error/);
    print "Get request took too long\n" if ($error == 500 && $msg =~ /read timeout/);
}

1;
