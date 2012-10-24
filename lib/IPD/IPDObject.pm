package IPD::IPDObject;

=head1 NAME 

IPD::IPDObject - a parent class used to initialize the individual IPDObject modules

=cut

use strict;
use warnings;
use IPD::Client;
use XML::Simple;
use Data::Dumper;
use IPD::Net;

sub new {
    my ($class, %args) = @_;
    my $self = {}; ###
    bless($self, $class);
    $self->_init(\%args);	###
    return $self;
}

sub _init {
    my ($self, $elements) = @_;
    foreach (keys %{$elements}) {
	my $cmd = 'set_' . $_;
	$self->$cmd(${$elements}{$_});
    } 
}

=item construct_xml($hash_ref, $obj_type) - private method

B<Description:> Builds simple XML code from a passed in hash_ref. 

B<Parameters:> --First argument is a hash_ref of key_value pairs from the object (essentially the object itself can be passed in)
 --Second argument is a name associated with the type of object that will served as the outermost XML code (for creating the correct type XML for the requested IPD entry type)

B<Returns:> The XML code itself

=cut

sub construct_xml {
    my ($class, $elem, $top) = @_;
    my $xml = '<?xml version="1.0" encoding="UTF-8"?>';
    $xml .= "\n<$top>\n";
    foreach (keys %{$elem}) {
	next if (ref ${$elem}{$_} eq 'HASH');		### deals with nested objects... can deal with in later updates
	if (defined (${$elem}{$_})) {
	    $xml .= "\t<$_>${$elem}{$_}</$_>\n";
	}
    }
    $xml .= "</$top>";
    return $xml;
}

=item XML2Object($xml) - private method

B<Description:> Converts XML code into key-value pairs using the XML::Simple module's XMLIn subroutine.  The SuppressEmpty parameter is enabled to prevent unused XML keys from placed as keys in the returned hash-ref.  The KeyAttr => [] ensures contents isn't folded into arrays by "name", "key", or "id" values.  The NoAttr option is enabled to make sure no XML tag attrbutes get placed into the hash.  If one wants to view the $hash_ref, the user can use the Data::Dumper module and 'print Dumper($hash_ref)' to view the data structure.

B<Parameters:> The XML that needs to be converted

B<Returns:> A hash reference that essentially becomes the attributes for an object

=back 

=cut


sub XML2Object {	#Converts XML representation into a hash-ref
    my ($class, $xml) = @_;
    my $elements = XMLin($xml, SuppressEmpty => 1, KeyAttr => [], NoAttr => 1);	# Convert XML->hashref
    # if there are multiple entries for an XML query, KeyAttr => [] will fold these entries into an array of hashes.  Otherwise it will remain a standard hash.

### Data::Dumper used for printing hash for debugging purposes
open FILE, ">/home/sadkins/logs/dump_ipd_xml.txt" or die;
print FILE Dumper($elements);
close FILE;
    return $elements;	#can use Data::Dumper to read hash
}

1;
