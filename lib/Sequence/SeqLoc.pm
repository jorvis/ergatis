package Sequence::SeqLoc;

use strict;

use constant
{
	PLUS_STRAND	=> "+",
	MINUS_STRAND	=> "-"
};

sub new
{
	my ($class, $data) = @_;
	my $this = {};
	bless $this, $class;
	$this->Init($data);
	return $this;
}

sub Init
{
	my ($this, $data) = @_;
	if ($data) {
		if (UNIVERSAL::isa($data, 'Sequence::SeqLoc')) {
			$this->SetId($data->GetId);
			$this->SetFrom($data->GetFrom);
			$this->SetTo($data->GetTo);
			$this->SetStrand($data->GetStrand);
			$this->SetSeqData($data->GetSeqData);
		}
		else {
			my @tokens = split /\t/, $data;
			$this->SetId($tokens[0]);
			$this->SetFrom($tokens[1]);
			$this->SetTo($tokens[2]);
			$this->SetStrand($tokens[3]);
		}
	}
}

sub GetId
{
	my $this = shift;
	return $this->{ID};
}

sub SetId
{
	my ($this, $id) = @_;
	$this->{ID} = $id;
}

sub GetFrom
{
	my $this = shift;
	return $this->{FROM};
}

sub SetFrom
{
	my ($this, $from) = @_;
	$this->{FROM} = $from;
}

sub GetTo
{
	my $this = shift;
	return $this->{TO};
}

sub SetTo
{
	my ($this, $to) = @_;
	$this->{TO} = $to;
}

sub GetStrand
{
	my $this = shift;
	return $this->{STRAND};
}

sub SetStrand
{
	my ($this, $strand) = @_;
	$this->{STRAND} = $strand;
}

sub GetLength
{
	my $this = shift;
	return $this->{LENGTH} if $this->{LENGTH};
	return $this->{TO} - $this->{FROM};
}

sub SetLength
{
	my ($this, $length) = @_;
	$this->{LENGTH} = $length;
}

sub GetSeqData
{
	my $this = shift;
	return $this->{SEQ_DATA};
}

sub SetSeqData
{
	my ($this, $seqdata) = @_;
	$this->{SEQ_DATA} = $seqdata;
}

sub IsSetSeqData
{
	my $this = shift;
	return defined $this->{SEQ_DATA};
}

sub Overlaps
{
	my ($this, $loc) = @_;
	return 0 if not UNIVERSAL::isa($loc, 'Sequence::SeqLoc');
	if ($this->GetId ne $loc->GetId) {
		return 0;
	}
	if ($this->GetStrand ne $loc->GetStrand) {
		return 0;
	}
	my $loc1 = $this;
	my $loc2 = $loc;
	($loc1, $loc2) = ($loc2, $loc1) if $loc1->GetFrom > $loc2->GetFrom;
	return $loc1->GetTo >= $loc2->GetFrom;
}

1;
