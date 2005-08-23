package Sequence::Sequence;

## class methods
sub new;
sub ReverseComplement;

## instance methods
sub Init;
sub GetId;
sub SetId;
sub GetTitle;
sub SetTitle;
sub GetSequence;
sub SetSequence;
sub GetSubSequence;
sub Print;

sub new
{
	my ($class, $id, $seq, $title) = @_;
	my $this = {};
	bless $this, $class;
	$this->Init($id, $seq, $title);
	return $this;
}

sub Init
{
	my ($this, $id, $seq, $title) = @_;
	$this->{ID} = $id if $id;
	$this->{SEQ} = $seq if $seq;
	$this->{TITLE} = $title if $title;
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

sub GetTitle
{
	my $this = shift;
	return $this->{TITLE};
}

sub SetTitle
{
	my ($this, $title) = @_;
	$this->{TITLE} = $title;
}

sub GetSequence
{
	my ($this, $minus) = @_;
	return $minus ? ReverseComplement($this->{SEQ}) : $this->{SEQ};
}

sub SetSequence
{
	my ($this, $seq) = @_;
	$this->{SEQ} = $seq;
}

sub GetSubSequence
{
	my ($this, $from, $to, $minus) = @_;
	return undef if !defined($from) || !defined($to);
	my $sub_seq = substr($this->{SEQ}, $from, ($to - $from + 1));
	return $minus ? &ReverseComplement($sub_seq) : $sub_seq;
}

sub ReverseComplement
{
	my $seq = reverse(shift);
	$seq =~ tr/acgtACGTrykmRYKM/tgcaTGCAyrmkYRMK/;
	return $seq;
}

sub Print
{
	my ($this, $out) = @_;
	if (!$out) {
		$out = STDOUT;
	}
	print $out ">", $this->GetId;
	print $out "\t", $this->GetTitle if $this->{TITLE};
	print $out "\n";
	print $out $this->GetSequence, "\n";
}

1;
