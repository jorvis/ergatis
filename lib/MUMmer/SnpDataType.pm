package MUMmer::SnpDataType;

use strict;
use warnings;

my @header	= qw(P1 SUB1 SUB2 P2 BUFF DIST R Q LEN_R LEN_Q CTX_R CTX_Q
		     FRM1 FRM2 TAG1 TAG2);

# class methods
sub GetHeader;
sub SetHeader;
sub GetNumberOfColumns;
sub new;

# instance methods
sub Init;
sub PrintAll;
sub GetQueryId;
sub SetQueryId;
sub GetSubjectId;
sub SetSubjectId;
sub GetQueryPosition;
sub SetQueryPosition;
sub GetSubjectPosition;
sub SetSubjectPosition;
sub GetQuerySubstitution;
sub SetQuerySubstitution;
sub GetSubjectSubstitution;
sub SetSubjectSubstitution;
sub GetDistanceToNearestMismatch;
sub SetDistanceToNearestMismatch;
sub GetDistanceToNearestEnd;
sub SetDistanceToNearestEnd;
sub GetNumberOfQueryRepeatAlignments;
sub SetNumberOfQueryRepeatAlignments;
sub GetNumberOfSubjectRepeatAlignments;
sub SetNumberOfSubjectRepeatAlignments;
sub GetQueryLength;
sub SetQueryLength;
sub GetSubjectLength;
sub SetSubjectLength;
sub GetQueryContext;
sub SetQueryContext;
sub GetSubjectContext;
sub SetSubjectContext;
sub GetQueryFrame;
sub SetQueryFrame;
sub GetSubjectFrame;
sub SetSubjectFrame;

sub GetHeader
{
	return join "\t", @header;
}

sub SetHeader
{
	@header = ();
	my $data = shift;
	chomp $data;
	my @tokens = split /\t/, $data;
	for (my $i = 0; $i < scalar(@tokens); ++$i) {
		my $col_type = $tokens[$i];
		$col_type =~ s/ +/_/g;
		$col_type =~ s/[\[\]]+//g;
		if ($col_type eq 'SUB') {
			push @header, 'SUB1', 'SUB2';
			++$i;
		}
		elsif ($col_type eq 'FRM') {
			push @header, 'FRM1', 'FRM2';
		}
		elsif ($col_type eq 'TAGS') {
			push @header, 'TAG1', 'TAG2';
		}
		else {
			push @header, $col_type;
		}
	}
}

sub GetNumberOfColumns
{
	return scalar(@header);
}

sub new
{
	my $class = shift;
	my $this = {};
	bless $this, $class;
	my $data = shift;
	$this->Init($data) if defined $data;
	return $this;
}

sub Init
{
	my ($this, $data) = @_;
	chomp $data;
	my @tokens = split /\t/, $data;
	for (my $i = 0; $i < scalar(@tokens); ++$i) {
		$this->{$header[$i]} = $tokens[$i];
	}
	--$this->{P1};
	--$this->{P2};
}

sub PrintAll
{
	my $this = shift;
	while (my ($key, $val) = each %{$this}) {
		print "$key\t= $val\n";
	}
}

sub GetQueryId
{
	my $this = shift;
	return $this->{TAG1};
}

sub SetQueryId
{
	my $this = shift;
	$this->{TAG1} = shift;
}

sub GetSubjectId
{
	my $this = shift;
	return $this->{TAG2};
}

sub SetSubjectId
{
	my $this = shift;
	$this->{TAG2} = shift;
}

sub GetQueryPosition
{
	my $this = shift;
	return $this->{P1};
}

sub SetQueryPosition
{
        my $this = shift;
        $this->{P1} = shift;
}

sub GetSubjectPosition
{
	my $this = shift;
	return $this->{P2};
}

sub SetSubjectPosition
{
        my $this = shift;
        $this->{P2} = shift;
}

sub GetQuerySubstitution
{
	my $this = shift;
	return $this->{SUB1};
}

sub SetQuerySubstitution
{
        my $this = shift;
        $this->{SUB1} = shift;
}

sub GetSubjectSubstitution
{
	my $this = shift;
	return $this->{SUB2};
}

sub SetSubjectSubstitution
{
        my $this = shift;
        $this->{SUB2} = shift;
}

sub GetDistanceToNearestMismatch
{
	my $this = shift;
	return $this->{BUFF};
}

sub SetDistanceToNearestMismatch
{
        my $this = shift;
        $this->{BUFF} = shift;
}

sub GetDistanceToNearestEnd
{
	my $this = shift;
	return $this->{DIST};
}

sub SetDistanceToNearestEnd
{
        my $this = shift;
        $this->{DIST} = shift;
}

sub GetNumberOfQueryRepeatAlignments
{
	my $this = shift;
	return $this->{R};
}

sub SetNumberOfQueryRepeatAlignments
{
        my $this = shift;
        $this->{R} = shift;
}

sub GetNumberOfSubjectRepeatAlignments
{
	my $this = shift;
	return $this->{Q};
}

sub SetNumberOfSubjectRepeatAlignments
{
        my $this = shift;
        $this->{Q} = shift;
}

sub GetQueryLength
{
	my $this = shift;
	return $this->{LEN_R};
}

sub SetQueryLength
{
        my $this = shift;
        $this->{LEN_R} = shift;
}

sub GetSubjectLength
{
	my $this = shift;
	return $this->{LEN_Q};
}

sub SetSubjectLength
{
        my $this = shift;
        $this->{LEN_Q} = shift;
}

sub GetQueryContext
{
	my $this = shift;
	return $this->{CTX_R};
}

sub SetQueryContext
{
        my $this = shift;
        $this->{CTX_R} = shift;
}

sub GetSubjectContext
{
	my $this = shift;
	return $this->{CTX_Q};
}

sub SetSubjectContext
{
        my $this = shift;
        $this->{CTX_Q} = shift;
}

sub GetQueryFrame
{
	my $this = shift;
	return $this->{FRM1};
}

sub SetQueryFrame
{
        my $this = shift;
        $this->{FRM1} = shift;
}

sub GetSubjectFrame
{
	my $this = shift;
	return $this->{FRM2};
}

sub SetSubjectFrame
{
        my $this = shift;
        $this->{FRM2} = shift;
}

1;
