package Blast::NcbiBlastHitDataType;

use strict;
use warnings;

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
		my @tokens = split /\t/, $data;
		$this->SetQueryId($tokens[0]);
		$this->SetSubjectId($tokens[1]);
		$this->SetPctId($tokens[2]);
		$this->SetAlignLen($tokens[3]);
		$this->SetNumMismatches($tokens[4]);
		$this->SetNumGaps($tokens[5]);
		$this->SetQueryStart($tokens[6]);
		$this->SetQueryEnd($tokens[7]);
		$this->SetSubjectStart($tokens[8]);
		$this->SetSubjectEnd($tokens[9]);
		$this->SetEValue($tokens[10]);
		$this->SetBitScore($tokens[11]);
	}
}

sub GetQueryId
{
	my $this = shift;
	return $this->{QUERY_ID};
}

sub SetQueryId
{
	my ($this, $query_id) = @_;
	$this->{QUERY_ID} = $query_id;
}

sub GetSubjectId
{
	my $this = shift;
	return $this->{SUBJECT_ID};
}

sub SetSubjectId
{
	my ($this, $subject_id) = @_;
	$this->{SUBJECT_ID} = $subject_id;
}

sub GetPctId
{
	my $this = shift;
	return $this->{PCT_ID};
}

sub SetPctId
{
	my ($this, $pct_id) = @_;
	$this->{PCT_ID} = $pct_id;
}

sub GetAlignLen
{
	my $this = shift;
	return $this->{ALIGN_LEN};
}

sub SetAlignLen
{
	my ($this, $align_len) = @_;
	$this->{ALIGN_LEN} = $align_len;
}

sub GetNumMismatches
{
	my $this = shift;
	return $this->{NUM_MISMATCHES};
}

sub SetNumMismatches
{
	my ($this, $num_mismatches) = @_;
	$this->{NUM_MISMATCHES} = $num_mismatches;
}

sub GetNumGaps
{
	my $this = shift;
	return $this->{NUM_GAPS};
}

sub SetNumGaps
{
	my ($this, $num_gaps) = @_;
	$this->{NUM_GAPS} = $num_gaps;
}

sub GetQueryStart
{
	my $this = shift;
	return $this->{QUERY_START};
}

sub SetQueryStart
{
	my ($this, $query_start) = @_;
	$this->{QUERY_START} = $query_start;
}

sub GetQueryEnd
{
	my $this = shift;
	return $this->{QUERY_END};
}

sub SetQueryEnd
{
	my ($this, $query_end) = @_;
	$this->{QUERY_END} = $query_end;
}

sub GetSubjectStart
{
	my $this = shift;
	return $this->{SUBJECT_START};
}

sub SetSubjectStart
{
	my ($this, $subject_start) = @_;
	$this->{SUBJECT_START} = $subject_start;
}

sub GetSubjectEnd
{
	my $this = shift;
	return $this->{SUBJECT_END};
}

sub SetSubjectEnd
{
	my ($this, $subject_end) = @_;
	$this->{SUBJECT_END} = $subject_end;
}

sub GetEValue
{
	my $this = shift;
	return $this->{EVALUE};
}

sub SetEValue
{
	my ($this, $evalue) = @_;
	$this->{EVALUE} = $evalue;
}

sub GetBitScore
{
	my $this = shift;
	return $this->{BITSCORE};
}

sub SetBitScore
{
	my ($this, $bitscore) = @_;
	$this->{BITSCORE} = $bitscore;
}

1;
