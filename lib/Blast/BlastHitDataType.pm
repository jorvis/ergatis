package Blast::BlastHitDataType;

use strict;
use warnings;

# class methods
sub new;

# instance methods
sub Init;
sub GetQueryName;
sub SetQueryName;
sub GetSubjectName;
sub SetSubjectName;
sub GetQueryBegin;
sub SetQueryBegin;
sub GetQueryEnd;
sub SetQueryEnd;
sub GetQueryLength;
sub SetQueryLength;
sub GetSubjectBegin;
sub SetSubjectBegin;
sub GetSubjectEnd;
sub SetSubjectEnd;
sub GetSubjectLength;
sub SetSubjectLength;
sub GetPercentId;
sub SetPercentId;
sub GetPercentSimilarity;
sub SetPercentSimilarity;
sub GetRawScore;
sub SetRawScore;
sub GetBitScore;
sub SetBitScore;
sub GetDescription;
sub SetDescription;
sub GetFrame;
sub SetFrame;
sub GetStrand;
sub SetStrand;
sub GetEValue;
sub SetEValue;
sub GetPValue;
sub SetPValue;
sub GetDatabaseName;
sub SetDatabaseName;
sub GetAlgorithm;
sub SetAlgorithm;
sub GetDate;
sub SetDate;
sub ToString;

sub new
{
	my ($class, $data) = @_;
	my $this = {};
	bless $this, $class;
	$this->Init($data) if $data;
	return $this;
}

sub Init
{
	my ($this, $data) = @_;
	if (UNIVERSAL::isa($data, 'Blast::BlastHitDataType')) {
		$data = $data->ToString();
	}
	my @tokens = split /\t/, $data;
	$this->SetQueryName($tokens[0]);
	$this->SetDate($tokens[1]);
	$this->SetQueryLength($tokens[2]);
	$this->SetAlgorithm($tokens[3]);
	$this->SetDatabaseName($tokens[4]);
	$this->SetSubjectName($tokens[5]);
	$this->SetQueryBegin($tokens[6]);
	$this->SetQueryEnd($tokens[7]);
	$this->SetSubjectBegin($tokens[8]);
	$this->SetSubjectEnd($tokens[9]);
	$this->SetPercentId($tokens[10]);
	$this->SetPercentSimilarity($tokens[11]);
	$this->SetRawScore($tokens[12]);
	$this->SetBitScore($tokens[13]);
	$this->SetDescription($tokens[15]);
	$this->SetFrame($tokens[16]);
	$this->SetStrand($tokens[17]);
	$this->SetSubjectLength($tokens[18]);
	$this->SetEValue($tokens[19]);
	$this->SetPValue($tokens[20]);

	if ($this->GetQueryBegin > $this->GetQueryEnd) {
		my $tmp = $this->GetQueryBegin;
		$this->SetQueryBegin($this->GetQueryEnd);
		$this->SetQueryEnd($tmp);
	}
	if ($this->GetSubjectBegin > $this->GetSubjectEnd) {
		my $tmp = $this->GetSubjectBegin;
		$this->SetSubjectBegin($this->GetSubjectEnd);
		$this->SetSubjectEnd($tmp);
	}
}

sub GetQueryName
{
	my $this = shift;
	return defined $this->{QUERY_NAME} ? $this->{QUERY_NAME} : "";
}

sub SetQueryName
{
	my ($this, $query_name) = @_;
	$this->{QUERY_NAME} = $query_name;
}

sub GetSubjectName
{
	my $this = shift;
	return defined $this->{SUBJ_NAME} ? $this->{SUBJ_NAME} : "";
}

sub SetSubjectName
{
	my ($this, $subj_name) = @_;
	$this->{SUBJ_NAME} = $subj_name;
}

sub GetQueryBegin
{
	my $this = shift;
	return defined $this->{QUERY_BEGIN} ? $this->{QUERY_BEGIN} : -1;
}

sub SetQueryBegin
{
	my ($this, $query_begin) = @_;
	$this->{QUERY_BEGIN} = $query_begin;
}

sub GetQueryEnd
{
	my $this = shift;
	return defined $this->{QUERY_END} ? $this->{QUERY_END} : -1;
}

sub SetQueryEnd
{
	my ($this, $query_end) = @_;
	$this->{QUERY_END} = $query_end;
}

sub GetQueryLength
{
	my $this = shift;
        return defined $this->{QUERY_LEN} ? $this->{QUERY_LEN} : -1;
}

sub SetQueryLength
{
	my ($this, $query_len) = @_;
	$this->{QUERY_LEN} = $query_len;
}

sub GetSubjectBegin
{
	my $this = shift;
        return defined $this->{SUBJ_BEGIN} ? $this->{SUBJ_BEGIN} : -1;
}

sub SetSubjectBegin
{
	my ($this, $subj_begin) = @_;
	$this->{SUBJ_BEGIN} = $subj_begin;
}

sub GetSubjectEnd
{
	my $this = shift;
        return defined $this->{SUBJ_END} ? $this->{SUBJ_END} : -1;
}

sub SetSubjectEnd
{
	my ($this, $subj_end) = @_;
	$this->{SUBJ_END} = $subj_end;
}

sub GetSubjectLength
{
	my $this = shift;
        return defined $this->{SUBJ_LEN} ? $this->{SUBJ_LEN} : -1;
}

sub SetSubjectLength
{
	my ($this, $subj_len) = @_;
	$this->{SUBJ_LEN} = $subj_len;
}

sub GetPercentId
{
	my $this = shift;
        return defined $this->{PCT_ID} ? $this->{PCT_ID} : -1;
}

sub SetPercentId
{
	my ($this, $pct_id) = @_;
	$this->{PCT_ID} = $pct_id;
}

sub GetPercentSimilarity
{
	my $this = shift;
        return defined $this->{PCT_SIM} ? $this->{PCT_SIM} : -1;
}

sub SetPercentSimilarity
{
	my ($this, $pct_sim) = @_;
	$this->{PCT_SIM} = $pct_sim;
}

sub GetRawScore
{
	my $this = shift;
        return defined $this->{RAW_SCORE} ? $this->{RAW_SCORE} : -1;
}

sub SetRawScore
{
	my ($this, $raw_score) = @_;
	$this->{RAW_SCORE} = $raw_score;
}

sub GetBitScore
{
	my $this = shift;
        return defined $this->{BIT_SCORE} ? $this->{BIT_SCORE} : -1;
}

sub SetBitScore
{
	my ($this, $bit_score) = @_;
	$this->{BIT_SCORE} = $bit_score;
}

sub GetDescription
{
	my $this = shift;
        return defined $this->{DESCR} ? $this->{DESCR} : "";
}

sub SetDescription
{
	my ($this, $descr) = @_;
	$this->{DESCR} = $descr;
}

sub GetFrame
{
	my $this = shift;
        return defined $this->{FRAME} ? $this->{FRAME} : -1;
}

sub SetFrame
{
	my ($this, $frame) = @_;
	$this->{FRAME} = $frame;
}

sub GetStrand
{
	my $this = shift;
        return defined $this->{STRAND} ? $this->{STRAND} : "";
}

sub SetStrand
{
	my ($this, $strand) = @_;
	$this->{STRAND} = $strand;
}

sub GetEValue
{
	my $this = shift;
        return defined $this->{EVAL} ? $this->{EVAL} : -1;
}

sub SetEValue
{
	my ($this, $eval) = @_;
	$this->{EVAL} = $eval;
}

sub GetPValue
{
	my $this = shift;
        return defined $this->{PVAL} ? $this->{PVAL} : -1;
}

sub SetPValue
{
	my ($this, $pval) = @_;
	$this->{PVAL} = $pval;
}

sub GetDatabaseName
{
	my $this = shift;
        return defined $this->{DB_NAME} ? $this->{DB_NAME} : "";
}

sub SetDatabaseName
{
	my ($this, $db_name) = @_;
	$this->{DB_NAME} = $db_name;
}

sub GetAlgorithm
{
	my $this = shift;
        return defined $this->{ALGO} ? $this->{ALGO} : "";
}

sub SetAlgorithm
{
	my ($this, $algo) = @_;
	$this->{ALGO} = $algo;
}

sub GetDate
{
	my $this = shift;
        return defined $this->{DATE} ? $this->{DATE} : "";
}

sub SetDate
{
	my ($this, $date) = @_;
	$this->{DATE} = $date;
}

sub ToString
{
	my $this = shift;
	return sprintf "%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t" .
		       "%d\t%.2f\t\t%s\t%s\t%s\t%d\t%.2e\t%.2e",
		       $this->GetQueryName(), $this->GetDate(),
		       $this->GetQueryLength(), $this->GetAlgorithm(),
		       $this->GetDatabaseName(), $this->GetSubjectName(),
		       $this->GetQueryBegin(), $this->GetQueryEnd(),
		       $this->GetSubjectBegin(), $this->GetSubjectEnd(),
		       $this->GetPercentId(),
		       $this->GetPercentSimilarity(),
		       $this->GetRawScore(), $this->GetBitScore(),
		       $this->GetDescription(), $this->GetFrame(),
		       $this->GetStrand(), $this->GetSubjectLength(),
		       $this->GetEValue(), $this->GetPValue();
}

1;
