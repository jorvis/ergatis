package BPlite;
use strict;
use overload '""' => '_overload';
###############################################################################
# BPlite
###############################################################################
sub new {
	my ($class, $fh) = @_;
	if (ref $fh !~ /GLOB/)
		{die "BPlite error: new expects a GLOB reference not $fh\n"}
	my $this = bless {};
	$this->{FH} = $fh;
	$this->{LASTLINE} = "";
	if ($this->_parseHeader) {
		# there are alignments
		$this->{REPORT_DONE} = 0;
		die if not defined $this->{QUERY};
		die if not defined $this->{DATABASE};
	}
	else {
		# empty report
		$this->{REPORT_DONE} = 1; 
	}
	return $this;
}
sub query       {shift->{QUERY}}
sub queryLength {shift->{QUERY_LENGTH}}
sub database    {shift->{DATABASE}}
sub nextSbjct {
	my ($this) = @_;
	$this->_fastForward or return 0;
	
	#######################
	# get all sbjct lines #
	#######################
	if ($this->{LASTLINE} =~ /^Parameters:|^\s+Database:/) {
		return 0;
	}
	my $def = $this->{LASTLINE};
	my $FH = $this->{FH};
	while(<$FH>) {
		if    ($_ !~ /\w/)            {next}
		elsif ($_ =~ /Strand HSP/)    {next} # WU-BLAST non-data
		elsif ($_ =~ /^\s{0,2}Score/) {$this->{LASTLINE} = $_; last}
		else                          {$def .= $_}
	}
	return 0 unless $def =~ /^>/;
	$def =~ s/\s+/ /g;
	$def =~ s/\s+$//g;
	my ($sbjct_length) = $def =~ /Length = ([\d,]+)$/;
	$sbjct_length =~ s/,//g;
	$def =~ s/Length = [\d,]+$//g;
	
	####################
	# the Sbjct object #
	####################
	my $sbjct = BPlite::Sbjct::new($def, $sbjct_length,
		$this->{FH}, $this->{LASTLINE}, $this);
	return $sbjct;
}
sub _parseHeader {
	my ($this) = @_;
	my $FH = $this->{FH};
	while(<$FH>) {
		if ($_ =~ /^Query=\s+(.+)/)    {
			my $query = $1;
			while(<$FH>) {
				last if $_ !~ /\S/;
				$query .= $_;
			}
			$query =~ s/\s+/ /g;
			$this->{QUERY} = $query;
			($this->{QUERY_LENGTH}) = $query =~ /([\d,]+) letters/;
			$this->{QUERY_LENGTH} =~ s/\D//g;
		}
		elsif ($_ =~ /^Database:\s+(.+)/) {$this->{DATABASE} = $1}
		elsif ($_ =~ /^>/)                {$this->{LASTLINE} = $_; return 1}
		elsif ($_ =~ /^Parameters|^\s+Database:/) {
			$this->{LASTLINE} = $_;
			$this->{REPORT_DONE} = 1;
			return 0; # there's nothing in the report
		}
	}
}
sub _fastForward {
	my ($this) = @_;
	return 0 if $this->{REPORT_DONE};
	return 1 if $this->{LASTLINE} =~ /^>/;
	if ($this->{LASTLINE} =~ /^Parameters|^\s+Database:/) {
		$this->{REPORT_DONE} = 1;
		return 1;
	}
	my $FH = $this->{FH};
	while(<$FH>) {
		die if /^BLAST/;
		if ($_ =~ /^>|^Parameters|^\s+Database:/) {
			$this->{LASTLINE} = $_;
			$this->{REPORT_DONE} = 1;
			return 1;
		}
	}
	warn "Possible parse error in _fastForward in BPlite.pm\n";
}
sub _overload {
	my ($this) = @_;
	return $this->{QUERY} . " vs. " . $this->{DATABASE};
}

###############################################################################
# BPlite::Sbjct
###############################################################################
package BPlite::Sbjct;
use overload '""' => 'name';
sub new {
	my $sbjct = bless {};
	($sbjct->{NAME}, $sbjct->{LENGTH}, $sbjct->{FH},$sbjct->{LASTLINE},
		$sbjct->{PARENT}) = @_;
	$sbjct->{HSP_ALL_PARSED} = 0;
	return $sbjct;
}
sub name {shift->{NAME}}
sub length {shift->{LENGTH}}
sub nextHSP {
	my ($sbjct) = @_;
	return 0 if $sbjct->{HSP_ALL_PARSED};
	
	############################
	# get and parse scorelines #
	############################
	my $scoreline = $sbjct->{LASTLINE};
	my $FH = $sbjct->{FH};
	my $nextline = <$FH>;
	return undef if not defined $nextline;
	$scoreline .= $nextline;
	my ($score, $bits);
	if ($scoreline =~ /\d bits\)/) {
		($score, $bits) = $scoreline =~
			/Score = (\d+) \((\S+) bits\)/; # WU-BLAST
	}
	else {
		($bits, $score) = $scoreline =~
			/Score =\s+(\S+) bits \((\d+)/; # NCBI-BLAST
	}
	
	my ($match, $length) = $scoreline =~ /Identities = (\d+)\/(\d+)/;
	my ($positive) = $scoreline =~ /Positives = (\d+)/;
	$positive = $match if not defined $positive;
	my ($p)       = $scoreline =~ /[Sum ]*P[\(\d+\)]* = (\S+)/;
	if (not defined $p) {
		($p) = $scoreline =~ /Expect\S* = (\S+)/; # thanks Keith Allen
	}
	$p =~ s/,//g;
	
	die "parse error $scoreline\n" if not defined $score;

	#######################
	# get alignment lines #
	#######################
	my @hspline;
	while(<$FH>) {
		die if /^BLAST/;
		if ($_ =~ /^WARNING:|^NOTE:|^ERROR:|^FATAL:/) {
			while(<$FH>) {last if $_ !~ /\S/}
		}
		elsif ($_ !~ /\S/)            {next}
		elsif ($_ =~ /^\s*Score/)     {$sbjct->{LASTLINE} = $_; last}
		elsif ($_ =~ /Strand HSP/)    {next} # not capturing
		elsif ($_ =~ /^\s*Frame/)     {next} # not capturing
		elsif ($_ =~ /^\s*Links/)     {next} # not capturing
		elsif ($_ =~ /^\s*Strand/)    {next} # not capturing
		elsif ($_ =~ /^EXIT CODE/) {last}
		elsif ($_ =~ /^>|^Parameters|^\s+Database:/)   {
			$sbjct->{LASTLINE} = $_;
			$sbjct->{PARENT}->{LASTLINE} = $_;
			$sbjct->{HSP_ALL_PARSED} = 1;
			last;
		}
		elsif ($_ =~ /^Query:/) {
			push @hspline, $_;           # store the query line
			my $l1 = <$FH>;              # either alignment line or sbjct line
			if ($l1 =~ /^Sbjct/) {
				push @hspline, "";  # dummy line, this is a -noseq option
				push @hspline, $l1; # so store a fake alignment and real sbjct
				next;
			}
			push @hspline, $l1;                 # grab/store the alignment line
			my $l2 = <$FH>; push @hspline, $l2; # grab/store the sbjct line
		}
	}
	
	#########################
	# parse alignment lines #
	#########################
	my ($ql, $sl, $as) = ("", "", "");
	my ($qb, $qe, $sb, $se) = (0,0,0,0);
	my (@QL, @SL, @AS); # for better memory management
			
	for(my $i=0;$i<@hspline;$i+=3) {
		#warn $hspline[$i], $hspline[$i+2];
		$hspline[$i]   =~ /^Query:\s+(\d+)\s+(\S+)\s+(\d+)/;
		return 0 unless defined $3;
		$ql = $2; $qb = $1 unless $qb; $qe = $3;
		
		my $offset = index($hspline[$i], $ql);
		$as = substr($hspline[$i+1], $offset, CORE::length($ql))
			if $hspline[$i+1];
		
		$hspline[$i+2] =~ /^Sbjct:\s+(\d+)\s+(\S+)\s+(\d+)/;
		return 0 unless defined $3;
		$sl = $2; $sb = $1 unless $sb; $se = $3;
		
		push @QL, $ql; push @SL, $sl; push @AS, $as;
	}

	##################
	# the HSP object #
	##################
	$ql = join("", @QL);
	$sl = join("", @SL);
	$as = join("", @AS);
	my $qgaps = $ql =~ tr/-/-/;
	my $sgaps = $sl =~ tr/-/-/;
	my $hsp = BPlite::HSP::new( $score, $bits, $match, $positive, $length, $p,
		$qb, $qe, $sb, $se, $ql, $sl, $as, $qgaps, $sgaps);
	return $hsp;
}

###############################################################################
# BPlite::HSP
###############################################################################
package BPlite::HSP;
use overload '""' => '_overload';
sub new {
	my $hsp = bless {};
	($hsp->{SCORE}, $hsp->{BITS},
		$hsp->{MATCH}, $hsp->{POSITIVE}, $hsp->{LENGTH},$hsp->{P},
		$hsp->{QB}, $hsp->{QE}, $hsp->{SB}, $hsp->{SE},
		$hsp->{QL}, $hsp->{SL}, $hsp->{AS}, $hsp->{QG}, $hsp->{SG}) = @_;
	$hsp->{PERCENT} = int(1000 * $hsp->{MATCH}/$hsp->{LENGTH})/10;
	return $hsp;
}
sub _overload {
	my $hsp = shift;
	return $hsp->queryBegin."..".$hsp->queryEnd." ".$hsp->bits;
}
sub score           {shift->{SCORE}}
sub bits            {shift->{BITS}}
sub percent         {shift->{PERCENT}}
sub match           {shift->{MATCH}}
sub positive        {shift->{POSITIVE}}
sub length          {shift->{LENGTH}}
sub P               {shift->{P}}
sub queryBegin      {shift->{QB}}
sub queryEnd        {shift->{QE}}
sub sbjctBegin      {shift->{SB}}
sub sbjctEnd        {shift->{SE}}
sub queryAlignment  {shift->{QL}}
sub sbjctAlignment  {shift->{SL}}
sub alignmentString {shift->{AS}}
sub queryGaps       {shift->{QG}}
sub sbjctGaps       {shift->{SG}}
sub qb              {shift->{QB}}
sub qe              {shift->{QE}}
sub sb              {shift->{SB}}
sub se              {shift->{SE}}
sub qa              {shift->{QL}}
sub sa              {shift->{SL}}
sub as              {shift->{AS}}
sub qg              {shift->{QG}}
sub sg              {shift->{SG}}


###############################################################################
# BPlite::Multi
###############################################################################
package BPlite::Multi;
sub new {
	my ($class, $fh) = @_;
	if (ref $fh !~ /GLOB/)
		{die "BPlite error: new expects a GLOB reference not $fh\n"}
	my $this = bless {};
	$this->{FH} = $fh;
	return $this;
}
sub nextReport {
	my ($this) = @_;
	my $FH = $this->{FH};
	while(<$FH>) {
		last if /^T?BLAST[NPX]/;
	}
	return 0 if not defined $_;
	my $blast = new BPlite($this->{FH});
	return $blast;
}

1;
__END__

=head1 NAME

BPlite - Lightweight BLAST parser

=head1 SYNOPSIS

 use BPlite;
 
 # single BLAST report
 
 my $report = new BPlite(\*STDIN);
 $report->query;
 $report->queryLength;
 $report->database;
 while(my $sbjct = $report->nextSbjct) {
     $sbjct->name;
     $sbjct->length;
     while (my $hsp = $sbjct->nextHSP) {
         $hsp->score;
         $hsp->bits;
         $hsp->percent;
         $hsp->P;
         $hsp->match;
         $hsp->positive;
         $hsp->length;
         $hsp->queryBegin;
         $hsp->queryEnd;
         $hsp->sbjctBegin;
         $hsp->sbjctEnd;
         $hsp->queryAlignment;
         $hsp->sbjctAlignment;
         $hsp->alignmentString;
         $hsp->queryGaps;
         $hsp->sbjctGaps;
     }
 }
 
 # multiple (concatenated) BLAST reports
 
 my $multiple_report = new BPlite::Multi(\*STDIN);
 while(my $blast = $multiple_report->nextReport) {
     while(my $sbjct = $blast->nextSbjct) {
         while(my $hsp = $sbjct->nextHSP) {
         }
     }
 }


=head1 DESCRIPTION

BPlite is a package for parsing BLAST reports. The BLAST programs are a family
of widely used algorithms for sequence database searches. The reports are
non-trivial to parse, and there are differences in the formats of the various
flavors of BLAST. BPlite parses BLASTN, BLASTP, BLASTX, TBLASTN, and TBLASTX
reports from both the high performance WU-BLAST, and the more generic
NCBI-BLAST.

Many people have developed BLAST parsers (I myself have made at least three).
BPlite is for those people who would rather not have a giant object
specification, but rather a simple handle to a BLAST report that works well
in pipes.

=head2 Object

BPlite has three kinds of objects, the report, the subject, and the HSP. To
create a new report, you pass a filehandle reference to the BPlite constructor.

 my $report = new BPlite(\*STDIN); # or any other filehandle

The report has two attributes (query and database), and one method (nextSbjct).

 $report->query;     # access to the query name
 $report->queryLength; 
 $report->database;  # access to the database name
 $report->nextSbjct; # gets the next subject
 while(my $sbjct = $report->nextSbjct) {
     # canonical form of use is in a while loop
 }

A subject is a BLAST hit, which should not be confused with an HSP (below). A
BLAST hit may have several alignments associated with it. A useful way of
thinking about it is that a subject is a gene and HSPs are the exons. Subjects
have one attribute (name) and one method (nextHSP).

 $sbjct->name;    # access to the subject name
 "$sbjct";        # overloaded to return name
 $sbjct->nextHSP; # gets the next HSP from the sbjct
 while(my $hsp = $sbjct->nextHSP) {
     # canonical form is again a while loop
 }

An HSP is a high scoring pair, or simply an alignment. HSP objects do not have
any methods, just attributes (score, bits, percent, P, match, positive, length,
queryBegin, queryEnd, sbjctBegin, sbjctEnd, queryAlignment, sbjctAlignment)
that should be familiar to anyone who has seen a blast report. For
lazy/efficient coders, two-letter abbreviations are available for the
attributes with long names (qb, qe, sb, se, qa, sa).

 $hsp->score;
 $hsp->bits;
 $hsp->percent;
 $hsp->P;
 $hsp->match;
 $hsp->positive;
 $hsp->length;
 $hsp->queryBegin;      $hsp->qb;
 $hsp->queryEnd;        $hsp->qe;
 $hsp->sbjctBegin;      $hsp->sb;
 $hsp->sbjctEnd;        $hsp->se;
 $hsp->queryAlignment;  $hsp->qa;
 $hsp->sbjctAlignment;  $hsp->sa;
 $hsp->alignmentString; $hsp->as;
 $hsp->queryGaps;       $hsp->qg;
 $hsp->sbjctGaps;       $hsp->sg;
 "$hsp"; # overloaded for qb..qe bits

I've included a little bit of overloading for double quote variable
interpolation convenience. A report will return the query and database, a
subject will return its name and an HSP will return its queryBegin, queryEnd,
and bits in the alignment. Feel free to modify this to whatever is most
frequently used by you.

So a very simple look into a BLAST report might look like this.

 my $report = new BPlite(\*STDIN);
 while(my $sbjct = $report->nextSbjct) {
     print "$sbjct\n";
     while(my $hsp = $sbjct->nextHSP) {
	 	print "\t$hsp\n";
     }
 }

The output of such code might look like this:

 >foo
     100..155 29.5
     268..300 20.1
 >bar
     100..153 28.5
     265..290 22.1

=head2 Concatenated BLAST reports

You can step through multiple BLAST reports if you have a file of concatenated
BLAST reports using the following construct.

 my $multiple_report = new BPlite::Multi(\*STDIN);
 while(my $blast = $multiple_report->nextReport) {
     while(my $sbjct = $blast->nextSbjct) {
         while(my $hsp = $sbjct->nextHSP) {
         }
     }
 }


=head1 AUTHOR

Ian Korf (ikorf@sapiens.wustl.edu, http://sapiens.wustl.edu/~ikorf)

=head1 ACKNOWLEDGEMENTS

This software was developed at the Genome Sequencing Center at Washington
Univeristy, St. Louis, MO.

=head1 COPYRIGHT

Copyright (C) 1999-2001 Ian Korf. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut










