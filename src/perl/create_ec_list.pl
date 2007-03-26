use strict;
use warnings;

use Blast::NcbiBlastHitDataType;

use File::Basename;
use IO::File;
use Getopt::Long qw(:config no_ignore_case bundling);

my %model_tu_map    = ();
my @hits_files       = ();
my $out             = *STDOUT;
my $max_eval        = 0.1;
my $rps             = 0;

&parse_options;
&create_list;

sub print_usage
{
	my $prog_name = basename($0);
	die << "END";
usage: $prog_name [--model_tu_map|-m <model_id_to_tu_id_map>]
        [ [--hits|-i <hits>] | [--hits_list|-l <hits_list>] ]
        [--output|-o <output>]
        [--max_eval|-e <max_eval_cutoff>] [--rps|-r] [--help|-h]

        --rps|-r:   hits come from rpsblast rather than blastpgp
                    [default = false]
        --model_tu_map|-m:  tab delimited list of model id to tu id
END
}

sub parse_options
{
	my %opts = ();
	GetOptions(\%opts, "model_tu_map|m=s", "hits|i=s", "hits_list|l=s",
            "output|o=s", "max_eval|e=f", "rps|r", "help|h");
	&print_usage if $opts{help};
	&init_map(\%model_tu_map, $opts{model_tu_map}) if $opts{model_tu_map};
	$max_eval = $opts{max_eval} if $opts{max_eval};
    parse_hits_list($opts{hits_list}) if $opts{hits_list};
	push @hits_files, $opts{hits} if $opts{hits};
	$out = new IO::File($opts{output}, "w") or
		die "Error writing to output file $opts{output}: $!"
		if ($opts{output});
	$rps = 1 if $opts{rps};
    if (!scalar(@hits_files)) {
        push @hits_files, "/dev/stdin";
    }
}

sub init_map
{
	my ($map, $fname) = @_;
	my $data = new IO::File($fname) or die "Error reading map data: $!";
	while (my $line = <$data>) {
		chomp $line;
		my @tokens = split /\t/, $line;
		$$map{$tokens[0]} = $tokens[1];
	}
}

sub parse_hits_list
{
    my ($hits_list) = @_;
    my $fh = new IO::File($hits_list) or die "Error reading hits list: $!";
    while (my $file = <$fh>) {
        chomp $file;
        push @hits_files, $file;
    }
}

sub create_list
{
    foreach my $hits_file (@hits_files) {
        my $hits = new IO::File($hits_file) or
            die "Error reading hits file $hits_file: $!";
        my %tus = ();
        while (my $line = <$hits>) {
            chomp $line;
            my $hit = new Blast::NcbiBlastHitDataType($line);
            my $profile_id = get_profile_id($hit);
            my $model = get_protein_id($hit);
            my $eval = $hit->GetEValue();
            my $bitscore = $hit->GetBitScore();
            next if $eval > $max_eval;
            my $tu = $model_tu_map{$model} ?
                $model_tu_map{$model} : $model;
            my $new_data = 1;
            foreach my $loaded_hit (@{$tus{$tu}}) {
                if ($profile_id eq get_profile_id($loaded_hit) &&
                        $eval >= $loaded_hit->GetEValue()) {
                    $new_data = 0;
                    last;
                }
            }
            push @{$tus{$tu}}, $hit if $new_data;
        }
        while (my ($tu, $hits) = each %tus) {

            my @written = ();
            my %written_ecs = ();

            foreach my $hit (sort hit_comparator @$hits) {
                if (is_valid($hit, \@written)) {
                    if (!$written_ecs{get_ec($hit)}) {
                        write_data($tu, $hit);
                        ++$written_ecs{get_ec($hit)};
                    }
                    push @written, $hit;
                }
            }
        }
    }
}

sub hit_comparator
{
	if ($a->GetEValue() == $b->GetEValue()) {
		return $b->GetBitScore() <=> $a->GetBitScore();
	}
	return $a->GetEValue() <=> $b->GetEValue();
}

sub overlap
{
	my ($hit1, $hit2) = @_;
	if (get_protein_id($hit1) ne get_protein_id($hit2)) {
		return 0;
	}
	my ($from1, $to1) = (get_protein_start($hit1),
		get_protein_end($hit1));
	my ($from2, $to2) = (get_protein_start($hit2),
		get_protein_end($hit2));
	($from1, $to1) = ($to1, $from1) if $from1 > $to1;
	($from2, $to2) = ($to2, $from2) if $from2 > $to2;
	if ($from1 < $from2) {
		return $to1 > $from2;
	}
	else {
		return $to2 > $from1;
	}
}

sub is_valid
{
	my ($hit, $written) = @_;
	foreach my $written_hit (@$written) {
		if (get_profile_id($hit) eq get_profile_id($written_hit)) {
			return 0;
		}
		if (overlap($hit, $written_hit)) {
			return 0;
		}
	}
	return 1;
}

sub get_ec
{
	my $hit = shift;
	my $ec = get_profile_id($hit);
	$ec =~ s/^\d+p//;
	return $ec;
}

sub get_profile_id
{
	my $hit = shift;
	return $rps ? $hit->GetSubjectId() : $hit->GetQueryId();
}

sub get_protein_id
{
	my $hit = shift;
	return $rps ? $hit->GetQueryId() : $hit->GetSubjectId();
}

sub get_ec_start
{
	my $hit = shift;
	return $rps ? $hit->GetSubjectStart() : $hit->GetQueryStart();
}

sub get_ec_end
{
	my $hit = shift;
	return $rps ? $hit->GetSubjectEnd() : $hit->GetQueryEnd();
}

sub get_protein_start
{
	my $hit = shift;
	return $rps ? $hit->GetQueryStart() : $hit->GetSubjectStart();
}

sub get_protein_end
{
	my $hit = shift;
	return $rps ? $hit->GetQueryEnd() : $hit->GetSubjectEnd();
}

sub write_data
{
	my ($tu, $hit) = @_;
	$out->printf("%s\t%s\t%s\t%s\t%.2f\t%d\t%d\n",
		$tu, get_ec($hit), get_protein_id($hit), $hit->GetEValue(),
		$hit->GetBitScore(), get_protein_start($hit),
		get_protein_end($hit));
}
