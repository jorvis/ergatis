#!/usr/local/bin/perl

use strict;
use warnings;
use Data::Dumper;
use XML::Twig;
use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;
use IO::File;

my $in                  = "/dev/stdin";
my $out                 = *STDOUT;
my $output_mrna_feats   = 0;
my $extract_all_ec      = 0;
my $percent_n_cutoff    = 10;
my $pid_feattype        = ["transcript"];
my $pid_field           = "id";
my $gene_feattype       = ["gene"];
my $gene_field          = "locus";
my $asmbl_seq = '';
my $asmbl_length;
## genes with equal or more than this % of Ns will be skiped (warning printed)


parse_options();
convert();

sub print_usage
{
    my $progname = basename($0);
    die << "END";
usage: $progname [--input|-i <input_bsml>]
        [--output|-o <output_tbl>] [--mrna|-m] [--ec_all|-e] [-h]

        -m: export mRNA features [ default - 0 (false) ]
        -e: extract all ec numbers (including incomplete e.g. #.#.#.-)
            [ default - false ]
END
}

sub parse_options
{
    my %opts = ();
    GetOptions(\%opts, 
               "input|i=s", 
               "output|o=s", 
               "pid_feattype=s",
               "pid_field=s",
               "gene_feattype:s",
               "gene_field:s",
               "mrna|m=i", 
               "help|h");
    print_usage() if $opts{help};
    $in = $opts{input} if $opts{input};
    $out = new IO::File($opts{output}, "w") or
        die "Error writing tbl to $opts{output}: $!"
        if $opts{output};
    $output_mrna_feats = 0 unless $opts{mrna};
    $extract_all_ec = 1 if $opts{ec_all};
    if ( defined $opts{percent_n_cutoff} ) {
        $percent_n_cutoff = $opts{percent_n_cutoff};
    }
    if ( defined $opts{pid_feattype} ) {
        @{$pid_feattype} = split(/,/,$opts{pid_feattype});
    }
    if( defined $opts{pid_field} ) {
        $pid_field = $opts{pid_field};
    }
    if ( defined $opts{gene_feattype} ) {
        @{$gene_feattype} = split(/,/,$opts{gene_feattype});
    }
    if( defined $opts{gene_field} ) {
        $gene_field = $opts{gene_field};
    }


}

sub convert
{
    my $twig = new XML::Twig();
    my %feats = ();

    my $output_features = [];
    my $organism_name;
    $twig->setTwigRoots({'Organism' => sub {
                         my ($twig,$elt) = @_;
                         $organism_name = $elt->att('genus')." ".$elt->att('species');
                         },
                         'Feature' => sub { process_feature(\%feats, @_); },
                         'Sequence' => \&process_seq });
    $twig->parsefile($in);
    $twig = new XML::Twig();
    $twig->setTwigRoots({'Feature-group' =>
            sub { process_feature_group(\%feats, $output_features, @_); } });
    $twig->parsefile($in);
    process_misc_feats(\%feats);

    my @srted_feats = sort {$a->{'start'} <=> $b->{'start'}} @$output_features;
    

    # First print out headers
    print $out "$organism_name\n";
    print $out scalar(@srted_feats)." proteins\n";
    print $out "Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n";
    my $seen_coords = {};
    foreach my $feat (@srted_feats) {
        if(!$seen_coords->{"$feat->{'start'}..$feat->{'stop'}"}) {
            print $out join("\t", ("$feat->{'start'}..$feat->{'stop'}",$feat->{'strand'},$feat->{'length'},$feat->{'pid'},$feat->{'gene'},'-','-','-',$feat->{'gene_product'}))."\n";
            $seen_coords->{"$feat->{'start'}..$feat->{'stop'}"} = 1;
        }
        else {
            print STDERR "Found duplicate genes at $feat->{'start'}..$feat->{'stop'}\n";
        }
    }
}

sub process_seq
{
    my ($twig, $elt) = @_;
    return if $elt->att('class') ne "assembly";
#    print $out ">Features ", $elt->att('id'), "\n";
    
    ## get the sequence length
    if ( $elt->has_child('Seq-data-import') ) {
        open(my $ifh, "<" . $elt->first_child('Seq-data-import')->att('source')) || die "can't open seq-data-import: $!";
        while (<$ifh>) {
            if (! /^\>/ ) {
                $asmbl_seq .= $_;
            }
        }
        $asmbl_seq =~ s/\s+//g;
        $asmbl_length = length $asmbl_seq;
    } else {
        die "error: sequence element has no seq-data-import\n";
    }
}

sub process_feature
{
    my ($feats, $twig, $elt) = @_;
    my $id = $elt->att('id');
    $feats->{$id} = $elt;
    $twig->purge();
}

sub process_feature_group
{
    my ($feats, $output_features, $twig, $elt) = @_;
    my @exons = ();
    my $cds = undef;
    my $class = undef;
    my $gene = undef;
    my $polypeptide = undef;
    my $transcript = undef;
    my @repeats = ();
    my $group_by_class = {};

    foreach my $feat_member ($elt->children('Feature-group-member')) {
        my $feat_type = $feat_member->att('feature-type');
        my $featref = $feat_member->att('featref');
        my $feat = $feats->{$featref};
        if(!$group_by_class->{$feat_type}) {
            $group_by_class->{$feat_type} = [];
        }
        push(@{$group_by_class->{$feat_type}}, $feat);
        
        if ($feat_type eq "exon") {
            push @exons, $feat;
        }
        elsif ($feat_type eq "transcript") {
            $class = "transcript";
            $transcript = $feat;
        }
        elsif ($feat_type eq "tRNA") {
            $class = "tRNA";
            $transcript = $feat;
        }
        elsif ($feat_type eq "ncRNA") {
            $class = "ncRNA";
        }
        elsif ($feat_type eq "gene") {
            $gene = $feat;
        }
        elsif ($feat_type eq "CDS") {
            $cds = $feat;
        }
        elsif ($feat_type eq "polypeptide") {
            $polypeptide = $feat;
        }
        elsif ($feat_type eq "repeat_region") {
            push @repeats, $feat;
        } elsif( $feat_type eq "signal_peptide" ) {
            $class = "signal_peptide";
        }
    }


    if (!$transcript) {
        return;
    }
    # Going to pull the location info here
    my $location = $transcript->first_child('Interval-loc');
    my $fmin = $location->att('startpos');
    my $fmax = $location->att('endpos');    
    my $strand = $location->att('complement') ==0 ? '+' : '-';
    my $aa_length = int(($fmax-$fmin)/3); # HACK - obviously this is not always right. The issue is that .ptt techincally wants amino acid length.

    # Pulling the gene product off of the transcript.
    my $gene_prod_feat_att = $transcript->first_child('Attribute[@name="gene_product_name"]');
    my $gene_product = '-';
    if($gene_prod_feat_att) {
        $gene_product = $gene_prod_feat_att->att('content');
    }

    my $gene_val = '-';
    foreach my $type (@$gene_feattype) {
        my $gene_feat = $group_by_class->{$type}->[0];
        if($gene_feat) {
            # Need to be able to handle additional fields here.
            if($gene_field eq 'id') {
                $gene = $gene_feat->att('id');
            }
            elsif($gene_field eq 'locus') {
                my @xrefs = $gene_feat->children('Cross-reference');
                map {
                    if($_->att('identifier-type') eq $gene_field) {
                        $gene_val = $_->att('identifier');
                    }
                }@xrefs;
            }
        }
        last if $gene_val ne '-';
    }

    # This is a HACK of course cause if there are more than one of these types we will have issues
    my $pid;

    foreach my $type (@$pid_feattype) {
        my $pid_feat = $group_by_class->{$type}->[0];
        
        if($pid_feat) {
            # Need to be able to handle additional fields here.
            if($pid_field eq 'id') {
                $pid = $pid_feat->att('id');
            }
        }
        last if $pid;
    }
    

    if($pid) {
        push(@$output_features, {'start'        => $fmin,
                                 'stop'         => $fmax,
                                 'strand'       => $strand,
                                 'length'       => $aa_length,
                                 'pid'          => $pid,
                                 'gene'         => $gene_val,
                                 'gene_product' => $gene_product
             });
    }
}

sub process_gene
{
    my ($gene, $transcript) = @_;
    my $locus = $gene->first_child('Cross-reference[@database="TIGR_moore"]');
    my $locus_tag = $gene->first_child('Attribute[@name="locus_tag"]');
    my $cross_reference = $gene->first_child('Cross-reference[@identifier-type="locus"]');
    my %attrs = ();
    $attrs{'locus_tag'} = $locus_tag->att('content') if $locus_tag;
    $attrs{'locus_tag'} = $locus->att('identifier') if $locus;
    $attrs{'locus_tag'} = $cross_reference->att('identifier') if $cross_reference;
    if ($transcript) {
        my $symbol =
            $transcript->first_child('Attribute[@name="gene_symbol"]');
        $attrs{'gene'} = $symbol->att('content') if $symbol;
    }
    print_feat([$gene], "gene", \%attrs);
}

sub process_cds
{
    my ($exons, $cds, $attrs) = @_;
    my @cdss;

    if (!defined($cds)){
	die "cds was not defined";
    }

    my $cds_loc = $cds->first_child('Interval-loc');
    my $cds_from = $cds_loc->att('startpos');
    my $cds_to = $cds_loc->att('endpos');
    foreach my $exon (@{$exons}) {
        my $exon_loc = $exon->first_child('Interval-loc');
        my $exon_from = $exon_loc->att('startpos');
        my $exon_to = $exon_loc->att('endpos');
        next if ($exon_to < $cds_from || $exon_from > $cds_to);
        my $elt = $exon->copy();
        if ($exon_from < $cds_from) {
            $elt->first_child('Interval-loc')->set_att('startpos', $cds_from);
        }
        if ($exon_to > $cds_to) {
            $elt->first_child('Interval-loc')->set_att('endpos', $cds_to);
        }
        push @cdss, $elt;
    }
    my $cog = $cds->first_child('Attribute[@name="top_cog_hit"]');
    if ($cog) {
        $attrs->{note} = $cog->att('content');
    }
    print_feat(\@cdss, "CDS", $attrs);
}

sub process_misc_feats
{
    my ($feats) = @_;
    while (my ($key, $val) = each %{$feats}) {
        my $class = $val->att('class');
        if ($class eq "repeat_region") {
            print_feat([$val], "repeat_region");
        }
    }
}

sub print_feat
{
    my ($feats, $class, $attrs) = @_;
    for (my $i = 0; $i < scalar(@{$feats}); ++$i) {
        my $feat = $feats->[$i];
        my $loc = $feat->first_child('Interval-loc');
        my $from = $loc->att('startpos') + 1;
        my $to = $loc->att('endpos');
        my $codon_start = 0;
        if ($loc->att('complement') == 1) {
            ($from, $to) = ($to, $from);
        }
        my $tmp_from = $from;
        my $tmp_to = $to;
        
        ## check to see if we're off the end of the sequence here
        #   I'm addressing the other cases as I see them.
        if ( $from > $asmbl_length && $from > $to ) {
            $tmp_from = "<$asmbl_length";
            $attrs->{codon_start} = 4 - ( $from - $asmbl_length );
        } elsif( $to > $asmbl_length && $to > $from ) {
            $tmp_to = ">$asmbl_length";
        }

        ## Deal with negative coordinates
        if ($from <= 0 && $from < $to) {
            if ($class eq "CDS") {
                my $frame = abs($from - 1) % 3;
                if ($frame == 1) {
                    $attrs->{codon_start} = 3;
                }
                elsif ($frame == 2) {
                    $attrs->{codon_start} = 2;
                }
            }
            $tmp_from = "<1";
        }

        if($to <= 0) {
            $tmp_to = ">1";
        }
       
        $from = $tmp_from;
        $to = $tmp_to;
        
        if ( $class eq 'gene' || $class eq 'CDS' || $class eq 'mRNA' ) {
            ## bail if this sequence has too many Ns

	    my $startpos =  $loc->att('startpos');
	    my $endpos = $loc->att('endpos');
	    
	    if ($startpos == $endpos){
		die "startpos == endpos ('$startpos') for <Feature> ".
		"with class '$class'";
	    }

            my $feat_seq = substr($asmbl_seq, $startpos, $endpos - $startpos );

	    my $ofeat_seq = $feat_seq;

	    my $feat_seq_length = length($feat_seq);

	    if ($feat_seq_length == 0){

		my $asmbl_seq_length = length($asmbl_seq);
		print STDERR "loc:". Dumper $loc;
		print STDERR "The Feature with class '$class' had feat_seq_length == '0'. ".
		"startpos '$startpos' endpos '$endpos' asmbl_seq_length ".
		"'$asmbl_seq_length' asmbl_length '$asmbl_length' original_feat_seq ".
		"'$ofeat_seq'";
		die;
	    }

            my $n_count = ( $feat_seq =~ tr/Nn// );
            my $n_perc = ($n_count / $feat_seq_length) * 100;

            if ( $n_perc >= $percent_n_cutoff ) {
                print STDERR "warning: skipping $class because it has too many Ns (", sprintf("%.1f",$n_perc), "\%)\n";
                return;
            }
        }
        
        print $out "$from\t$to";
        print $out "\t$class" if $i == 0;
        print $out "\n";
    }

    while (my ($key, $val) = each %{$attrs}) {
        if (ref($val) eq "ARRAY") {
            foreach my $v (@{$val}) {
                print_attribute($out, $key, $v);
            }
        }
        else {
            print_attribute($out, $key, $val);
        }
    }
}

sub print_attribute
{
    my ($out, $key, $val) = @_;
    print $out "\t\t\t$key\t$val\n";
}
