#! /local/perl/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BSML::BsmlBuilder;
use Data::Dumper;
use File::Basename;

my %options = ();
my $results = GetOptions (\%options, 'ali_dir|a=s', 'output|o=s', 'verbose|v', 'help|h');

###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $output    = $options{'output'};
my $ali_dir   = $options{'ali_dir'};
$ali_dir =~ s/\/+$//;
my $molecule_type = $options{'mol_type'} || 'protein';

###-------------------------------------------------------###

my $builder = new BSML::BsmlBuilder;

opendir(DIR, $ali_dir) or die "Unable to access $ali_dir due to $!";
while( my $file = readdir(DIR)) {
    next if ($file =~ /^\.{1,2}$/);  #skip  "." ,  ".."
    if($file =~ /^(.+?)(\.msf\.dtab)?\.msf$/) {
	my $fam_name = $1;
	my $MSF_alignments = process_MSF_file("$ali_dir/$file");
	next if(keys %$MSF_alignments < 1);   #skip empty msf files
	my $table = $builder->createAndAddMultipleAlignmentTable('molecule-type' => $MSF_alignments->{'mol_type'});
	my $summary = $builder->createAndAddAlignmentSummary( 'multipleAlignmentTable' => $table,
							      'seq-type' =>  $MSF_alignments->{'mol_type'},
							      'seq-format' => 'MSF' 
							      );
	my $aln = $builder->createAndAddSequenceAlignment( 'multipleAlignmentTable' => $table );
	my $seqnum=0;
	my $sequences_tag;
	foreach my $seq (keys %{ $MSF_alignments->{'proteins'} }) {
	    $seqnum++;
	    my $alignment = join ('', @{ $MSF_alignments->{'proteins'}->{$seq}->{'alignment'} });
	    my $start  = $MSF_alignments->{'proteins'}->{$seq}->{'end5'}; 
	    my $end    = $MSF_alignments->{'proteins'}->{$seq}->{'end3'};
	    my $domain_length = $end - $start + 1;
	    $builder->createAndAddAlignedSequence( 'alignmentSummary' => $summary,
						   'seqnum' => $seqnum,
						   'length' => $domain_length,
						   'name'   => $seq,
                                                   'start'  => $start
                                                 );
	    $builder->createAndAddSequenceData( 'sequenceAlignment' => $aln,
						'seq-name' => $seq,
						'seq-data' => $alignment
                                              );
	    $sequences_tag .= "$seqnum:";
	}
	$aln->addattr( 'sequences', $sequences_tag );

    }

}

$builder->write( $output );





sub process_MSF_file {

    my $file = shift;


    my $MSF_alignments ={};
    open(MSF, "$file") or die "Unable to open $file due to $!";
    my $line;
    my $msf_type;
    while(defined($line = <MSF>) and $line !~ /^\/\//) {
	if( $line =~ /MSF:\s+([\S]+)\s+Type:\s+([\S]+)\s+Check/) {
	    my $msf_length = $1;
	    if($2 eq 'P') {
		$msf_type = 'protein';
	    }elsif($2 eq 'N') {
		$msf_type = 'nucleotide';
	    }else {
		$msf_type = 'protein';
	    }
	    $MSF_alignments->{'mol_type'} = $msf_type;
	}

	if($line =~ /Name:\s+([\S]+)\s+Len:\s+([\S]+)\s+Check:\s+([\S]+)\s+Weight:\s+([\S]+)/) {
	    my $name    = $1;
	    my $ali_len = $2;
	    my $check   = $3;
	    my $weight  = $4;
	    my ($protein, $coord) = split(/\//, $name);
	    my ($start5, $end3) = split("-", $coord);
	    
	    $MSF_alignments->{'proteins'}->{$protein}->{'length'} = $ali_len;
	    $MSF_alignments->{'proteins'}->{$protein}->{'check'}  = $check;
	    $MSF_alignments->{'proteins'}->{$protein}->{'weight'} = $weight;
	    $MSF_alignments->{'proteins'}->{$protein}->{'end5'}   = $start5;
	    $MSF_alignments->{'proteins'}->{$protein}->{'end3'}   = $end3;
	    $MSF_alignments->{'proteins'}->{$protein}->{'alignment'} = [];
	}
    }

    my $replacements;
    my $spaces;
    while($line = <MSF>) {
	if($line =~ /^([\S]+)/) {
	    my $name = $1;
	    if($name =~ /(\/.+)/) {
		$replacements = $1;
		$spaces = " " x length($replacements);
		$line =~ s/$replacements/$spaces/;
		$name =~ s/$replacements//g;
	    }
	    if(exists($MSF_alignments->{'proteins'}->{$name})) {
		push( @{ $MSF_alignments->{'proteins'}->{$name}->{'alignment'} }, $line );
            } else {
		print STDERR "ERROR, $name is not valid protein name\n";
            }
	}
    }

    return $MSF_alignments;

}
	
	
	    
