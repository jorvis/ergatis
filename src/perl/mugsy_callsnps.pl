#!/usr/bin/env perl

=head1 NAME

mugsy_callsnp.pl - Will call SNPs from whole genome alignment MAF

=head1 SYNOPSIS

 USAGE: mugsy_callsnp.pl
       --input_indexed_maf=/path/to/file.maf
       --output_file=/path/to/output.txt
       --fasta=/path/to/multi.fasta
       --mugsy_mapping_dir=/path/to/mugsy/mapping
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_indexed_maf,-i>
    Path to input indexed maf file, indexed with mafindex.pl 

B<--output_file,-o>
    Path to output tab file

B<--fasta,-f>
    Path to all sequences in one multi fasta file.

B<--mugsy_mapping_dir,-m>
    Path to mugsy mapping install directory.

B<--log,-l>
    Logfile.

B<--debug,-d>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  INPUT

    Reads in an indexed MAF (multiple alignment format) file. Expects sequence identifiers to be in the format of
    organism.molecule. Should be indexed with mafindex.pl.

=head1 OUTPUT

    Output is tab delimited file. Three columns for each query organism. Will report any location where there is
    a difference in base that is not a gap.

=head1  CONTACT

    Kevin Galens
    kgalens@gmail.com

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

## BioPerl 
use Bio::Perl;
use Bio::DB::Fasta;
use Bio::Seq;

## Stolen from here
## http://stackoverflow.com/questions/1672782/fastest-way-to-find-mismatch-positions-between-two-strings-of-the-same-length
use Inline C => << '...';                                                       
  void find_diffs(char* x, char* y) {                                           
    int i;                                                                      
    Inline_Stack_Vars;                                                          
    Inline_Stack_Reset;                                                         
    for(i=0; x[i] && y[i]; ++i) {                                               
      if(x[i] != y[i] && x[i] != '-' && y[i] != '-') {                                                        
        Inline_Stack_Push(sv_2mortal(newSViv(i)));                              
      }                                                                         
    }                                                                           
    Inline_Stack_Done;                                                          
  }
...

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
####################################################

my %options;
my $results = GetOptions (\%options,
			  "input_indexed_maf|i=s",
			  "output_file|o=s",
			  "fasta|f=s",
			  "mugsy_mapping_dir|m=s",
			  "log|l=s",
			  "debug|d=s",
			  "help|h"
                          );

&check_options(\%options);

my $atree = AlignmentTree::deserialize($options{'input_indexed_maf'});
my $fh = Bio::DB::Fasta->newFh($options{'fasta'},'-reindex'=>1);
my $db = {};
$db->{$_->id} = $_ while( <$fh> );

my %molecules;
my $snps = [];

# Grab all the alignments in the alignment tree
foreach my $alnname ( sort {$a cmp $b} keys %{$atree->{_alignments}} ) {
    &_log($DEBUG, "Starting alignment named: $alnname");
    
    # Grab the alignment information
    my ($alnobj,$aln_bv,$align_width) = @{$atree->{_alignments}->{$alnname}}; 
    my ($mmatrix,$seqmatrix,$names) = $atree->getAlignmentMatrix($alnname,1,$align_width,$db);

    # Make sure we've stored all the names
    for( @{$names} ) {
	my ($org, $mol) = &get_org_and_mol( $_ );
        $molecules{$org}->{$mol} = 1 unless( exists( $molecules{$org}->{$mol} ) );
    }
    
    &find_snps( $seqmatrix, $names, $atree, $alnname, $snps );
}

&print_snps( \%molecules, $options{'output_file'}, $snps );
&_log($DEBUG, "Output: $options{'output_file'}");

sub print_snps {
    my ($molecules, $output_file, $snps) = @_;

    my @ordered_names = sort keys %{$molecules};
    open( OUT, "> $output_file") or &_log($ERROR, "Can't open $output_file: $!");
    print OUT "$_:mol\t$_:pos\t$_:base\t" for( @ordered_names );
    print OUT "\n";
    foreach my $variant ( @{$snps} ) {
	my @line;
	foreach my $name ( @ordered_names ) {
	    if( exists( $variant->{$name} ) ) {
		push( @line, @{$variant->{$name}} );
	    } else {
		push( @line, "", "", "" );
	    }
	}
	print OUT join("\t", @line)."\n";
    }
    close(OUT);
}

sub find_snps {
    my ($seqmatrix, $names, $atree, $alnname, $data) = @_;
    
    my @all_diffs;
    my $start = time;
    for( my $i = 0; $i < @{$seqmatrix} - 1; $i++ ) {
        for( my $j = $i + 1; $j < @{$seqmatrix}; $j++ ) { 
           push(@all_diffs, find_diffs( $seqmatrix->[$i], $seqmatrix->[$j] ) );
        }
    }
    
    my %alnis;
    foreach my $name ( @{$names} ) {
        $alnis{$name} = $atree->getAlignedInterval( $alnname, $name );
    }
    
    #Coordinates stored are zero based coordinates in relation to the aligned sequence.
    my %uniq_pos;
    foreach my $d ( @all_diffs ) {
        if( exists( $uniq_pos{$d} ) ) {
            next;
        } else {
            $uniq_pos{$d} = 1;
        }
        
        my $tmp = {};
        for( my $i = 0; $i < @{$names}; $i++ ) {
		    
	    # columntocoords takes 1 based coords, $d is zero based so add 1.
	    # returns zero based coordinates
            my ($start, $end) = AlignmentTree::columntocoords( $alnis{$names->[$i]}, $d+1, $d+1 );
            
            # Should store first by organism, then by molecule
	    my ($org, $mol) = &get_org_and_mol( $names->[$i] );
	    # SNP position returned should be 1 based in the final output so added 1 to $start
            $tmp->{$org} = [$mol, ($start+1), substr( $seqmatrix->[$i], $d, 1 )];
        }
        push(@{$data}, $tmp);
    }
}

sub get_org_and_mol {
    my ($string) = @_;
    my ($o,$m) = ($1, $2) if( $string =~ /^([^\.]+)\.(.*)$/ );
    die("Could not parse organism or name from $string") unless( $o && $m );
    ($o, $m);
}



sub check_options {
    my $opts = shift;
    
    if( $opts->{'help'} ) {
	&_pod;
    }
    
    if( $opts->{'log'} ) {
	open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
    }
    
    $debug = $opts->{'debug'} if( $opts->{'debug'} );

    foreach my $req ( qw(input_indexed_maf fasta output_file mugsy_mapping_dir) ) {
	&_log($ERROR, "Option $req is required") unless( $opts->{$req} );
    }

    ## This should point to installation of mugsy
    ## for AlignmentTree
    eval("use lib('$opts->{mugsy_mapping_dir}')");
    die($@) if( $@ );
    eval("use AlignmentTree");# or die($@);
    die($@) if( $@ );
}

sub _log {
   my ($level, $msg) = @_;
   if( $level <= $debug ) {
      print STDOUT "$msg\n";
   }
   print $logfh "$msg\n" if( defined( $logfh ) );
   exit(1) if( $level == $ERROR );
}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
