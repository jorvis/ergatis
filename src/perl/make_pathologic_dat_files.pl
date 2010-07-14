#!/usr/local/bin/perl

=head1 NAME

make_pathologic_dat_files.pl - Will create files needed to use pathologic (pathway tools)

=head1 SYNOPSIS

 USAGE: make_pathologic_dat_files.pl
       --genbank_file_list=/path/to/file.gbk.list
       --fasta_file_list=/path/to/assembly.fsa.list
       --ncbi_taxon_id=122587
       --parent_taxon_id=487
       --organism_name="Neisseria meningitidis"
       --project=gnmm04
       --output_directory=/path/to/dir
    [  --codon_table=11
       --strain=M04
       --log=/path/to/file.log
       --debug=4
       --help
     ]

=head1 OPTIONS

B<--genbank_file_list,-g>
    List of genbank annotaiton files

B<--fasta_file_list,-f>
    List of fasta files.  Basename of files should match those in the genbank list

B<--ncbi_taxon_id,-t>
    The taxon id for the organism.  If it's not in NCBI's taxonomy tree (www.ncbi.nlm.nih.gov/Taxonomy)
    then use the --parent_taxon_id option. One of these two should be provided

B<--parent_taxon_id,-a>
    The taxon id of the closest parent to where this organism would fall
    in the NCBI taxonomy.

B<--organism_name,-n>
    The organism name

B<--strain,-s>
    Optional.  Strain name

B<--project,-p>
    Used for the dbname and how the short title for the database in pathway tools will look

B<--output_directory,-o>
    Output directory to put output file in

B<--codon_table,-c>
    Optional, default = 11

B<--log,-l>
    Logfile.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

 DESCRIPTION
 
=head1  INPUT
    Describe the input

=head1 OUTPUT
    Describe the output

=head1  CONTACT

    Kevin Galens
    kgalens@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Copy;
use Data::Dumper;

my $data_files = {};
my ($p_tax, $tax);
my $options = &check_options;

print &make_organism_dat_file($options)."\n";
print &make_genetic_elements_file( $data_files, $options )."\n";

sub make_genetic_elements_file {
    my ($files, $options) = @_;
    my $gedat_file = $options->{'output_directory'}."/genetic-elements.dat";
    my $chrm_counter = 1;

    open(OUT, "> $gedat_file") or die("Could not open $gedat_file for writing: $!");
    
    foreach my $base ( keys %{$files} ) {
        print OUT "ID\tCHROM-".$chrm_counter++."\n";
        print OUT "NAME\t$base\n";
        print OUT "TYPE\t:CHRSM\n";
        print OUT "CIRCULAR?\tNO\n";
        print OUT "ANNOT-FILE\t".$files->{$base}->{'gbk'}->{'basename'}."\n";
        print OUT "SEQ-FILE\t".$files->{$base}->{'fsa'}->{'basename'}."\n";
        print OUT "//\n";

        #copy the annot and seq files to the directory
        copy( $files->{$base}->{'gbk'}->{'file'}, $options->{'output_directory'}."/".$files->{$base}->{'gbk'}->{'basename'} );
        copy( $files->{$base}->{'fsa'}->{'file'}, $options->{'output_directory'}."/".$files->{$base}->{'fsa'}->{'basename'} );
    }

    close(OUT);
    return $gedat_file;
}

sub make_organism_dat_file {
    my ($options) = @_;
    my $org_dat_file = $options->{'output_directory'}."/organism-params.dat";
    
    open(OUT, "> $org_dat_file") or die("Can't open $org_dat_file for writing: $!");
    my $id = uc($options->{'project'});
    my $name = $options->{'organism_name'};
    my $dbname = ucfirst( lc($options->{'project'}) );
    
    print OUT "ID\t$id\nSTORAGE\tFILE\nNAME\t$name\nDBNAME\t$dbname\n";
    if( $options->{'ncbi_taxon_id'} ) {
        print OUT "NCBI-TAXON-ID\t".$options->{'ncbi_taxon_id'}."\n";
    } elsif( $options->{'parent_taxon_id'} ) {
        print OUT "DOMAIN\tTAX-".$options->{'parent_taxon_id'}."\n";
    }

    if( $options->{'strain'} ) {
        print OUT "STRAIN\t".$options->{'strain'}."\n";
        print OUT "RANK\t|strain|\n";
    }

    if( $options->{'codon_table'} ) {
        print OUT "CODON_TABLE\t".$options->{'codon_table'}."\n";
    }
    
    close(OUT);
    return $org_dat_file;
}

sub check_options {
    my %options;
    my $results = GetOptions (\%options,
                              'genbank_file_list|g=s',
                              'fasta_file_list|f=s',
                              'ncbi_taxon_id|t=s',
                              'parent_taxon_id|a=s',
                              'organism_name|n=s',
                              'strain|s=s',
                              'project|p=s',
                              'output_directory|o=s',
                              'codon_table|c=s',
                              'log|l=s',
                              'help'
                              );
    
    if( $options{'help'} ) {
        &_pod;
    }

    foreach my $req ( qw(genbank_file_list fasta_file_list organism_name project output_directory ) ) {
        die("Option --$req is required") unless( $options{$req} );
    }

    $data_files = &parse_input_lists( $options{'genbank_file_list'}, $options{'fasta_file_list'} );

    unless( $options{'ncbi_taxon_id'} || $options{'parent_taxon_id'} ) {
        die("Either the --ncbi_taxon_id or the --parent_taxon_id options is required") ;
    }
       
    return \%options;
}

sub parse_input_lists {
    my ($gbks, $fsas) = @_;
    my %retval;

    foreach my $input ( [$gbks, "gbk"], [$fsas, "fsa"] ) {
        my ($list,$type) = @{$input};
        open(IN, "< $list") or die("Could not open $list: $!");
        chomp( my @files = <IN> );
        close(IN);
        
        foreach my $file ( @files ) {
            my $base = $2 if( $file =~ /.*\/((.*)\..*)/ );
            $retval{$base}->{$type} = { 'file' => $file,
                                        'basename' => $1 };
        }
    }
    
    #check to make sure each entry in the hash has a gbk and fsa entry
    map { 
        die("genbank_file_list and fasta_file_list should have an equal number of ".
            "entries and the base names should match. They do not.") 
            unless( exists($retval{$_}->{'fsa'} ) && exists( $retval{$_}->{'gbk'} ) );
    } keys %retval;

    return \%retval;
}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
