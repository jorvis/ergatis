#!/user/local/bin/perl
=head1 NAME

btab2taxon.pl

=head1 SYNOPSIS

 USAGE: btab2taxon.pl
       --input=/path/to/file/list
       --ncbitax=/path/to/taxonomy
       --gi2tax=/path/to/gi2tax
       --dbhost=database.host.com
       --taxondb=taxon_db_name
       --taxoncoll=taxon_collection
      [--help]
     
=cut
use strict;
use GiTaxon;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
$| = 1;

my %options;
my $results = GetOptions (\%options,
                          'input=s', # Comma separated list of files
                          'input_file_list=s',
                          'ncbitax=s',
                          'gitax=s',
                          'dbhost=s',
                          'taxondb=s',
                          'taxoncoll=s',
                          'mapfile=s',
                          'overwrite=s',
                          'idx_dir=s',
                          'load_db_only=s',
                          'help|h'
                          );

my $ncbitax = $options{ncbitax} ? $options{ncbitax} : '/local/db/repository/ncbi/blast/20110831_122720/taxonomy/taxdump/';
my $gi2tax = $options{gitax} ? $options{gitax} : '/local/db/repository/ncbi/blast/20110831_122720/taxonomy/gi_taxid_nucl.dmp';
my $dbhost = $options{dbhost} ? $options{dbhost} : 'tettelin-lx.igs.umaryland.edu';
my $taxondb = $options{taxondb} ? $options{taxondb} : 'gi2taxon';

my $taxoncoll = $options{taxoncoll};
if(!$taxoncoll) {
    $ncbitax =~ /(\d+\_\d+)/;
    my $date = $1;
    if($gi2tax =~ /nuc/) {
        $taxoncoll = "gi2taxonnuc_$date";
    }
    else {
        $taxoncoll = "gi2taxonprot_$date";
    }

}
my $idx_dir = $options{idx_dir};

if(!$idx_dir && -e "$ncbitax/names") {
    $idx_dir = $ncbitax;
}
elsif(!$idx_dir) {
    $idx_dir='/tmp/';
}

my $gi2tax = GiTaxon->new(
    {'nodes' => "$ncbitax/nodes.dmp",
     'names' => "$ncbitax/names.dmp",
     'gi2tax' => $gi2tax,
     'chunk_size' => 100000,
     'idx_dir' => $idx_dir,
     'host' => $dbhost,
     'gi_db' => $taxondb,
     'gi_coll' => $taxoncoll
    });    

if($options{load_db_only}) {
    exit(0);
}

my $files = [];

if($options{input}) {
    my @f = split(/,/,$options{input});
    $files = \@f;
}
elsif($options{input_file_list}) {
    open IN, "<$options{input_file_list}" or die "Unable to open list";
    
    while(<IN>) {
        chomp;
        push(@$files,$_);
    }
}

my $trace_lookup = {};

if($options{mapfile}) {
    &load_map_file();
}


foreach my $file(@$files) {
    &process_file($file);
}

sub process_file {
    my $file = shift;

    my $out_fh = *STDOUT;
    if($options{overwrite}) {
        open $out_fh, ">$file\_tmp";
    }
    open IN, "<$file" or die "Unable to open $file\n";
    while(<IN>) {
        chomp;
        my $line = $_;
        my @new_fields = split(/\t/,$line);
        my @fields = @new_fields;
        if($options{mapfile}) {
            @fields = @{&add_trace_info(\@new_fields)};
        }
            
        my $tax = $gi2tax->getTaxon($fields[1]);
        if($tax->{'taxon_id'}) {
            push(@fields,$tax->{'taxon_id'});
            
            if($tax->{'name'}) {
                push(@fields,($tax->{'name'},$tax->{'lineage'}));
            }
        }
        print $out_fh join("\t",@fields);
        print $out_fh "\n";
    }
    if($options{overwrite}) {
        print `mv $file\_tmp $file`;
    }
}
sub add_trace_info {
    my $list = shift;

    my $id = $list->[0];
    if($trace_lookup->{$list->[0]}) {
        push(@$list,($trace_lookup->{$id}->{'template_id'},$trace_lookup->{$id}->{'trace_end'}));
        return $list;
    }
    else {
        print STDERR "couldn't find trace info for $list->[0]\n";
        return [];
    }
}

sub load_map_file {
    open MAP, "<$options{map_file}" or die "Unable to open $options{map_file}\n";
    while(<MAP>) {
        my @fields = split;
        $trace_lookup->{$fields[0]} = {
            'template_id' => $fields[2],
            'trace_end' => $fields[1]
        };
    }
    close MAP;
}
