#!/usr/local/bin/perl

=head1 NAME

create_batch_prokpipe.pl - Creates a prokaryotic annotation pipeline for each row of information

=head1 SYNOPSIS

 USAGE: perl_template.pl
       --input_file=/path/to/some/db_info.tab
       --output_dir=/path/to/output_dir
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>
    Metadata file.  See INPUT section

B<--config_file,-c>
    The pipeline config template file.  Values for parameters are added and saved as a separate config file

B<--log,-l>
    Logfile.

B<--debug,-d>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

 DESCRIPTION

=head1  INPUT
	A tab-separated or comma-separated file containing the following columns:
 1. Repository Root path
 2. Abbreviation (name that will be typically tagged onto output files and headers)
 3. Genome (Two space-separated names)
 4. Taxon ID
 5. Input Fasta File
 6. Input Fasta List (if you want to pass that in instead of a file)
 7. Input BSML List
 8. Accession File Path
 9. Chado DB name
 10. Chado Host server (defaults to manatee-db)
 11. Chado server username
 12. Chado server password
 13. IPD Study Stage (put Study Stage ID)

 For more information about the fields, please consult http://ergatis.igs.umaryland.edu/prokbuilder/help.php#batch

 Each pipeline must also have the exact same template options enabled from Step 1

=head1 OUTPUT
	A list file containing paths to individual comma-separated or tab-delimited files.
	Each file will be in the form of:
	dbname	id_prefix	/path/to/rules.txt	/path/to/gene_syms.txt	Other_db_info	...


=head1  CONTACT
    Shaun Adkins
    sadkins@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use CGI;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;

my $host = "manatee-db";
my $trans_table = 11;

my $input_file;
my $config_file;
my $layout_file;
my %repo_hash = ();
my $stdout;

####################################################

my %options;
my $results = GetOptions (\%options,
                         "input_file|i=s",
                         "config_file|c=s",
                         "log|l=s",
                         "debug|d=s",
                         "help|h"
                          );

&check_options(\%options);

parse_delimited_file($input_file);
foreach my $repo (keys %repo_hash){
    $stdout .= $repo. "=" .$repo_hash{$repo} . "\n";
}

print STDOUT $stdout;
exit(0);

# Parse the input file and write each row to individual files
sub parse_delimited_file {
    my $in = shift;
    my $count = 0;

    open IN, $in or die "Cannot open input file $in for reading: $!\n";
    while (<IN>) {
        my $line = $_;
        chomp $line;
        next if $line =~ /^\s*$/;	# Skip blank lines
        my @row = split (/,|\t/, $line);

        # Use specified value if defined, otherwise use default
        $host = $row[9] if ($row[9]);

        # Using count number, create file path for new config file to be used.
        my $new_config;
        ($new_config = $config_file) =~ s/\.config/_$count\.config/;
        replace_config_values(\@row, $config_file, $new_config);

        # Instantiate our new pipeline
        my $repo_root = $row[0];
        $repo_hash{$repo_root}++;
        my $result = eval { system("/usr/local/bin/perl ./perl/run_prok_pipeline.pl --layout $layout_file --config $new_config --repository_root $repo_root 1>/dev/null") };
        print STDERR $@ unless ($result);	# Catch error if eval fails
        $count++;
    }

    close IN;
    return;
}

# Add parameter values to a config file
sub replace_config_values {
    my ($pipe_data, $config_file, $new_config) = @_;
    open CONF, $config_file or die "Cannot open $config_file for reading: $!\n";
    open NEW, ">".$new_config or die "Cannot open $new_config for writing: $!\n";

    my $line;
    while (<CONF>){
        chomp;
        $line = $_;
        if (/\$;ABBREVIATION\$;=/) {$line .= $pipe_data->[1];}
        if (/\$;GENOME\$;=/) {$line .= $pipe_data->[2];}
        if (/\$;TAXID\$;=/) {$line .= $pipe_data->[3];}
        # Sometimes INPUT FASTA shows up... sometimes INPUT FSA FILE does
        # Also since either the FSA file or list is used, want to suppress any uninitialized value errors from the other
        if (/\$;INPUT_FASTA\$;=/) {$line .= $pipe_data->[4] if (defined $pipe_data->[4]);}
        if (/\$;INPUT_FSA_FILE\$;=/) {$line .= $pipe_data->[4] if (defined $pipe_data->[4]);}
        if (/\$;INPUT_FSA_LIST\$;=/) {$line .= $pipe_data->[5] if (defined $pipe_data->[5]);}
        if (/\$;INPUT_BSML_LIST\$;=/) {$line .= $pipe_data->[6];}
        if (/\$;ACCESSION_FILE\$;=/) {$line .= $pipe_data->[7];}
        if (/\$;DB\$;=/) {$line .= $pipe_data->[8];}
        if (/\$;HOST\$;=/) {$line .= $host;}    # Field 9
        if (/\$;USER\$;=/) {$line .= $pipe_data->[10];}
        if (/\$;PASS\$;=/) {$line .= $pipe_data->[11];}
        if (/\$;IPD_STUDY_STAGE\$;=/) {$line .= $pipe_data->[12];}

        print NEW $line . "\n";
    }

    close CONF;
    close NEW;
    return $new_config;
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

   foreach my $req ( qw(input_file config_file) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }

   $input_file = $opts->{'input_file'};
   die ("Metadata file $input_file does not exist.  Please check filepath.") if (! -e $input_file);
   $config_file = $opts->{'config_file'};
   ($layout_file = $config_file) =~ s/config/layout/;
   die ("Config file $config_file does not exist.") if (! -e $config_file);
   die ("Layout file $layout_file does not exist.") if (! -e $layout_file);
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
