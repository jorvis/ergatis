#!/usr/bin/perl

use lib $ENV{'PERL_MOD_DIR'};

=head1 NAME

sort_snoscan_output.pl - run sort-snos on snoscan output.

=head1 SYNOPSIS

    USAGE: sort_snoscan_output.pl 
                --input_file=/path/to/some_file.snoscan.raw
                --output_file=/path/to/some_file.snoscan.sort.raw
                --generate_out=0 or 1
              [ --position_complement
		--position_hit
		--position_both
		--position_start
		--only_mapped
		--only_unmapped
		--only_top_hits=50
		--require_minimum=score
		--exclude_hits_over=score
		--only_header=expression
		--no_sort
              ]

=head1 OPTIONS

B<--input_file,-i>
    

B<--output_file,-o>
    The file of sorted predicted snos.

B<--generate_out,-g>
    If 1, sort-snos is run and output file is generated.
    If 0, sort-snos will not run.

B<--position_complement,-P>
    Sort snos by Position of complementarity on rRNA.

B<--position_hit, -H>
    Sort snos by position of hit in query sequence.

B<--position_both,-R>
    Sort snos by position & Remove lower-scoring duplicate hits (both start & end bounds must match.

B<--position_start,-r>
    Same as --position_both,-R but only start bound must match to count as duplicate.

B<--only_mapped,-M>
    Sort snos, output only hits to Mapped sites.

B<--only_unmapped,-U>
    Sort snos, output only hits to Unmapped sites.

B<--only_top_hits,-T>
    Sort snos, output only top <int> hits per meth site (def=50).

B<--require_minimum,-S> 
    Sort snos, require minimum score.

B<--only_expression_header,-e> 
    Extract only snos with <expr> in header line.

B<--no_sort,-F>
    Don't sort -- just filter & output in same order.

B<--debug,-d>
    Debug level.  Use a large number to turn on verbose debugging.

B<--log,-l>
    Log file

B<--help,-h>
    This help message.

=head1  DESCRIPTION

This script takes in a snoscan output file and sorts the results according to the options specified above.  It simply runs the sort-snos tool found within the snoscan package.  

=head1  INPUT

The input is defined with --input_file and should be an output file from the program snoscan.  For an example, see 

http://lowelab.ucsc.edu/snoscan/snoscanReadme.html

=head1  OUTPUT

The output is as described in snoscan readme.

http://lowelab.ucsc.edu/snoscan/snoscanReadme.html

If --generate_out,-g is specified as 0(zero), no output is generated and sort-snos is not run.


=head1  CONTACT

    Kevin Galens
    kgalens@tigr.org

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Ergatis::Logger;
use Switch;

#The location of the program
my $SORTPROG = 
    "/usr/local/bin/sort-snos";

my %options = ();
my $results = GetOptions (\%options, 
                          'input_file|i=s',
                          'output_file|o=s',
                          'generate_out|g=i',
                          'position_complement|P',
                          'position_hit|H',
                          'position_both|R',
                          'position_start|r',
                          'only_mapped|M',
                          'only_unmapped|U',
                          'only_top_hits|T=i',
                          'require_minimum|S=i',
                          'only_expression_header|e=s',
                          'no_sort|F',
                          'log|l=s',
                          'debug|d=s',
                          'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});

$logger = $logger->get_logger();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}


## make sure everything passed was peachy
&check_parameters(\%options);

open(OUT, ">$options{output_file}") || 
    $logger->logdie("Could not open outfile");
    #die ("Unable to open output file ($!)");

##If we shouldn't run the sort, just stop.
if( !($options{generate_out}) ) {
    print OUT "No output generated\n";
    close OUT;
    exit(0);
}

##take inputs and build command string
my $cmd_str = "$SORTPROG ";
my $input = $options{input_file};


foreach my $option (keys %options) {
    if(defined($options{$option})) {
	switch($option) {
	    case "position_complement" {$cmd_str.="-P "}
	    case "position_hit" {$cmd_str.="-H "}
	    case "position_both" {$cmd_str.="-R "}
	    case "position_start" {$cmd_str.="-r "}
	    case "only_mapped" {$cmd_str.="-M "}
	    case "only_unmapped" {$cmd_str.="-U "}
	    case "only_top_hits" {$cmd_str.="-T $options{only_top_hits} "}
	    case "require_minimum" {$cmd_str.="-S $options{require_minimum} "}
	    case "only_expression_header" {$cmd_str.="-e $options{only_expression_header}"}
	    case "no_sort" {$cmd_str.="-F "}
	    
	}	
    } 
}

$cmd_str.="$input >$options{output_file}";

exec($cmd_str);

close OUT;

exit;

#########################################################################
#                             SUB ROUTINES                              #
#########################################################################

sub check_parameters {
    my $options = shift;
    
    ## make sure input_file and output_dir were passed
    unless ( $options{input_file} && $options{output_file} 
	     && defined($options{generate_out})) {
        #$logger->logdie("Required options are:\n\t--input_file\n".
	#"\t--output_dir\n\t--fragment_");
        pod2usage({-exitval => 2,  -message => "error message", 
		   -verbose => 1, -output => \*STDERR}); 
	
    }
    
    ## make sure input_file exists
    if (! -e "$options{input_file}") {
	print "Can't find input file\n$options{input_file}\n";
        $logger->logdie("the input file passed ($options{input_file})".
			"cannot be read or does not exist");
    }
        
}
