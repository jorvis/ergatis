#! /local/perl/bin/perl
use strict;
use warnings;

# input: an obo file
# output: the number of non is_obsolete terms
#         eg count([Term]) - count(is_obsolete)
#         <#term>  <# obsolete>  <# non-obsolete> <# relationships>

defined ($ARGV[0]) or print_usage('No input file specified.');
my $ifile = $ARGV[0];
(-r $ifile) or print_usage("Input file $ifile is not a readable file.");

my $n_term = 0;
my $n_obs = 0;
my $n_isa = 0;
open(my $FIN, $ifile) || die "Unable to open $ifile: $!";
while (my $line = <$FIN>) {
    if ($line =~ /\[Term\]/) {
	++$n_term;
    }
    elsif ($line =~ /is_obsolete:\s*true/) {
	++$n_obs;
    }
    elsif ($line =~ /is_a:/) {
	++$n_isa; #regardless of obsolete status
    }
}
close($FIN);

print "Input file: $ifile\n";
print "Total terms\t$n_term\n";
print "  Obsolete\t$n_obs\n";
print "  Non-obsolete\t".($n_term - $n_obs)."\n";
print "is_a relationships\t$n_isa\n";

#
# subroutines
#
sub print_usage
{
        my $byebye = shift;
        my $progname = $0;
        die << "END";
$byebye
usage: $progname <input file>
output: <#terms>  <# obsolete>  <# non-obsolete> <# is_a relationships>
END
}
