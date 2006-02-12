#! /local/perl/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

my $str =  "cat $ARGV[0]".' | '."$ARGV[1]".'>'."$ARGV[2]";
system $str;
