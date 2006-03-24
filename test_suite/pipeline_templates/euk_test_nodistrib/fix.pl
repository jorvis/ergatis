#!/usr/bin/perl
$fname = shift @ARGV;

open (IN, $fname);
open (OUT, ">out");
while (<IN>) {
	print OUT $_; 
	if (/^\[workflowdocs/) { 
		print OUT "\$;TAG\$;                 = \$Name$\n\$;NODISTRIB\$;    = 0\n\$;REVISION\$;            = \$Revision$\n";
	}
}
close IN;
close OUT;
system("mv out $fname");
