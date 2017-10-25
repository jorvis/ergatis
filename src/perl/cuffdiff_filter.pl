#!/usr/bin/env perl

eval 'exec /usr/bin/env perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1 NAME

cuffdiff_filter.pl -  filter cuffdiff output on 
                      [FDR, log_fold_change, FPKM, P value]

=head1 SYNOPSIS

 USAGE: filter_on_fc.pl
	[REQUIRED]
        --cuff_list    = list file of .txt output files from Cuffdiff
	--output_dir   = output directory to write the UP and DOWN regulated gene lists	
	--filters       = string containing combinations of various filters <FDR=0.5,FPKM=10:FDR=0.3,P=0.01:FDR=0.05,FPKM=15>
	--project_name = this will be the file name for xls document

	[OPTIONAL]
        --help

=head1  CONTACT

    Priti Kumari
    pkumari@som.umaryland.edu

=cut

use strict;
use warnings;
use File::Spec;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Spreadsheet::WriteExcel;

my %options = ();
my $results = GetOptions (\%options,
						  'cuff_list=s',
						  'output_dir=s',
						  'filters=s',
						  'project_name=s',
			                          'help|h') || pod2usage();

# display documentation
if( $options{'help'} ) {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

# make sure parameters are correct
&check_parameters();

if( !-d $options{output_dir} ) {
	print STDERR "Output directory doesn't exist.  Creating...\n";
	mkdir $options{output_dir};
}
$options{output_dir} = File::Spec->canonpath($options{output_dir});

my %worksheets;
my ($file,$IN, $fh, $SA, $out, $summary_all, $i, $p, $stub, $sample_name, $sample1, $sample2, $total_output, $line, $workbook, $cmd);
my (@vals, @arr,  @parameters, @filters, @results_files, @files, @samples, @up, @down);
my %param;
my %outfiles = ();
my $f = 0;


open( $IN, "<", $options{cuff_list} ) or die "Error cannot open the list file";
@results_files = <$IN>;
close $IN;


@filters = split (/\:/,$options{'filters'});

$summary_all = $options{output_dir}."/".$options{project_name}."_summary.txt";
open($SA, ">", $summary_all) or die "Error Cannot open the summary output file";
print $SA "Filters\tSample\tUP\tDOWN\tTotal\n";



foreach $file (@results_files) {
	chomp $file;
	print STDERR "Processing $file...\n";
	# set sample names
	($i, $p, $stub) = File::Spec->splitpath($file);
	print "Stub $stub\n";
	($sample1, $i) = split(/vs/, $stub);
	$sample2 = (split(/\./,$i))[0];
	# output for up and down fc
	(@files, @samples, @up, @down) = ();
	
	for ($i = 1 ;$i <= scalar @filters; $i++) {
	    $total_output = $options{output_dir} . "/" . $sample1 . "_vs_" . $sample2 . "_".$i.".txt";
	    if (! exists $outfiles{$i}) {$outfiles{$i} = [$total_output];}
	    else {push (@{$outfiles{$i}}, $total_output);}
	    my $FH ;
	    open($FH, ">", $total_output) or die "Error Cannot open the output file";
	    push (@samples , $sample1 . "_vs_" . $sample2 . "_".$i);
	    push (@files , $FH);
	    push (@up ,0);
	    push (@down ,0);
	}

	open($fh, "<", $file );

	$line = <$fh>;
	chomp($line);
	if ($line =~ /^test_id/) {
	    @arr = split (/\t/,$line);
	}
	foreach (@arr) {
	    foreach $i (@files) {     
		print $i "$_\t";
	    }
	}
	foreach $i (@files) {
	    print $i "\n";
	}

	while (<$fh>) {
	    $f = 0;
	    chomp ($_);
	    @vals = split(/\t/,$_);

	    for ($i = 0; $i< scalar @filters ;$i++) {
		@parameters = split (/,/,$filters[($i)]);
		%param = ();
		foreach $p (@parameters) {
		    @arr = split(/=/,$p);
		    $param{$arr[0]} = $arr[1];
		}
		if (! exists $param{'UFC'}) {
		    $param{'UFC'} = 0;
		}
		if (! exists $param{'DFC'}) {
		    $param{'DFC'} = 0;
		}

		$f = run_filter_checks(\@vals,\%param);
		if ($f == 0) {
		    $out = $files[$i];
		    compare_expression(\@vals, $out, \%param, \$up[$i], \$down[$i]);
		}
		
	    }
	    
	    # report UP regulated genes in sample1 compared to sample2

	}
	for( $i = 0; $i < scalar @filters ;$i++) {
	    close $files[$i];
	    $p = $up[$i] + $down[$i];
	    next if ($p == 0);
	    print "$samples[$i]\tUpregulated: $up[$i]\tDownregulated: $down[$i]\tTotal: $p\n";	
	    print $SA "$filters[$i]\t$samples[$i]\t$up[$i]\t$down[$i]\t$p\n";
	    
	}
	close $fh;

}

close $SA;

foreach (keys %outfiles) {
   foreach $i (@{$outfiles{$_}}) {
      $p = 0;
      open ($fh , "<$i") or die "Cannot open file";
      while(<$fh>) {
	  $p++;
	  last if ($p>3)
      }
      close $fh;
      if ($p == 1) {
	  $cmd = "rm ". $i;
	  exec_command($cmd);
      }
      next if (! -e $i);
      $cmd = "head -1 ".$i." >".$options{output_dir}."/temp.txt";
      exec_command($cmd);
      $cmd = "grep -v test_id ".$i."| sort -k10 -g -r >>". $options{output_dir}."/temp.txt";
      exec_command($cmd);
      $cmd = "mv ".$options{output_dir}."/temp.txt ".$i; 
      exec_command($cmd);
  }
}

for ($i = 1; $i <= scalar @filters; $i++) {
    $workbook = Spreadsheet::WriteExcel->new($options{output_dir} . "/$options{project_name}.cuffdiff_$i.xls");
    &write_to_excel($workbook,\@{$outfiles{$i}});
}




print STDERR "\n";



###############################

#
# sample1 is up-regulated vs sample2 or down-regulated
#
sub compare_expression {
	my ($vals, $total_out, $param, $up_reg, $down_reg) = @_;
	my $i;
	my $data_out = "";
	
	foreach ($i = 0; $i< (scalar @{$vals}) ;$i++) {
		$data_out.= $$vals[$i]."\t";
	}
	$data_out.= "\n";
	print $total_out $data_out;
	
	if( $$vals[9] > $param->{'UFC'} ) {
	    $$up_reg ++;
	}
	elsif( $$vals[9] < $param->{'DFC'} ) {
	    $$down_reg ++;
	}

}


#
# return = 1 means that the feature will be filtered out
#
sub run_filter_checks {
	my $vals = shift;
	my $param = shift;
	my $significance_test;
	
	
        # FDR cutoff not satisfied
	if (exists $param->{'FDR'}) {
	    if( $$vals[12] < $param{'FDR'} ) {
		$significance_test= 0;
	    }
	    else {
		$significance_test = 1;
	    }
	}
	elsif (exists $param->{'P'}) {
	    if( $$vals[11] < $param{'P'} ) {
		$significance_test= 0;
	    }
	    else {
		$significance_test = 1;
	    }
	}
	else {
	    if ( $$vals[13] ne "yes" ) {
		$significance_test = 1;
	    }
	    else {
		$significance_test = 0;
	    }
	}
	if ($$vals[6] ne "OK" || $significance_test) {
	    return 1;
	}

	if (exists $param->{'P'}) {
	    if ($$vals[11] > $param->{'P'} ) {
		return  1;
	    }
	}
	if (exists $param->{'FPKM'}) {
	    if ($$vals[7] < $param->{'FPKM'} && $$vals[8] < $param->{'FPKM'} ) {
		return 1;
	    }
	}	
	
        
	if( $$vals[9] < $param->{'UFC'}  && $$vals[9] > $param->{'DFC'} ) {
	    return 1;
	}
	
       
	
	return 0;
}


sub write_to_excel {
    my ($workbook) = shift ;
    my $files = shift ;
    my ($i,$cmd,$prefix,@sheets, $fh, $j, $len, $t);
    my ($worksheet, $format, $format1);
    my $c = 0;
    foreach(@{$files}) {
	($i,$i,$prefix) = File::Spec->splitpath($_);
	$prefix = (split (/\./,$prefix))[0];
	if (length $prefix > 31) {
	    $prefix = substr($prefix,0,25);
	    $prefix.="_".$c;
	    $c++;
	    
	}
	push(@sheets,$prefix);
    }
    $format = $workbook->add_format(bold=>1,size=>14);
    $format->set_num_format('00.E+00');
    for($i = 0 ; $i< scalar @{$files}; $i++) {
	$c = 0;
	next if (! -e $$files[$i] );
        $worksheet = $workbook->add_worksheet($sheets[$i]);
	open($fh, "<$$files[$i]") or die "Cannot open the output txt file";
	while (<$fh>) {
	    chomp($_);
	    @arr = split (/\t/, $_);
	    if ($_ =~ /^test_id/) {
		$worksheet->set_column("A:J",30);
	    }
	    $t = 0;
	    for ($j = 0; $j< scalar @arr ; $j++) {
		next if ($j == 3 || $j == 4 || $j == 5 || $j == 6 || $j == 10 || $j == 13);
		if ($t == 6 ) {
		    if ($_ =~ /^test_id/) {
			$worksheet->write($c,$t,"Abs_log_fold",$format);
		    }
		    else {
			$len = abs($arr[9]);
			$worksheet->write($c,$t,$len,$format);
		    }
		    $t++;
		}
		$worksheet->write($c,$t,$arr[$j],$format);
		$t++;
	    }
	    $c++;
	}
	close $fh;
    }
    
}
	    
    
sub exec_command {
	my $sCmd = shift;
	
	if ((!(defined $sCmd)) || ($sCmd eq "")) {
		die "\nSubroutine::exec_command : ERROR! Incorrect command!\n";
	}
	
	my $nExitCode;
	
	print STDERR "\n$sCmd\n";
	$nExitCode = system("$sCmd");
	if ($nExitCode != 0) {
		die "\tERROR! Command Failed!\n\t$!\n";
	}
	print STDERR "\n";
	
	return;
}    


sub check_parameters {
	
	my @required = qw( cuff_list output_dir project_name );
	for my $option ( @required ) {
        unless  ( defined $options{$option} ) {
            die "--$option is a required option";
        }
	if (! defined $options{'filters'} ){
	    $options{'filters'} = "FDR=0.05,FPKM=10,UFC=1,DFC=-1";
	}
    }

}


