
#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1 NAME

filter_deseq.pl -   filter DESeq output on [FDR, log fold change, Read Count, P value]

=head1 SYNOPSIS

 USAGE: filter_on_fc.pl
	[REQUIRED]
	--deseq_list = list file of .txt output files from DESeq
	--output_dir = output directory to write the UP and DOWN regulated gene lists	
	--filters = string containing combinations of various filters <FDR=0.5,RC=10:FDR=0.3,RC=5:FDR=0.05,RC=15>
	--project_name = this will be the file name for xls document

	[OPTIONAL]
	--map_file = file that has mapping information of [gene_id  TAB  gene_name TAB   gene_desc]
	--abbrev_file = file that has abbreviations for sample names (used for tab names in xls sheet which cannot be > 31 characters)
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
						  'deseq_list=s',
						  'output_dir=s',
						  'filters=s',
						  'project_name=s',
						  'map_file=s',		  # optional
						  
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
my $map_info = {};
my ($file,$IN, $fh, $SA, $out, $summary_file, $summary_all, $i, $p, $stub, $sample_name, $sample1, $sample2, $total_output, $line, $workbook, $cmd, $t);
my (@vals, @arr,  @parameters, @filters, @results_files, @files, @samples, @up, @down);
my %param;
my @su =();
my %outfiles = ();
my $f = 0;




open( $IN, "<", $options{deseq_list} ) or die "Error cannot open the list file";
@results_files = <$IN>;
close $IN;


if( defined $options{ 'map_file' } ) {
	print STDERR "\nGetting map information of gene id and gene name\n";
	$map_info = get_map_info( $options{ 'map_file' } );
}


@filters = split (/\:/,$options{'filters'});



$summary_all = $options{output_dir}."/".$options{project_name}."_summary.txt";
open($SA, ">", $summary_all) or die "Error Cannot open the summary output file";
print $SA "Filters\tSample\tUP\tDOWN\tTotal\n";

foreach $file (@results_files) {
	chomp $file;
	print STDERR "Processing $file...\n";
	# set sample names
	($i, $p, $stub) = File::Spec->splitpath($file);
	($sample1, $i) = split(/\_vs\_/, $stub);
	$sample2 = (split(/\./,$i))[0];
	# output for up and down fc
	(@files, @samples, @up, @down) = ();
	
	for ($i = 1 ;$i <= scalar @filters; $i++) {
	    $total_output = $options{output_dir} . "/" . $sample1 . "_vs_" . $sample2 . "_".$i.".filtered.txt";
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
	$t = 0 ;
	if ($line =~ /ID/) {
	    $t = 1;
	    @arr = split (/\t/,$line);
	
	    foreach (@arr) {
		foreach $i (@files) {     
		    print $i "$_\t";
		}
	    }
	}

	if (exists $options{'map_file'}) {
	    foreach $i (@files) {
		print $i "gene_symbol\tgene_name";
	    }
	}
	if ($t == 1) {
	    foreach $i (@files) {
		print $i "\n";
	    }
	}

	while (<$fh>) {
	    $f = 0;
	    chomp ($_);
	    @vals = split(/\t/,$_);
	    
	    push (@vals,0);
	    push (@vals,0);
	    
	    
		# if user mentioned map_file, we want to print name in output - Mahesh Vangala
	    if( $options{ 'map_file' } ) {
		if( !defined $$map_info{ $vals[ 0 ] }{ 'gene_symbol' } ) {
			$$map_info{ $vals[0] }{ 'gene_symbol' } = "NOT FOUND";
		}
		
		if( !defined $$map_info{ $vals[ 0 ] }{ 'gene_name' } ) {
			$$map_info{ $vals[ 0 ] }{ 'gene_name' } = "NOT FOUND";
		}
		
		push (@vals , $$map_info{ $vals[ 0 ] }{ 'gene_symbol' });
		push (@vals , $$map_info{ $vals[ 0 ] }{ 'gene_name' });
	    }
	    
	    for ($i = 0; $i< scalar @filters ;$i++) {
		
		@parameters = split (/,/,$filters[($i)]);
		%param = ();
		foreach $p (@parameters) {
		    @arr = split(/=/,$p);
		    $param{$arr[0]} = $arr[1];
		}
		if (! defined $param{'UFC'}) {
		    $param{'UFC'} = 1;
		}
		
		if (! defined $param{'DFC'}) {
		    $param{'DFC'} = -1;
		}

		$f = run_filter_checks(\@vals,\%param);
		if ($f == 0) {
		    $out = $files[$i];
		    find_up_or_down_regulated(\@vals, $out, \%param, \$up[$i], \$down[$i]);
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

close $SA ;

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
      $cmd = "grep -v ID ".$i."| sort -k6 -g -r >>". $options{output_dir}."/temp.txt";
      exec_command($cmd);
      $cmd = "mv ".$options{output_dir}."/temp.txt ".$i; 
      exec_command($cmd);
  }
}

for ($i = 1; $i <= scalar @filters; $i++) {
    $workbook = Spreadsheet::WriteExcel->new($options{output_dir} . "/$options{project_name}.RNA_Seq_$i.xls");
    &write_to_excel($workbook,\@{$outfiles{$i}});
}




print STDERR "\n";



###############################

#
# sample1 is up-regulated vs sample2 or down-regulated
#
sub find_up_or_down_regulated {
	my ($vals, $total_out, $param, $up_reg, $down_reg) = @_;
	my $i;
	my $l = 2;
	my $data_out = "";
	
	if (exists $options{'map_file'}) {
	    $l = 4;
	}
	
	foreach ($i = 0; $i< (scalar @{$vals}) - $l ;$i++) {
	    $data_out.= $$vals[$i]."\t";
	}
	
	
	if (exists $options{'map_file'}) {
	    $data_out.= $$map_info{$$vals[0]}{'gene_symbol'}."\t";
	    $data_out.= $$map_info{$$vals[0]}{'gene_name'};
	}
	$data_out.= "\n";
	
	if( $$vals[8] ) {
	    print $total_out $data_out;
	    $$up_reg ++;
	} 

	elsif( $$vals[5] > $param->{'UFC'} ) {
	    print $total_out $data_out;
	    $$up_reg ++;
	}
	if( $$vals[9] ) {
	    print $total_out $data_out;	
	    $$down_reg ++;
	} 
	elsif( $$vals[5] < $param->{'DFC'} ) {
	    print  $total_out $data_out;
	    $$down_reg ++;
	}

}


#
# return = 1 means that the feature will be filtered out
#
sub run_filter_checks {
	my $vals = shift;
	my $param = shift;
	my $filtered = 0;
	my $nMax = 0;

		
        # FDR cutoff not satisfied
	if (exists $param->{'FDR'}) {
	    if( $param->{'FDR'} != 0 && $$vals[7] > $param{'FDR'} ) {
		$filtered = 1;
	    }
	}
	if (exists $param->{'P'}) {
	    if ($$vals[6] > $param->{'P'} ) {
		$filtered = 1;
	    }
	}
       
	# log_fc = Inf; read count for sample2 is 0 and read count for sample1 is greater than the read_count_cutoff
	if( $$vals[2] == 0 && $$vals[3] > $param->{'RC'} ) {
		$$vals[8] = 1;
	}
	
	# log_fc = -Inf; read count for sample1 is 0 and read count for sample2 is greater than the read_count_cutoff
	if( $$vals[3] == 0 && $$vals[2] > $param->{'RC'} ) {                
		$$vals[9] = 1;
	}
	
	$nMax = (($$vals[3] > $$vals[2]) ? $$vals[3] : $$vals[2]);
	if( $nMax < $param->{'RC'} ) {
		$filtered = 1;
	}
		
	return $filtered;
}


sub write_to_excel {
    my ($workbook) = shift ;
    my $files = shift ;
    my ($i,$cmd,$prefix,@sheets, $fh, $len, $j, $t);
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
	    if ($_ =~ /^ID/) {
		$len = chr((64 + scalar(@arr) - 3));
		$worksheet->set_column("A:$len",30);
	    }
	    $t = 0;
	    for ($j = 0; $j< scalar @arr ; $j++) {
		next if ($j == 1 || $j == 4 || $j == 8 || $j == 9);
		if ($t == 4 ) {
		    if ($_ =~ /ID/) {
			$worksheet->write($c,$t,"Abs_log_fold",$format);
		    }
		    else {
			$len = abs($arr[5]);
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
	
	my @required = qw( deseq_list output_dir project_name );
	for my $option ( @required ) {
        unless  ( defined $options{$option} ) {
            die "--$option is a required option";
        }
	if (! defined $options{'filters'} ){
	    $options{'filters'} = "FDR=0.05,RC=10,UFC=1,DFC=-1";
	}
    }

}

sub get_map_info {
	my( $file ) = @_;
	my $info = {};
	
	open( FH, "<$file" ) or die "Error in opening the file, $file, $!\n";
	
	while( my $line = <FH> ) {
		chomp $line;
		next if( $line =~ /^\#/ );
		my( $gene_id, $gene_symbol, $gene_name ) = split( "\t", $line );

		$$info{ $gene_id }{ 'gene_symbol' } = $gene_symbol;
		$$info{ $gene_id }{ 'gene_name' } = $gene_name;
	}
	
	close FH or die "Error in closing the file, $file, $!\n";

	return $info;
}
