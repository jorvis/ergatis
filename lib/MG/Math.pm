package MG::Math;

use vars qw(@ISA @EXPORT $VERSION);
use Exporter;

@ISA = qw(Exporter);
@EXPORT = qw(&min &max &telltime &currdate &mean_var &corr &create_matrix_jaccard &create_matrix_lnjaccard &create_matrix_spearman &create_matrix_euclidean &create_matrix_rmsd &create_matrix_nrmsd &create_matrix_abssum);


sub min {
    @ary = sort {$a <=> $b} @_;
    return $ary[0];
}
sub max {
    @ary = sort {$b <=> $a} @_;
    return $ary[0];
}
sub telltime {
    @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
    ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    $year = 1900 + $yearOffset;
    $thetime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
    return $thetime;
}
sub currdate {
    ($second, $minute, $hour, $dayOfMonth, $month, $year, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    $theday = join("-",$year+1900,sprintf("%02s",$month+1),sprintf("%02s",$dayOfMonth));
    return $theday;
}
sub create_matrix_jaccard {
    my @files = @{shift @_};
    my %clusters = %{shift @_};
    my %mhash = ();
    my @temp = @files;
 F1:foreach my $f1 (@files) {
	shift @temp;
	@ogrps = keys %{$clusters{$f1}};
    F2:foreach my $f2 (@temp) {
	    #my $a = 0;
	    #my $b = 0;
	    my $total = 0;
	    my $share = 0;
	F3:foreach my $orgrp (@ogrps) {
		next F3 if $orgrp =~ m/single/;
		if ($clusters{$f1}{$orgrp} > 0 || $clusters{$f2}{$orgrp} > 0) {
		    $total ++;
		}if ($clusters{$f1}{$orgrp} > 0 && $clusters{$f2}{$orgrp} > 0) {
		    $share ++;
		}
	    }
	    $mhash{$f1}{$f2} = 1 - ($share/$total);
	    $mhash{$f2}{$f1} = 1 - ($share/$total);
	}
    }
    return \%mhash;
}
sub create_matrix_lnjaccard {
    my @files = @{shift @_};
    my %clusters = %{shift @_};
    my %mhash = ();
    my @temp = @files;
 F1:foreach my $f1 (@files) {
	shift @temp;
	@ogrps = keys %{$clusters{$f1}};
    F2:foreach my $f2 (@temp) {
	    my $total = 0;
	    my $share = 0;
	F3:foreach my $orgrp (@ogrps) {
		next F3 if $orgrp =~ m/single/;
		if ($clusters{$f1}{$orgrp} > 0 || $clusters{$f2}{$orgrp} > 0) {
		    $total ++;
		}if ($clusters{$f1}{$orgrp} > 0 && $clusters{$f2}{$orgrp} > 0) {
		    $share ++;
		}
	    }
	    $mhash{$f1}{$f2} = -log($share/$total);
	    $mhash{$f2}{$f1} = -log($share/$total);
	}
    }
    return \%mhash;
}
sub create_matrix_spearman { #sum of squares normalized
    my @files = @{shift @_};
    my %clusters = %{shift @_};
    my %mhash = ();
    my %rank = ();
    foreach my $ds (@files) {
	my @r = ();
	foreach $og (keys %{$clusters{$ds}}) {
	    push @r, [$og, $clusters{$ds}{$og}];
	}
	@r = sort {$b->[1] <=> $a->[1]} @r;
	foreach $i (0..$#r) {
	    $rank{$ds}{$r[$i]->[0]} = $i+1;
	}
    }
    my @temp = @files;
    foreach my $f1 (@files) {
	shift @temp;
	@ogrps = keys %{$clusters{$f1}};
	foreach my $f2 (@temp) {
	    my $sum = 0;
	    foreach my $orgrp (@ogrps) {
		$sum += abs($rank{$f1}{$orgrp} - $rank{$f2}{$orgrp}) ** 2;
	    }
	    $mhash{$f1}{$f2} = 1 - ((6*$sum)/(scalar(@ogrps)*((scalar(@ogrps)**2)-1)));
	    $mhash{$f2}{$f1} = 1 - ((6*$sum)/(scalar(@ogrps)*((scalar(@ogrps)**2)-1)));
	}
    }
    return \%mhash;
}
sub create_matrix_euclidean { #sum of squares normalized
    my @files = @{shift @_};
    my %clusters = %{shift @_};
    my %mhash = ();
    my @temp = @files;
    foreach my $f1 (@files) {
	shift @temp;
	@ogrps = keys %{$clusters{$f1}};
	foreach my $f2 (@temp) {
	    my $sum = 0;
	    foreach my $orgrp (@ogrps) {
		$sum += ($clusters{$f1}{$orgrp} - $clusters{$f2}{$orgrp}) ** 2;
	    }
	    $mhash{$f1}{$f2} = sqrt($sum);
	    $mhash{$f2}{$f1} = sqrt($sum);
	}
    }
    return \%mhash;
}
sub create_matrix_rmsd { #sum of squares normalized
    my @files = @{shift @_};
    my %clusters = %{shift @_};
    my %mhash = ();
    my @temp = @files;
 F1:foreach my $f1 (@files) {
	shift @temp;
	@ogrps = keys %{$clusters{$f1}};
    F2:foreach my $f2 (@temp) {
	    my $n = 0;
	    my $sum = 0;
	F3:foreach my $orgrp (@ogrps) {
		#next F3 if $orgrp =~ m/single/;
		if ($clusters{$f1}{$orgrp} > 0 || $clusters{$f2}{$orgrp} > 0) {
		    $n ++;
		}
		$sum += ($clusters{$f1}{$orgrp} - $clusters{$f2}{$orgrp}) ** 2;
	    }
	    $mhash{$f1}{$f2} = sqrt($sum/$n);
	    $mhash{$f2}{$f1} = sqrt($sum/$n);
	}
    }
    return \%mhash;
}
sub create_matrix_abssum { #sum of squares normalized
    my @files = @{shift @_};
    my %clusters = %{shift @_};
    my %mhash = ();
    my @temp = @files;
    foreach my $f1 (@files) {
	shift @temp;
	@ogrps = keys %{$clusters{$f1}};
	foreach my $f2 (@temp) {
	    my $sum = 0;
	    foreach my $orgrp (@ogrps) {
		$sum += abs($clusters{$f1}{$orgrp} - $clusters{$f2}{$orgrp});
	    }
	    $mhash{$f1}{$f2} = $sum;
	    $mhash{$f2}{$f1} = $sum;
	}
    }
    return \%mhash;
}
sub create_matrix_nrmsd { #sum of squares normalized
    my @files = @{shift @_};
    my %clusters = %{shift @_};
    my @temp = @files;
    my %mhash = ();
 F1:foreach my $f1 (@files) {
	shift @temp;
	@ogrps = keys %{$clusters{$f1}};
    F2:foreach my $f2 (@temp) {
	    my $n = 0;
	    my @values = ();
	    my $sum = 0;
	F3:foreach my $orgrp (@ogrps) {
		#next F3 if $orgrp =~ m/single/;
		if ($clusters{$f1}{$orgrp} > 0 || $clusters{$f2}{$orgrp} > 0) {
		    $n ++;
		    push @values, ($clusters{$f1}{$orgrp},$clusters{$f2}{$orgrp});
		}
		$sum += ($clusters{$f1}{$orgrp} - $clusters{$f2}{$orgrp}) ** 2;
	    }
	    $mhash{$f1}{$f2} = sqrt($sum/$n)/(max(@values) - min(@values));
	    $mhash{$f2}{$f1} = sqrt($sum/$n)/(max(@values) - min(@values));
	}
    }
    return \%mhash;
}
sub mean_var  {
    my @array = @_;
    $sum = 0;
    $total = $#array + 1;
    foreach $element (@array) {
        $sum += $element;
    }
    $mean = $sum/$total;
    $prevar = 0;
    foreach $element (@array) {
        $prevar += ($element - $mean)**2
    }
    $variance = $prevar/$#array;
    return ($mean, $variance);
}
sub corr {
    my @v1 = @{shift @_};
    my @v2 = @{shift @_};
    if (scalar(@v1) ne scalar(@v2)) {
	return("error");
    }
    my $n = scalar(@v1);
    my $sumx = 0;
    my $sumy = 0;
    my $sumyx = 0;
    my $sumx2 = 0;
    my $sumy2 = 0;
    foreach my $i(0..$#v1) {
	$x = $v1[$i];
	$y =  $v2[$i];
	$sumyx += ($y * $x);
	$sumx += $x;
	$sumy += $y;
	$sumx2 += ($x * $x);
	$sumy2 += ($y * $y);
    }
    $top = ($n * $sumyx) - ($sumx * $sumy);
    $bottom = sqrt((($n * $sumx2)-($sumx * $sumx))*(($n * $sumy2)-($sumy * $sumy)));
    if ($bottom > 0) {
	$corr = sprintf("%.2f", $top/$bottom);
    }else {
	$corr = 1;
    }
    return($corr)
}

return 1;
