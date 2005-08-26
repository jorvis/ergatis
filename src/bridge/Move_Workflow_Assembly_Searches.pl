#! /usr/local/bin/perl

=head2
	Move_Workflow_Assembly_Searches.pl
		
=text

Paolo Amedeo, March 28th 2005

This script is used to copy over the results from worflow searches on scaffolds
and other long molecules that need to be recombined.
It has been necessary to separate this script from Move_Workflow_Results.pl due to
practical reasons.

=cut

use strict;
use warnings;
use Getopt::Std;
use File::Copy;
use File::Path;

our ($opt_L, $opt_D, $opt_A, $opt_B, $opt_a, $opt_h, $opt_l, $opt_t, $opt_T, $opt_z);
getopts('L:D:ABahfl:t:T:z');

MAIN:{
	my ($prN) = ($0 =~ /\/?([^\/]+)$/);
	my $fragmented = 1;
	my $message = "\n\nUsage:   $prN   <Options>\n\n";
	my $bad = 0;
	my ($proj_dir, $source_dir, $target_dir, $file_end, $btab, $db);
	my $frg_ln = 50000; # default fragment length
	my $copied = 0;
	my %files = ();
	
	if (defined $opt_A && defined $opt_B || defined $opt_A && defined $opt_a || defined $opt_B && defined $opt_a){
		$message .= "Options -A, -B and -a are incompatible among them: use one and oly one of them\n\n" unless $opt_h;
		++$bad;
	}
	elsif (!defined $opt_A &! defined $opt_B &! defined $opt_a){
		$message .= "you must specify one and only one among options -A, -B and -a\n\n" unless $opt_h;
		++$bad;
	}
	
	if (defined $opt_l && $opt_l =~ /\d/ && $opt_l !~ /\D/){
		if ($opt_l){
			$frg_ln = $opt_l;
		} else {
			$fragmented = 0;
		}
	}
	elsif (defined $opt_l){
		$message .= "Bad value for option -l (Length of the fragment)\n\n" unless $opt_h;
		++$bad;
	}
	
	if (defined $opt_D){
		$proj_dir = "$ENV{ANNOTATION_DIR}/" . uc($opt_D) . '/asmbls';
		$db = $opt_D;
	} else {
		$message .= "Option -D (Source directory) is required\n\n" unless $opt_h;
		++$bad;
	}

	if ($opt_B){
		$btab =  1;
	} else {
		$btab = 0;
	} 
	
	if ($opt_a){
		$frg_ln = 0;
		$fragmented = 0;
	}

	my $file_kwd = $fragmented ? qr/($opt_D\.assembly\.(\d+))\.(\d+)\./ : qr/($opt_D\.assembly\.(\d+))\./ if defined $opt_D;

	if (defined $opt_L && open(my $filelist, $opt_L) &! $bad){
		my ($good, $crap) = (0) x 2;
		while (<$filelist>){
			chomp();
			unless  (/$file_kwd/){
				warn "File: \"$_\" is not recognized by the search pattern\n\n";
				++$crap;
				next;
			}
			my ($mol_name, $asmbl, $segment) = $fragmented ? ($1, $2, $3) : ($1, $2, 0);
			push(@{$files{$asmbl}[$segment]}, [$mol_name, $_]);
			++$good;
		}
		close($filelist);
		die "\n\nThe program has aborted because it was able to recognize $good files but not $crap other ones\n\n" if $crap;
	}
	elsif (defined $opt_L){
		$message .= "Bad value for -L or impossible to open the file $opt_L\n\n";
		++$bad;
	} else {
		$message .= "Option -L (List of files) is required\n\n" unless $opt_h;
		++$bad;
	}

	if (defined $opt_t){
		($file_end = $opt_t) =~ s/^\.//;
		$file_end =~ s/\.gz$//;   # eliminating the possible double gz extension
		$file_end .= '.btab' unless $file_end =~ /\.btab$/;  # Appending btab, in the case it has been forgotten...
		
	} else {
		$message .= "Option -t (Target file 'ending') is required\n\n" unless $opt_h;
		++$bad;
	}
	
	if (defined $opt_T){
		($target_dir = $opt_T) =~ s/^[\/.]+//;
		$target_dir =~ s/\/$//;
		
	} else {
		$message .= "Option -T not specified.\nAssuming as target directory \$ANNOTATION_DIR/DB/asmbls/\$asmbl_id\n\n" unless $opt_h;
	}
	
	
	

	$message .= "
################################### Options ##################################
#
# -L List of files
#
# -D Annotation database
#
# -B btab file
#
# -A Alignment file
#
# -a assembly file
#
# -h print this option menu and quit
#
# -l Length of the fragment (default: $frg_ln nt)
#
# -t Target file 'ending' (i.e. nr.btab, everything after asmbl_id)
#
# -T Targret directory inside $ENV{ANNOTATION_DIR}/DB/asmbls/\$asmbl_id (i.e. blastp)
#
# -z Compress the target file
#
###############################################################################

Note: The program assumes by default that the sequences have been fragmented.
      If the sequence is instead in a single piece, use option -l giving 0 as argument
";
	
	die $message if $bad || $opt_h;

		
	while (my ($asmbl, $asm_files) = each %files){
		my $target_path = "$proj_dir/$asmbl";
		
		system("$ENV{EGC_SCRIPTS}/ensure_asmbl_dir.dbi -D $opt_D -p $ENV{EGC_SCRIPTS}/egc_password -a $asmbl") && die "\n\nImpossible to find the project directory $target_path\n\n" unless -d $target_path;
		
		if (defined $target_dir && $target_dir =~ /\S/){
			$target_path .= "/$target_dir";
			
			unless (-d $target_path){
				mkpath($target_path) || die "\n\nImpossible to create the directory $target_path\n\n";
				chmod(0777, $target_path) || warn "Impossible to change permissions to directory $target_path\n";
			}
		}
		
		foreach my $n (0..$#{$asm_files}){
			unless (defined $asm_files->[$n]){
				warn "No results for Assembly $asmbl - Segment $n\n";
				next;
			}
			foreach my $info (@{$asm_files->[$n]}){
				my ($mol_name, $file, $target_file) = (@{$info}, $target_path);
			
				$target_file .= $btab  ? "/$asmbl.$file_end" : 
				                $fragmented ? "/$asmbl\_s" . $frg_ln * $n . ".$file_end" : "/$asmbl.$file_end";

				if (-e $target_file){
					chmod(0666, "$target_file") || warn "Impossible to change permissions to the pre-existing file $target_file\n";
					unlink("$target_file") || warn "Impossible to delete the pre-existing file $target_file\n" unless $btab && $fragmented;
				}
				my $mode = $file =~ /gz$|gzip$/ ? '<:gzip' : '<';
				
				if (open(my $srcfile, $mode, "$file")){
					my $open_string = $btab && $fragmented ? ">>$target_file" : ">$target_file";
			
					if (open(my $tgt, $open_string)){
						while (<$srcfile>){
							if ($opt_a){
								s/$mol_name\S*/$asmbl $db/g;
							} else {
								s/$mol_name\S*/$asmbl/g;
							}
							print {$tgt} $_;
					
						}
						close($tgt);
					} else {
						warn "Impossible to copy the file $file to $target_file (\"$open_string\")\n";
						next;
					}
					close($srcfile);
				} else {
					warn "Impossible to access to the source file $source_dir/$file\n";
					next;
				}
		
				if ($opt_z){ # requirested to compress the file...
					if (-e "$target_file.gz"){
						chmod(0666, "$target_file.gz") || warn "Impossible to change permissions to the pre-existing file $target_file.gz\n";
						unlink("$target_file.gz") || warn "Impossible to delete the pre-existing file $target_file.gz\n";
					}
					system("gzip $target_file") && warn "Errors compressing the file $target_file\n";
					chmod(0666, "$target_file.gz") || warn "Impossible to change permissions to the file $target_file.gz\n";
				} else {
					chmod(0666, $target_file) || warn "Impossible to change permissions to the file $target_file\n";
				}
			}		
		}
	}
}
