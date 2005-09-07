#! /usr/local/bin/perl

=head2
	Move_Workflow_InterPro_Results.pl
		
=text

Paolo Amedeo, May 26th 2005

This script is designed for copying InterPro v. 4.0 search results to the 
respective assembly directories.

It takes a list of files and, if they don't already exist, it generates the proper
directories in each assembly directory.

=cut

use strict;
use warnings;
use Getopt::Std;
use File::Copy;
use File::Path;

our ($opt_L, $opt_D, $opt_R, $opt_r, $opt_h);
getopts('L:D:hRr');

MAIN:{
	my ($prN) = ($0 =~ /\/?([^\/]+)$/);
	my $message = "\n\nUsage:   $prN   <Options>\n\n";
	my $safetime = 5; # Time for which the program would pause (allowing killing) in the case option -R has been used
	my $bad = 0;
	my ($proj_dir, $source_dir, $target_dir, $file_end);
	my $summary_tool = 'iprscan';
	my %files = ();


		
	
	if (defined $opt_D){
		$proj_dir = "$ENV{ANNOTATION_DIR}/" . uc($opt_D) . '/asmbls';
	} else {
		$message .= "Option -D (Annotation database) is required\n\n" unless $opt_h;
		++$bad;
	}


	my $old_file_kwd =  qr/($opt_D\.model\.(\d+)(\d{5}))\.(\w+)/ unless $bad;  #avoids complains if $opt_D is not specified
	my $new_file_kwd =  qr/($opt_D\.model\.(\d+)_(\d+))\.(\w+)/ unless $bad;  #avoids complains if $opt_D is not specified
	
	my %prog_to_ev_type = ( 'blastprodom' => "BlastProDom",
				'coils'       => "Coil",
				'hmmpir'      => "HMMpir",
				'hmmsmart'    => 'HMMSmart',
				'hmmtigr'     => 'HMM_TIGR',
				'fprintscan'  => 'FPrintScan',
				'scanregexp'  => 'ScanRegExp',
				'profilescan' => 'ProfileScan',
				'superfamily' => 'superfamily',
				'seg'         => 'Seg');

	if (defined $opt_L && open(my $filelist, $opt_L) &! $bad){
		my ($good, $crap) = (0) x 2;
		while (<$filelist>){
			chomp();
			unless  (/$new_file_kwd/ || $old_file_kwd){
				warn "File: \"$_\" is not recognized by the search pattern\n\n";
				++$crap;
				next;
			}
			my ($mol_name, $asmbl, $model, $tool) = ($1, $2, $3, $4);
			
			next if $tool eq $summary_tool || $tool eq 'hmmpfam' || $tool eq 'hmmtigr';
			unless (exists $prog_to_ev_type{$tool} && defined $prog_to_ev_type{$tool}){
				warn "Tool not found: \"$tool\"\tFile: $_\n";
				next;
			}
			
			push(@{$files{$asmbl}{$prog_to_ev_type{$tool}}}, [$mol_name, $model, $_]);
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

	
	
	

	$message .= "
################################### Options ##################################
#
# -L List of files
#
# -D Annotation database
#
# -R Remove source files
#
# -h Help: prints this option menu and quit
#
###############################################################################

";
	
	die $message if $bad || $opt_h;

	if ($opt_R){
		print "\n\nYou have chosen to remove the source files. Hit Ctrl-C within $safetime secs to stop.";
		sleep($safetime);
	}
	
	while (my ($asmbl, $mdl_files) = each %files){
		my $printed = 0;
		system("$ENV{EGC_SCRIPTS}/ensure_asmbl_dir.dbi -D $opt_D -p $ENV{EGC_SCRIPTS}/egc_password -a $asmbl") && die "\n\nImpossible to find the project directory $proj_dir/$asmbl\n\n" unless -d "$proj_dir/$asmbl";
		
		foreach my $tool (keys %{$mdl_files}){
			my $target_path = "$proj_dir/$asmbl/$tool";
			
			unless (-d $target_path){
				mkpath($target_path) || die "\n\nImpossible to create the directory $target_path\n\n";
				chmod(0777, $target_path); # || warn "Impossible to change permissions to directory $target_path\n";
			}
		
			foreach my $info (@{$mdl_files->{$tool}}){
				my ($mol_name, $model, $file, $target_file) = (@{$info}, $target_path);
				$model = sprintf("%d.m%05d", $asmbl, $model);
				$target_file .= "/$model.$tool";
			
				
				if (-e $target_file){
					chmod(0666, "$target_file"); # || warn "Impossible to change permissions to the pre-existing file $target_file\n";
					unlink("$target_file") || warn "Impossible to delete the pre-existing file $target_file\n";
				}
				
				my $mode = $file =~ /gz$|gzip$/ ? '<:gzip' : '<';
				
				if (open(my $srcfile, $mode, "$file")){
					if (open(my $tgt, ">$target_file")){
						while (<$srcfile>){
							s/$mol_name\S*/$model/g;
							print {$tgt} $_;
						}
						close($tgt);
					} else {
						warn "Impossible to copy the file $file to $target_file\n";
						next;
					}
					close($srcfile);
				} else {
					warn "Impossible to access to the source file $file\n";
					next;
				}
		
				chmod(0666, $target_file); # || warn "Impossible to change permissions to the file $target_file\n";
				
				if ($opt_R ||  $opt_r){
					print "Deleting source file: '$file'..";
					chmod(0666, $file);
					my $outcome = unlink($file) ? "OK" : "Impossible to delete the source file: '$!'\n";
					print $outcome;
				}	
			}
		}
	}
}
