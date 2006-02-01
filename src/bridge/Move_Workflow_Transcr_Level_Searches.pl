#! /usr/local/bin/perl

=head2
	Move_Workflow_Transcr_Level_Searches.pl
		
=text

Paolo Amedeo, March 29th 2005

This script is used to copy over the results from worflow transcript-level
searches to the proper project directories, with the proper name format.

It is a modification of the original Move_WorkFlow_Results.pl: there are 
too many differences between genomic and transcript level computes that
it was better to split the original script into this and
Move_Workflow_Assembly_Searches.pl

=cut

use strict;
use warnings;
use Getopt::Std;
use File::Copy;
use File::Path;

our ($opt_L, $opt_D, $opt_h, $opt_p, $opt_t, $opt_T, $opt_z, $opt_R, $opt_r);
getopts('L:D:hp:t:T:zRr');

MAIN:{
	my ($prN) = ($0 =~ /\/?([^\/]+)$/);
	my $safetime = 5; # Time for which the program would pause (allowing killing) in the case option -R has been used
	my $message = "\n\nUsage:   $prN   <Options>\n\n";
	my $bad = 0;
	my ($proj_dir, $source_dir, $target_dir, $file_end, $file_pfix);
	my %files = ();
#	my ($proj_dir, $source_dir, $target_dir, $file_end, $btab);

	if (defined $opt_D){
		$proj_dir = "$ENV{ANNOTATION_DIR}/" . uc($opt_D) . '/asmbls';
	} else {
		$message .= "Option -D (Source directory) is required\n\n" unless $opt_h;
		++$bad;
	}


# 	my $old_file_kwd =  qr/($opt_D\.model\.(\d+)(\d{5})(?:\.\d[\d.]+)*)\./ unless $bad;  #avoids complains if $opt_D is not specified
	my $new_file_kwd =  qr/($opt_D\.model\.(\d+)_(\d+)(?:\.\d[\d.]+)*)\./ unless $bad;  #avoids complains if $opt_D is not specified

	if (defined $opt_L && open(my $filelist, $opt_L) &! $bad){
		my ($good, $crap) = (0) x 2;
		while (<$filelist>){
			chomp();
			unless  (/$new_file_kwd/){       # || /$old_file_kwd/){
				warn "File: \"$_\" is not recognized by the search pattern\n\n";
				++$crap;
				next;
			}
			my ($mol_name, $asmbl, $model) = ($1, $2, $3);
			push(@{$files{$asmbl}}, [$mol_name, $model, $_]);
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
		$file_end = $opt_t;
                $file_end = ".$file_end" unless $file_end =~ /^\./;
		$file_end =~ s/gz$//;
        } else {
                $file_end = '';
        }
	
        
        if (defined $opt_p){
		$file_pfix = $opt_p;
                $file_pfix .= '.'  unless $file_pfix =~ /\.$/;
        } else {
                $file_pfix = '';
        }
        
        unless (defined $opt_p || defined $opt_t){
		$message .= "One among options -t and -p (Target file 'ending' or prefix) is required\n\n" unless $opt_h;
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
# -h Help: prints this option menu and quit
#
# -p Prefix (i.e. something to put before the model name in the name of the target file)
#
# -t Target file 'ending' (i.e. nr.btab, everything after asmbl_id or model name)
#
# -T Targret directory inside $ENV{ANNOTATION_DIR}/DB/asmbls/\$asmbl_id (i.e. blastp)
#
# -z Compress the target file
#
# -R Remove source files
#
###############################################################################

";
	
	die $message if $bad || $opt_h;
	
	if ($opt_R){
		print "\n\nYou have chosen to remove the source files. Hit Ctrl-C within $safetime secs to stop.";
		sleep($safetime);
	}

	while (my ($asmbl, $mdl_files) = each %files){
		my $target_path = "$proj_dir/$asmbl";
		
		system("$ENV{EGC_SCRIPTS}/ensure_asmbl_dir.dbi -D $opt_D -p $ENV{EGC_SCRIPTS}/egc_password -a $asmbl") && die "\n\nImpossible to find the project directory $target_path\n\n" unless -d $target_path;
		
		if (defined $target_dir && $target_dir =~ /\S/){
			$target_path .= "/$target_dir";
			
			unless (-d $target_path){
				mkpath($target_path) || die "\n\nImpossible to create the directory $target_path\n\n";
				chmod(0777, $target_path); # || warn "Impossible to change permissions to directory $target_path\n";
			}
		}
	
		foreach my $info (@{$mdl_files}){
			my ($mol_name, $model, $file, $target_file) = (@{$info}, $target_path);

			$model = "$asmbl.m$model";

			#$model = $opt_N ? sprintf("%d.m%06d", $asmbl, $model) : sprintf("%d.m%05d", $asmbl, $model);
			
                        $target_file .= "/$file_pfix$model$file_end";
			
			print STDERR "Moving $file to $target_file.. ";
			
			if (-e $target_file){
				chmod(0666, "$target_file"); # || warn "Impossible to change permissions to the pre-existing file $target_file\n";
				unlink("$target_file") || warn "Impossible to delete the pre-existing file $target_file\n";
			}
			
			my $mode = $file =~ /gz$|gzip$/ ? '<:gzip' : '<';
			
			if (open(my $srcfile, $mode, "$file")){
				if (open(my $tgt, ">$target_file")){
					while (<$srcfile>){
						s/$mol_name/$model/g;
						print $tgt $_;
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
		
			if ($opt_z){ # requirested to compress the file...
				if (-e "$target_file.gz"){
					chmod(0666, "$target_file.gz"); # || warn "Impossible to change permissions to the pre-existing file $target_file.gz\n";
					unlink("$target_file.gz") || warn "Impossible to delete the pre-existing file $target_file.gz\n";
				}
				system("gzip $target_file") && warn "Errors compressing the file $target_file\n";
				chmod(0666, "$target_file.gz"); # || warn "Impossible to change permissions to the file $target_file.gz\n";
			} else {
				chmod(0666, $target_file); # || warn "Impossible to change permissions to the file $target_file\n";
			}
			print STDERR "OK\n";
			
			if ($opt_R ||  $opt_r){
				print "Deleting source file: '$file'..";
				chmod(0666, $file);
				my $outcome = unlink($file) ? "OK" : "Impossible to delete the source file: '$!'\n";
				print $outcome;
			}
		}		
	}
}
