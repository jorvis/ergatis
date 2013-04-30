#!/usr/bin/perl

#########################################################################
#                                                                       #
# Copyright 2012 Arthur Brady (abrady@umiacs.umd.edu).                  #
#                                                                       #
# Last revision 2012.10.26.1404.                                        #
#                                                                       #
# This program is free software: you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by  #
# the Free Software Foundation, either version 3 of the License, or     #
# (at your option) any later version.                                   #
#                                                                       #
# This program is distributed in the hope that it will be useful,       #
# but WITHOUT ANY WARRANTY; without even the implied warranty of        #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
# GNU General Public License for more details.                          #
#                                                                       #
# You should have received a copy of the GNU General Public License     #
# along with this program.  If not, see <http://www.gnu.org/licenses/>. #
#                                                                       #
#########################################################################

use strict;

$| = 1;

my $filePrefix = shift;

my $libraryArg = shift;

if ( $libraryArg eq '' ) {
   
   system("touch $filePrefix.failed");

   die("Usage: $0 <input file prefix> <library name>\n");

} elsif ( not -e "$filePrefix.fna" or not -e "$filePrefix.qual" ) {
   
   system("touch $filePrefix.failed");

   die("Can't find one or both of (\"$filePrefix.fna\" and \"$filePrefix.qual\"); exiting.\n");
}

foreach my $inFile ( "$filePrefix.fna", "$filePrefix.qual" ) {
   
   if ( not open IN, "<$inFile" ) {
      
      system("touch $filePrefix.failed");
      
      die("Can't open $inFile for reading.\n");
   }

   my $outFile = $inFile;

   $outFile =~ s/\.([^\.]+)$/.newbler-ready.$1/;

   if ( not open OUT, ">$outFile" ) {
      
      system("touch $filePrefix.failed");
      
      die("Can't open $outFile for writing.\n");
   }

   while ( my $line = <IN> ) {
      
      if ( $line !~ /^>/ ) {
         
         print OUT $line;

      } else {
         
         chomp $line;

         $line =~ /^>(\S+)/;

         my $readID = $1;

         my $dir = '';

         my $template = '';

         if ( $readID =~ /\/1$/ ) {
            
            $dir = 'F';

            $template = $readID;

            $template =~ s/\/\d$//;

         } elsif ( $readID =~ /\/2$/ ) {
            
            $dir = 'R';

            $template = $readID;

            $template =~ s/\/\d$//;

         } elsif ( $line =~ /^>(\S+)\s+(\S+)/ ) {
            
            $template = $1;

            my $secondBlock = $2;

            if ( $secondBlock =~ /([12]):.:\d+:\S+/ ) {
               
               my $dirNum = $1;

               if ( $dirNum == 1 ) {
                  
                  $dir = 'F';

               } else {
                  
                  $dir = 'R';
               }
            }
         }

         if ( $dir eq '' ) {
            
            system("touch $filePrefix.failed");

            die("Couldn't parse mate ID from read header: \"$line\".  Exiting.\n");

         } else {
            
            print OUT ">$readID template=$template dir=$dir library=$libraryArg\n";
         }

      } # end if ( $line !~ /^> )

   } # end while ( input-line iterator )

   close OUT;

   close IN;

} # end foreach ( file selector )

system("touch $filePrefix.complete");

