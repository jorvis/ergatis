#!/usr/bin/perl

#########################################################################
#                                                                       #
# Copyright 2012 Arthur Brady (abrady@umiacs.umd.edu).                  #
#                                                                       #
# Last revision 2012.10.26.1302.                                        #
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

my $inFile = shift;

my $qualRange = shift;

if ( $qualRange ne 'phred33' and $qualRange ne 'phred64' ) {
   
   system("touch $inFile.failed");

   die("Usage: $0 <fastq file> <phred33|phred64>\n");
}

my $qualHash = {};

if ( $qualRange eq 'phred64' ) {
   
   $qualHash = {
      
      ';' => -5,   '<' => -4,   '=' => -3,   '>' =>  -2,   '?' => -1,
      '@' =>  0,   'A' =>  1,   'B' =>  2,   'C' =>  3,   'D' =>  4,
      'E' =>  5,   'F' =>  6,   'G' =>  7,   'H' =>  8,   'I' =>  9,
      'J' => 10,   'K' => 11,   'L' => 12,   'M' => 13,   'N' => 14,
      'O' => 15,   'P' => 16,   'Q' => 17,   'R' => 18,   'S' => 19,
      'T' => 20,   'U' => 21,   'V' => 22,   'W' => 23,   'X' => 24,
      'Y' => 25,   'Z' => 26,   '[' => 27,   '\\' => 28,   ']' => 29,
      '^' => 30,   '_' => 31,   '`' => 32,   'a' => 33,   'b' => 34,
      'c' => 35,   'd' => 36,   'e' => 37,   'f' => 38,   'g' => 39,
      'h' => 40
   };

} else {
   
   $qualHash = {
      
      '!' =>  0,   '"' =>  1,   '#' =>  2,   '$' =>  3,   '%' =>  4,
      '&' =>  5,   '\'' =>  6,   '(' =>  7,   ')' =>  8,   '*' =>  9,
      '+' => 10,   ',' => 11,   '-' => 12,   '.' => 13,   '/' => 14,
      '0' => 15,   '1' => 16,   '2' => 17,   '3' => 18,   '4' => 19,
      '5' => 20,   '6' => 21,   '7' => 22,   '8' => 23,   '9' => 24,
      ':' => 25,   ';' => 26,   '<' => 27,   '=' => 28,   '>' => 29,
      '?' => 30,   '@' => 31,   'A' => 32,   'B' => 33,   'C' => 34,
      'D' => 35,   'E' => 36,   'F' => 37,   'G' => 38,   'H' => 39,
      'I' => 40,   'J' => 41,   'K' => 42,   'L' => 43,   'M' => 44
   };
}

my $outFNA = $inFile;

if ( $outFNA =~ /\./ ) {
   
   $outFNA =~ s/\.([^\.]+)$/.fna/;

   if ( $outFNA eq $inFile ) {
      
      system("touch $inFile.failed");

      die("You can't run this script on a FASTA file.\n");
   }

} else {
   
   $outFNA .= '.fna';
}

my $outQual = $inFile;

if ( $outQual =~ /\./ ) {
   
   $outQual =~ s/\.([^\.]+)$/.qual/;

} else {
   
   $outQual .= '.qual';
}

if ( not open IN, "<$inFile" ) {
   
   system("touch $inFile.failed");

   die("Can't open $inFile for reading.\n");
}

if ( not open FNA, ">$outFNA" ) {
    
   system("touch $inFile.failed");

   die("Can't open $outFNA for writing.\n");
}

if ( not open QUAL, ">$outQual" ) {
   
   system("touch $inFile.failed");

   die("Can't open $outQual for writing.\n");
}

my $lineCount = 0;

my $phase = 0;

while ( my $line = <IN> ) {
   
   chomp $line;
   
   $lineCount++;

   if ( $phase == 0 ) {
      
      if ( $line !~ /^\@/ ) {
         
         system("touch $inFile.failed");

         die("Expected header on line $lineCount; got \"$line\" instead.  Exiting.\n");

      } else {
         
         $line =~ s/^\@/>/;

         print FNA "$line\n";

         print QUAL "$line\n";
      }

   } elsif ( $phase == 1 ) {
      
      print FNA "$line\n";

   } elsif ( $phase == 3 ) {
      
      print QUAL &convert($line) . "\n";
   }

   $phase++;

   $phase %= 4;
}

system("touch $inFile.complete");

sub convert {
   
   my $inLine = shift;

   my $outLine = '';
   
   my $index = 0;
   
   while ( $index < length($inLine) ) {
      
      my $newChar = substr($inLine, $index, 1);

      $outLine .= $qualHash->{$newChar} . ' ';

      $index++;
   }

   $outLine =~ s/\s*$//;

   return $outLine;
}
