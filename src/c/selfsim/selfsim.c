/****************************************************************
  Copyright (c) 1999 The Institute for Genomic Research
  Author: Steven L. Salzberg
  Take a genome in FASTA format, a window size, and a shift, and
  compute the chi-square statistic for each window.  Print
  out the results as <center of window> <chisquare>.
  Input is:  chi_square_windows <genome> windowsize shift
  Defaults are window=2000 and shift=1000
****************************************************************/

#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  <string.h>

FILE *GENOME_FILE;

main(argc, argv)
int argc;
char *argv[];
{
  char one_line[120];
  static char genome[28000000];   /* max genome size is 8 Mbases */
  char Ch;
  short int Base;
  long int window_size=2000, shift=1000, this_window;
  long int start, total_codons=0, GC_count;
  long int seq_length, i=0, j, k;
  int codon;
  long int codon_counts[64];
  long int observed[64];
  double expected[64];
  float codon_freqs[64];
  double chi_square;

  if (argc < 2) {
    fprintf(stderr, "Computes a chi-squared statistic on subsequences of DNA in a genome, comparing\neach window's trinucleotide composition to the rest of the genome\nUsage: chi_square <genome> <optional window-size> <optional shift> \n  Default window size is 2000.\n  The shift is how much to shift over for the next window;\n  the default is 1000, which mean the windows overlap.\n");
    exit(1);
  }
  if ((GENOME_FILE = fopen(argv[1], "r")) == NULL ) {
    fprintf(stderr, "Can't open genome file.\n");
    exit(1);
  } 
  if (argc >= 3) 
    window_size = atoi(argv[2]);
  if (argc == 4) 
    shift = atoi(argv[3]);

  /* read in the genome from argv[1].  The genome is in FASTA
     format - the first line is an identifier and other misc info,
     and each subsequence line has 60-80 DNA bases.  */
  fgets(one_line, 100, GENOME_FILE); /* ignore first line */
  while((Ch = fgetc(GENOME_FILE)) != EOF) {
    if  (Ch == '\n' || Ch == ' ' || Ch == '\t')
      continue;
    Ch = tolower (Ch);
    switch  (Ch)  {
      case  'a' : Base = 0;
	break;
      case  'c' : Base = 1;
	break;
      case  'g' : Base = 2;
	break;
      case  't' : Base = 3;
	break;
      default : Base = 4;
      }
    genome[i] = Base;
    i++;
  }
  genome[i] = '\0';  /* mark end of sequence */
  seq_length = i;

  for (i=0; i<64; i++) {
    codon_counts[i] = 0;
  }
  /* count codons in all 6 frames, both strands, for whole genome */
  for (i=0; i<seq_length-2; i++) {
    codon = 16*genome[i] + 4*genome[i+1] + genome[i+2];
    /* forward strand */
    if (codon < 64) /* only a,c,g,t codons counted */
      codon_counts[codon]++;
    /*reverse strand */
    codon = 16*(3-genome[i+2]) + 4*(3-genome[i+1]) + 3-genome[i];
    if (codon < 64)
      codon_counts[codon]++;
  }
  for (i=0; i<64; i++) {
    total_codons += codon_counts[i];
  }
  for (i=0; i<64; i++) {
    codon_freqs[i] = (float)codon_counts[i]/total_codons;
    /* number of expected occurrences of a codon is that codon's
       frequency times number of codons in the window.  This is
       window_size * 6 / 3, for six frames and 3 bases per codon. */
    expected[i] = (double) codon_freqs[i] * window_size * 2;
  }
  /* now compute the chi-square statistic for each window.  The
     statistic is simply Sum[observed-expected)^2/expected]
     summed over all 64 codons.  */
  printf("# Pos is the position of the window center\n");
  printf("# Pos \t Chi2\t GC\%\n");
  for  (start = 0; start <= seq_length - window_size;  start += shift) {
    /* for each window within the genome */
    for (i=0; i<64; i++) 
	observed[i] = 0;
    GC_count = 0;
    for  (i = 0; i < window_size;  i++) {
      j = i+start;
      codon = 16*genome[j] + 4*genome[j+1] + genome[j+2];
      if (codon < 64) /* only a,c,g,t codons counted */
	observed[codon]++;
      codon = 16*(3-genome[j+2]) + 4*(3-genome[j+1]) + 3-genome[j];
      if (codon < 64) /* only a,c,g,t codons counted */
	observed[codon]++;
      /* keep track of GC content too */
      if (genome[j] == 1 || genome[j] == 2) 
	GC_count++;
    }
    chi_square = 0.0;
    for (j=0; j<64; j++) {
      chi_square += pow((double)observed[j] - expected[j],2)/expected[j];
    }
    /* print out results */
    printf("%ld\t %4.0f\t %1.1f\n",
	   start+(window_size/2),chi_square,(float)GC_count*100.0/window_size);
  }
}
