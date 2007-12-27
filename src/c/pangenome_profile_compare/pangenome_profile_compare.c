/*
 * This script performs the pangenome analysis for all genes of a given genome.
 *
 * Input is two binary files, one encoding gene presence/absence profiles,
 * and another containing profile masks representing the comparison genome set.
 *
 * Format is:
 *
 *  unsigned long PROFILE_ARRAY_SIZE   //number of unsigned long values in array
 *  unsigned long PROFILE_FRAME_COUNT  //number of unsigned longs per profile (eg: profiles > 1 unsigned long are split between 2 or more)
 *  unsigned long PROFILE_ELEMENT_1_FRAME_1_OF_x    // where x = PROFILE_FRAME_COUNT
 *  unsigned long PROFILE_ELEMENT_?_FRAME_?_OF_? 
 *  ...
 * 
 * Contact:
 *  Brett Whitty
 *  bwhitty@jcvi.org
 *
 */

#include <stdlib.h>
#include <stdio.h>

#define MAX_ARRAY_SIZE 1000000
#define TWO(c) (0x1u << (c))
#define MASK(c) (((unsigned int)(-1)) / (TWO(TWO(c)) + 1u))
#define COUNT(x,c) ((x) & MASK(c)) + (((x) >> (TWO(c))) & MASK(c))

// bit count function
int parallel_bitcount (unsigned int n) {
    n = COUNT(n, 0) ;
    n = COUNT(n, 1) ;
    n = COUNT(n, 2) ;
    n = COUNT(n, 3) ;
    n = COUNT(n, 4) ;
//    n = COUNT(n, 5) ;   // for 64-bit integers 
    return (int) n ;
}

int main (int argc, char *argv[]) {

    FILE *infile;   //pointer to profile/mask input file
    unsigned int profiles[MAX_ARRAY_SIZE];  //array for storing profiles
    unsigned int masks[MAX_ARRAY_SIZE];     //array for storing masks
    unsigned long *profile_array_size;      //profile array size
    unsigned long *profile_frames;          //number of frames per profile
    unsigned long *mask_array_size;         //mask array size
    unsigned long *mask_frames;             //number of frames per mask
    unsigned long pa_val = 0;
    unsigned long pf_val = 0;
    unsigned long ma_val = 0;
    unsigned long mf_val = 0;    
    unsigned long gene_count = 0;   //gene count for genome
    unsigned long core_count = 0;   //core gene count
    unsigned long shared_count = 0; //shared gene count
    unsigned long new_count = 0;    //new gene count
    unsigned long bitsum;           //profile bit sum
    //unsigned long selfsum;
    unsigned long subject_genome_count;
    unsigned long i, j, k, m;       //iterators
    
    profile_array_size = &pa_val;
    profile_frames = &pf_val;
    mask_array_size = &ma_val;
    mask_frames = &mf_val;
    //unsigned int *profiles, *masks;
    
    if (argc < 3) {
        fprintf(stderr, "%s\n", "Pan-genome analysis profile comparison tool");
        fprintf(stderr, "%s\n\n", "Brett Whitty, JCVI, 2007 (bwhitty@jcvi.org)");
        fprintf(stderr, "%s\n", "usage:");
        fprintf(stderr, "%s ", argv[0]);
        fprintf(stderr, "%s\n", "/path/to/profiles.X.out /path/to/masks.X.out");
        return 1;
    } 
   
    //read profiles
    if ((infile = fopen(argv[1], "rb")) == NULL) { 
        fprintf(stderr, "ERROR: Couldn't open profiles file.\n");
        exit(EXIT_FAILURE);
    }
    
    fread(profile_array_size, sizeof(unsigned long), 1, infile);    
    fprintf(stderr, "%d %s\n", (int) *profile_array_size, "profile elements");

    fread(profile_frames, sizeof(unsigned long), 1, infile);
    fprintf(stderr, "%d %s\n", (int) *profile_frames, "profile frames");

//    profiles = (unsigned int *) calloc(((int) profile_array_size * (int) profile_frames), sizeof(unsigned int));

    fread(profiles, sizeof(unsigned int), (*profile_array_size * *profile_frames), infile);

    if (! fclose(infile) == 0) {
        fprintf(stderr, "ERROR: Couldn't close profiles file.\n");
        exit(EXIT_FAILURE);
    }
   
    infile = NULL;
    
    //read masks
    if ((infile = fopen(argv[2], "rb")) == NULL) { 
        fprintf(stderr, "ERROR: Couldn't open mask file.\n");
        exit(EXIT_FAILURE);
    }
   
    fread(mask_array_size, sizeof(unsigned long), 1, infile);    
    fprintf(stderr, "%d %s\n", (int) *mask_array_size, "mask elements");
    
    fread(mask_frames, sizeof(unsigned long), 1, infile);
    fprintf(stderr, "%d %s\n", (int) *mask_frames, "mask frames");
    
    if (*profile_frames != *mask_frames) { 
        fprintf(stderr, "ERROR: Frame count mismatch between profiles and masks.\n");
        exit(EXIT_FAILURE);
    }
    
    fread(masks, sizeof(unsigned int), (*mask_array_size * *mask_frames), infile);

    if (! fclose(infile) == 0) {
        fprintf(stderr, "ERROR: Couldn't close masks file.\n");
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "%s\n", "running analysis...");
    
    gene_count = *profile_array_size / *profile_frames;
    
    for (m = 0; m < *mask_array_size; m = m + *mask_frames) {
        
        subject_genome_count = 0;
        core_count = 0;
        shared_count = 0;
        new_count = 0;
        
        for (i = 0; i < *mask_frames; i++) {
            subject_genome_count += parallel_bitcount(masks[m + i]);
        }
       
        for (j = 0; j < *profile_array_size; j = j + *profile_frames) {
            bitsum = 0;
            for (k = 0; k < *profile_frames; k++) {
                bitsum += parallel_bitcount(masks[m + k] & profiles[j + k]);
            }
            if (bitsum == subject_genome_count) {
                core_count++;
            }
            if (bitsum > 1) {
                shared_count++;
            }
            if (bitsum == 1) {
                new_count++;
            }
        }
        fprintf(stdout, "%lu\t", subject_genome_count);
        fprintf(stdout, "%lu\t%lu\t%lu\n", core_count, shared_count, new_count);
    }
    fprintf(stderr, "%s\n", "done.");

    return 0;
}

