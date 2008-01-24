#!/bin/awk -f

BEGIN   {
         if  (ARGC < 3)
             Usage_Exit();
 
         if  (MAX_GENE_LEN == 0)
             MAX_GENE_LEN = 100000;
 
         len = ARGV [1];
         delete ARGV [1];
 
         sep = ARGV [2];
         delete ARGV [2];
        }
 
 
        {
         if  (1 * $2 < $3)
             {
              gene_len = 1 + $3 - $2;
              dir = 1;
             }
           else
             {
              gene_len = 1 + $2 - $3;
              dir = -1;
             }
         if  (gene_len > MAX_GENE_LEN)
             dir *= -1;
 
         printf "%s %8d %8d\n", $1, $2 - dir * (sep + len),
              $2 - dir * (sep + 1);
        }
