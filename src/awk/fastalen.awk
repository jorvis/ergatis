#!/bin/awk -f
# Usage:  fastalen.awk
#   Print the lengths of the sequences in the multifasta file
#   read from stdin.

        {
         if  (substr ($1, 1, 1) == ">")
             {
              if  (len > 0)
                  printf "%s\t%d\n", tag, len;
              tag = substr ($1, 2);
              len = 0;
             }
           else
             {
              len += length ($1);
             }
        }

END     {
         if  (len > 0)
             printf "%s\t%d\n", tag, len;
        }
