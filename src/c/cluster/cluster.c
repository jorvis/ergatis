#include <stdio.h>
#include <string.h>

#include "cluster_defs.h"
#include "cluster_lists_defs.h"
#include "cluster_vars.h"

int main ( argc, argv )
     int  argc;
     string argv[];
{
   void build_pair_list(), build_all_clusters();
   bool input();

   FIRST_GO_ROUND = 1;
   PAIR_LIST_FIRST = PAIR_LIST_LAST = NULL;

   if (input (argc, argv) == TRUE)
   {
      build_pair_list (&PAIR_LIST_FIRST, &PAIR_LIST_LAST);
      
      if (PAIR_LIST_FIRST != NULL && PAIR_LIST_LAST != NULL)
	build_all_clusters (PAIR_LIST_FIRST, PAIR_LIST_LAST);
      else printf ("\nEmpty list of pairs.\n");
      printf("\n");
   }
   return(0);
}

/*############################### BUILDING CLUSTERS */


/*
  build_all_clusters:
  cluster elements using list of pairs

  Input:
  filename - input file name
*/
void
build_all_clusters ( list_first, list_last )
     HITPTR list_first;
     HITPTR list_last;
{
   void    build_one_cluster();
   HITPTR  ptr;

   while (list_first != NULL)
   {
      /* select element */
      ptr = list_first; 

      /* add to cluster all elements that match the selected element */
      if (ptr->valid == 1)
	build_one_cluster (ptr, list_first);

      list_first = list_first-> next; /* take next pair */
      free (ptr);
   }
}


/*
  build_one_cluster:

  build cluster beginning from a given element: all elements matching
  the given one will be added to the cluster list

  For the more stringent version (-f option) the new element will
  be aded if it matches all other elements in the cluster as well
  (every element has to match every other one).
*/
void
build_one_cluster ( elem, pair_l_first )
     HITPTR elem, pair_l_first;
{
   HITPTR new_hit();
   HITPTR cluster_first, cluster_last, curr_el, ptr, biglist;
   int cluster_count=0;

   if (PRINTDATA) printf ("\n\n--> Cluster for: %s", elem->query);
   if (!PRINTDATA) {
       if (!FIRST_GO_ROUND) 
       {
	   printf("\n");
       }
       FIRST_GO_ROUND = 0;
       printf ("%s", elem->query);
   }
 
  /* initialize current cluster list */
   cluster_first = cluster_last = NULL;

   /* create first element, initiate the cluster */
   ptr = new_hit (elem-> query, "");
   add_obj (ptr, &cluster_first, &cluster_last);
   cluster_count++;

   curr_el = cluster_first;

   /* do it for every element in cluster */
   while (curr_el != NULL)
   {
      if (PRINTDATA) printf ("\nexpanding %s", curr_el->query); 
      
      /* select all sequences from big list which curr_el hits */
      biglist = pair_l_first;
      while (biglist != NULL)
      {
	 if  (biglist->valid == 1)
	   if (strcasecmp (biglist-> query, curr_el-> query) == 0) 
	   {
	      if (elem_belongs (biglist-> hit, cluster_first)== FALSE)
	      {
		 if (FAMILY_MODE==TRUE && cluster_count>1)
		   if (pairs_up_with_all (biglist-> hit,
					  cluster_first,
					  pair_l_first)==FALSE)
		     break;

		 if (PRINTDATA) printf ("\nadd: %s", biglist->hit);  
		 else printf (" %s", biglist->hit);  
		 
		 /* create element, add new element to the cluster */
		 ptr = new_hit (biglist-> hit, "");
		 add_obj (ptr, &cluster_first, &cluster_last);
		 cluster_count++;
	      }
	      /* ... mark element as "used" */
	      biglist-> valid = 0;
	   }
	   else /* check the hit name */
	     if (strcasecmp (biglist-> hit, curr_el-> query) == 0) 
	     {
		if (elem_belongs (biglist-> query, cluster_first)== FALSE)
		{
		   if (FAMILY_MODE==TRUE && cluster_count>1)
		     if (pairs_up_with_all (biglist-> query,
					    cluster_first,
					    pair_l_first)==FALSE)
		       break;

		   if (PRINTDATA) printf ("\nadd: %s", biglist->query); 
		   else printf (" %s", biglist->query); 
		   /* create element, add new element to the cluster */
		   ptr = new_hit (biglist-> query, "");
		   
		   add_obj (ptr, &cluster_first, &cluster_last);
		   cluster_count++;
		}
		/* ... mark element as "used" */
		biglist-> valid = 0;
	     }
	 biglist = biglist-> next;
      }
      /* move to the next element in current cluster */
      curr_el = curr_el-> next;
   }

   /* at this point our list contains names of sequences in cluster */
   if (PRINTDATA)
   {
      ptr = cluster_first;
      printf ("\n");
      while (ptr != NULL)
      {
	 printf ("%s ", ptr-> query);
	 ptr = ptr-> next;
      }
      printf ("\n");
   }
   /* free up current cluster list */
   free_list (&cluster_first, &cluster_last); 
}


/*
  buld_pair_list:

  build a list of pairs reading pairs from stdin
*/

void
build_pair_list ( list_first, list_last )
     HITPTR *list_first;
     HITPTR *list_last;
{
   HITPTR  new_hit();
   char    line[ LINELEN ], queryname[ ACCLEN ], hitname[ ACCLEN ];
   HITPTR  ptr;

   /**** read pairs of query-hit from stdin */
   while (fgets (line, LINELEN, stdin) != NULL ) 
   {
      sscanf (line,"%s%s", queryname, hitname); 
      queryname[ strlen(queryname)] = EOS;
      hitname[ strlen(hitname)] = EOS;

      ptr = new_hit (queryname, hitname);
      add_obj (ptr, list_first, list_last);
   }
}


