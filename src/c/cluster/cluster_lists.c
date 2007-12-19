#include <stdio.h>
#include <string.h>

#include "cluster_defs.h"
#include "cluster_lists_defs.h"

/*###################### procedures for LISTS */

/* ----------------- new_hit ()
   returns a pointer to a new hit structure
*/

HITPTR
new_hit ( queryname, hitname )
     char * queryname;
     char * hitname;
{
   HITPTR obj;
   
   obj = (HITPTR) malloc (sizeof (HitStruct));
   if  (obj != NULL)
   {
      strcpy (obj-> query, queryname);
      strcpy (obj-> hit, hitname);
      obj-> valid = 1;
      obj-> next = NULL;
   }
   else
   {
      printf ("new_hit (malloc) Not enough memory.");
      exit(1);
   }
   return (obj);
}

/* ----------------- add_obj
   add a new Hit structure to a one-way list lbeg...lend
*/

void  
add_obj ( obj, lbeg, lend )
     HITPTR obj, *lend, *lbeg;
{
   if ( (*lend == NULL) && (*lbeg == NULL) )
     *lbeg = obj;
   else
   {
      (*lend)-> next = obj;
   }

   *lend = obj;
}

/* ----------------- free_list
   delete entire list lbeg...lend
*/

void  
free_list ( lbeg, lend )
     HITPTR *lend, *lbeg;
{
   HITPTR tmp;

   tmp = *lbeg;
   while  (tmp != NULL)
   {
      *lbeg = tmp-> next;
      free (tmp);
      tmp = *lbeg;
   }
   *lend = NULL;
}



/*
  elem_belongs:
  check if a given name is in the query field in a given list

  Scan the whole list from lbeg to NULL.
*/
bool
elem_belongs ( name, lbeg )
     char * name;
     HITPTR lbeg;
{
   bool out = FALSE;
   HITPTR tmp;

   tmp = lbeg;
   while (tmp != NULL)
   {
      if (strcasecmp (tmp-> query, name) == 0) /* found */
      {
	 out = TRUE;
	 break;
      }
      tmp = tmp-> next;
   }
   return (out);
}



/*
  pair_belongs:
  check if a given pair is in the list of pairs

  Scan the whole list from lbeg to NULL.
*/
bool
pair_belongs ( name1, name2, lbeg )
     char * name1;
     char * name2;
     HITPTR lbeg;
{
   bool out = FALSE;
   HITPTR tmp;

   tmp = lbeg;
   while (tmp != NULL)
   {
      if ((strcasecmp (tmp-> query, name1) == 0 &&
	   strcasecmp (tmp-> hit, name2) == 0) ||
	  (strcasecmp (tmp-> query, name2) == 0 &&
	   strcasecmp (tmp-> hit, name1) == 0))
      {
	 out = TRUE;
	 break;
      }
      tmp = tmp-> next;
   }
   return (out);
}


/*
  pairs_up_with_all:
  check if a given name is in a pair with each element in the list
  cluster_lbeg

  Scan the whole list from cluster_lbeg to NULL.
*/
bool
pairs_up_with_all ( name, cluster_lbeg, all_pairs_beg )
     char * name;
     HITPTR cluster_lbeg, all_pairs_beg;
{
   bool belongs=TRUE, out = TRUE;
   HITPTR tmp;

   tmp = cluster_lbeg;
   while (tmp != NULL && belongs)
   {
/*      printf ("\npair: %s, %s\n", tmp-> query, name); */
      if ((belongs=pair_belongs(tmp-> query, name, all_pairs_beg))==FALSE) 
      {
	 out = FALSE;
	 break;
      }
      tmp = tmp-> next;
   }
   return (out);
}




