#include <stdio.h>
#include <string.h>
#include <errno.h>

void error(), syserr ();
   
/*  common error messages */
char  *errors[] = 
{
  "",
  "End of file encountered too soon.",
  "Not enough memory.",
  "Invalid option.",
  "Required option was not supplied."
};

int Last_error = 4;    /* the last index in  errors */


/****  error

  display error messages on  stderr  and terminate

  Input:
  msg  -  extra message
  ind  -  index in array  errors with common error messages

  for system error use   ind < 0
  for message only use   ind = 0
  for common error message use   0 < ind <= Last_error
  
  Format:
  ind < 0    ERROR:  msg  (errno;  system message)
  ind = 0    ERROR:  msg  
  0 < ind <= Last_error   ERROR:  msg  (ind:  common error message )
  ind > Last_error    ERROR: error - error index outside the range!
  
  error exits program with 1
*/

void
error( msg, ind )
     char msg[];
     int  ind;
{
   if ( ind < 0 )
     syserr( msg );
   else 
   { 
      if( ind == 0 )
	fprintf (stderr, "\nERROR: %s\n", msg );
      else {
	 if( ind <= Last_error )
	   fprintf (stderr,
		    "\nERROR: %s (%d: %s )\n", msg, ind, errors[ ind ] );
	 else
	   fprintf (stderr,
		    "\nERROR: error - error index outside the range!\n" );
      }
   }
   
   exit( 1 );
}


/***
  syserr:
  print system call error message and terminate 

*/
 
void  
syserr( msg )  
     char *msg;
{
      
   fprintf( stderr, "\nERROR: %s (%d", msg, errno );
   fprintf (stderr, "; %s)\n", strerror(errno));
}

