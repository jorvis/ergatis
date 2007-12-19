#include <stdio.h>
#include <string.h>
#include <malloc.h>

#include "cluster_defs.h"
#include "cluster_inp_defs.h"
#include "cluster_inp_vars.h"

void error(), syserr ();

bool input();
FILE *openfile();

void
usage( argc, argv )
     int  argc;
     char  ** argv;
{
   printf("\nUsage:  %s %s\n", argv[ 0 ], PARAMS);
   fflush (stdout);
}

bool
input ( argc, argv )
     int  argc;
     char  ** argv;
{
   int  out = FALSE;

   if (process_command_line( argc, argv ) == TRUE)
     out = TRUE;

   return( out );
}


bool
process_command_line( argc, argv )
     int  argc;
     char  ** argv;
{
   int  out = FALSE;
   
   if (argc == 1) 
     if (input_stream() == TRUE)
     {
	out = TRUE;
     }
     else
       usage (argc, argv);
   else
     if (options (argc, argv) == TRUE)
       out = TRUE;
   
   return (out);
}

/*  
  options:

  options  processes options from the command line

  Input:
  argc  - number of parameters in the command line
  argv  - parameters

  returns TRUE if analysis was successful, FALSE in the opposite case
  places in optind the argv index of the next argument to be processed   
*/

options( argc, argv )
     int  argc;
     char  ** argv;
{
   int  c, out = TRUE;
   char *optstring = "fho:";
   extern char *optarg;
   extern int  optind, opterr;
 
   opterr = 0;   /* disable error messages from  getopt  */
   while (((c = getopt( argc, argv, optstring ) ) != -1 ) && out == TRUE)
     switch( c )
     {
      case 'o': 
	strcpy (OUTFILE, optarg);
	OUTPUT_TO_FILE = TRUE;
	break;

      case 'f':  
	FAMILY_MODE = TRUE;
	break;

      case 'h':
	usage( argc, argv );
	out = FALSE;
	break;

      case '?':
	usage( argc, argv );
	error( "", 3 );
	out = FALSE;
	break;
     }

   return( out );
}



/* 
   input_stream:
 
   see if there is anything given on stdin
*/
bool
input_stream ()
{
   int  out = FALSE;
   int c;

   if ((c=fgetc(stdin))!=EOF) 
   { 
      out = TRUE;
      ungetc (c, stdin);
   }
   
   return (out);
}

