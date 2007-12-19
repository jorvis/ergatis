typedef char * string;

#define EOS '\0'
#define PRINTDATA FALSE

#define PARAMS " [-f] [-o <outfile>] \
\n\t -f              every element of the cluster has to \
\n\t                 match every other one; this is an more stringent \
\n\t                 version than the transitive closure which is the default \
\n\t -o <outfile>    report file, stdout by default (not yet implemented)\
\n\t -h              print help message \
\n\t \
\nThe input comes from stdin in the format of space-separated pairs \
\nof elements, for example: \
\nD G \
\nD E \
\nA B \
\nA C \
\nG E \
\nshould produce the output:\
\nD G E \
\nA B C \
\nwhen run without the -f option, and: \
\nD G E \
\nA B \
\nA C \
\nwhen run with the -f option."
