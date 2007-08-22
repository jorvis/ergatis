#include <stdio.h>
#include <errno.h>

#include <vector>
#include <set>
#include <string>
#include <unistd.h>

using namespace std;

#define BUF_SIZE 4096

typedef vector<FILE *> TFiles;

void uniq(TFiles &files, FILE *out, int dbtype)
{
    typedef set<string> TWritten;
    TWritten written;
    bool newRecord = true;
    bool validRecord = true;
    size_t buf_size = BUF_SIZE;
    char *buf = reinterpret_cast<char *>(malloc(sizeof(char) * buf_size));
    ssize_t num_read = 0;


    /*
      This is a hack on top of a hack.
      The program has been merely utilizing the size of the delimiter
      string value to measure and slice up the records (regardless
      what that particular value might be).

      This new hack works because all sybase delimiters are two strings
      long i.e.: field \0\t and record \0\n.  Whereas, all postgresql
      delimiters are one character long i.e.: field \t and record \n.

      Need to consider revision in the future.
     */

    char sybaseDelim[] = {0, '\n'};
    char postgresqlDelim[] = {'\n'};

    if (!buf) {
        perror("Cannot allocate buffer");
        exit(1);
    }
    for (TFiles::iterator i = files.begin(); i != files.end(); ++i) {
        while ((num_read = getline(&buf, &buf_size, *i)) >= 0) {
            if (newRecord) {
                string cksum(buf, buf + 32);
                if (!written.count(cksum)) {
                    written.insert(cksum);
                    fwrite(buf, sizeof(char), num_read, out);
                    validRecord = true;
                }
                else {
                    validRecord = false;
                }
            }
            else {
                if (validRecord) {
                    fwrite(buf, sizeof(char), num_read, out);
                }
            }

	    if (dbtype){
	      // postgresql
	      if (num_read >= sizeof(postgresqlDelim)) {
                newRecord = memcmp((buf + num_read) - sizeof(postgresqlDelim), 
                                   postgresqlDelim, sizeof(postgresqlDelim)) == 0;
	      }
	      else {
                newRecord = 0;
	      }
	    }
	    else {
	      // sybase
	      if (num_read >= sizeof(sybaseDelim)) {
                newRecord = memcmp((buf + num_read) - sizeof(sybaseDelim), 
                                   sybaseDelim, sizeof(sybaseDelim)) == 0;
	      }
	      else {
                newRecord = 0;
	      }
	    }
        }
    }
    free(buf);
}

int main(int argc, char **argv)
{
    TFiles files;
    FILE *out = stdout;
    bool help = false;
    int opt;

    /* Declaring dbtype in order to support TIGR-specific delimiters:
       sybase field delimiter is  \0\t
       sybase record delimiter is \0\n
       postgresql field delimiter is  \t
       postgresql record delimiter is \n
       Default dbtype=0 for sybase.
       If -d option specified, dbtype=1 and postgresql support is enabled.
       Default support is for sybase.
    */
    int dbtype=0;
    const char *OPTSTR = "i:o:h:d";
    while ((opt = getopt(argc, argv, OPTSTR)) != EOF) {
        switch (opt) {
        case 'i':
        {
            FILE *fp = fopen(optarg, "r");
            if (!fp) {
                perror(("Error accessing " + string(optarg)).c_str());
            }
            else {
                files.push_back(fp);
            }
            break;
        }
        case 'o':
            out = fopen(optarg, "w");
            if (!out) {
                perror(("Error writing " + string(optarg)).c_str());
            }
            break;
        case 'h':
            help = true;
            break;
	case 'd':
	  dbtype=1;
	  break;

        }
    }
    if (help) {
        fprintf(stderr, "usage: %s [-i <input_record_file> ...]\n"
                "\t[-o <output_uniqued_record_file>] [-h]\n",
                basename(argv[0]));
        return 1;
    }
    if (files.empty()) {
        files.push_back(stdin);
    }
    uniq(files, out, dbtype);
    for (TFiles::iterator i = files.begin(); i != files.end(); ++i) {
        if (*i != stdin) {
            fclose(*i);
        }
    }
    if (out && out != stdout) {
        fclose(out);
    }
}
