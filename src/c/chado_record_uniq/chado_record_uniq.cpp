#include <stdio.h>
#include <errno.h>

#include <cstdlib>
#include <cstring>

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

    char sybaseDelim[] = {0, '\n'};
    char postgresqlDelim[] = {'\n'};
    char mysqlDelim[] = {'\n'};
    
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

            if (dbtype == 1){
                // postgresql
                if (num_read >= sizeof(postgresqlDelim)) {
                    newRecord = memcmp((buf + num_read) - sizeof(postgresqlDelim), 
                                   postgresqlDelim, sizeof(postgresqlDelim)) == 0;
                }
                else {
                    newRecord = 0;
                }
            } 
            else if (dbtype == 2) {
                // mysql
                if (num_read >= sizeof(mysqlDelim)) {
                    newRecord = memcmp((buf + num_read) - sizeof(mysqlDelim), 
                                   mysqlDelim, sizeof(mysqlDelim)) == 0;
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
       mysql field delimiter is  \t
       mysql record delimiter is \n
       dbtypes:
            0 - sybase (default)
            1 - postgres
            2 - myqsl
    */
    int dbtype=0;
    const char *OPTSTR = "i:o:d:h";
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
	case 'd':
	  {

	    char sybase[] = "sybase";
	    char postgresql[] = "postgresql";
        char mysql[] = "mysql";
	    
	    if (strcmp(optarg,sybase) == 0){
	      dbtype=0;
	    }
	    else if (strcmp(optarg,postgresql) == 0){
	      dbtype=1;
	    }
	    else if (strcmp(optarg,mysql) == 0){
	      dbtype=2;
	    }
	    else {
	      fprintf(stderr, "Invalid database_type '%s'.  Supported types are sybase, mysql and postgresql.\n", optarg);
	      exit(10);
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
