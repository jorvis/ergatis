#include <stdio.h>
#include <errno.h>

#include <vector>
#include <set>
#include <string>
#include <unistd.h>

using namespace std;

#define BUF_SIZE 4096

typedef vector<FILE *> TFiles;

void uniq(TFiles &files, FILE *out)
{
    typedef set<string> TWritten;
    TWritten written;
    bool newRecord = true;
    bool validRecord = true;
    size_t buf_size = BUF_SIZE;
    char *buf = reinterpret_cast<char *>(malloc(sizeof(char) * buf_size));
    ssize_t num_read = 0;
    char delim[] = {0, '\n'};
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
            if (num_read >= sizeof(delim)) {
                newRecord = memcmp((buf + num_read) - sizeof(delim), 
                                   delim, sizeof(delim)) == 0;
            }
            else {
                newRecord = 0;
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
    const char *OPTSTR = "i:o:h";
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
    uniq(files, out);
    for (TFiles::iterator i = files.begin(); i != files.end(); ++i) {
        if (*i != stdin) {
            fclose(*i);
        }
    }
    if (out && out != stdout) {
        fclose(out);
    }
}
