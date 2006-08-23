#include <stdio.h>
#include <errno.h>

#include <vector>
#include <set>
#include <string>

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
    auto_ptr<vector<char> > prevBuf(new vector<char>());
    auto_ptr<vector<char> > currBuf(new vector<char>());
    char delim[] = {0, '\n'};
    if (!buf) {
        perror("Cannot allocate buffer");
        exit(1);
    }
    for (TFiles::iterator i = files.begin(); i != files.end(); ++i) {
        while ((num_read = getline(&buf, &buf_size, *i)) >= 0) {
            if (num_read < static_cast<ssize_t>(sizeof(delim))) {
                currBuf = prevBuf;
            }
            copy(buf, buf + num_read, back_inserter(*currBuf));
            if (newRecord) {
                string cksum(currBuf->begin(), currBuf->begin() + 32);
                if (!written.count(cksum)) {
                    written.insert(cksum);
                    fwrite(&currBuf->at(0), sizeof(char),
                           currBuf->size(), out);
                    validRecord = true;
                }
                else {
                    validRecord = false;
                }
            }
            else {
                if (validRecord) {
                    fwrite(&currBuf->at(0), sizeof(char),
                           currBuf->size(), out);
                }
            }
            newRecord = memcmp(&*(currBuf->end() - sizeof(delim)),
                               delim, sizeof(delim)) == 0;
            prevBuf = currBuf;
            currBuf = auto_ptr<vector<char> >(new vector<char>());
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
