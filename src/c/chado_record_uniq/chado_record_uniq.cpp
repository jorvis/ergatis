#include <stdio.h>
#include <errno.h>

#include <vector>
#include <set>
#include <string>

using namespace std;

#define BUF_SIZE 4096
#define EOR_LEN 5

typedef vector<FILE *> TFiles;

void uniq(TFiles &files, FILE *out)
{
    typedef set<string> TWritten;
    TWritten written;
    bool newRecord = true;
    bool validRecord = true;
    char buf[BUF_SIZE];
    auto_ptr<string> prevBuf(new string());
    auto_ptr<string> currBuf(new string());

    //int counter = 0;

    for (TFiles::iterator i = files.begin(); i != files.end(); ++i) {

        //printf("Processing file %d\n", ++counter);

        while (fgets(buf, sizeof(buf), *i)) {
            if (strlen(buf) < EOR_LEN) {
                currBuf = prevBuf;
            }
            currBuf->append(buf);

            if (newRecord) {
                string cksum = currBuf->substr(0, 32);
                if (!written.count(cksum)) {
                    written.insert(cksum);
                    fprintf(out, "%s", buf);
                    validRecord = true;
                }
                else {
                    validRecord = false;
                }
            }
            else {
                if (validRecord) {
                    fprintf(out, "%s", buf);
                }
            }
            newRecord = currBuf->substr(currBuf->size() - 5) == "????\n";
            prevBuf = currBuf;
            currBuf = auto_ptr<string>(new string());
        }
    }
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
