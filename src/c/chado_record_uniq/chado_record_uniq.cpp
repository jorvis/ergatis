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
    for (TFiles::iterator i = files.begin(); i != files.end(); ++i) {
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

string create_read_pipe(const char *file)
{
    string pipe = "perl -pe 's/\\0\\n/\?\?\?\?\\n/g; "
        "s/\\0\\t/\?\?\\t\?\?/g;' ";
    return pipe + file;
}

string create_write_pipe(const char *file)
{
    string pipe = "perl -pe 's/\\?\\?\\?\\?\\n/\\0\\n/g; "
        "s/\\?\\?\\t\\?\\?/\\0\\t/g;'";
    return pipe + " > " + file;
}

int main(int argc, char **argv)
{
    TFiles files;
    FILE *out = 0;
    bool help = false;
    int opt;
    const char *OPTSTR = "i:o:h";
    while ((opt = getopt(argc, argv, OPTSTR)) != EOF) {
        switch (opt) {
        case 'i':
        {
            FILE *fp = popen(create_read_pipe(optarg).c_str(), "r");
            if (!fp) {
                perror(("Error accessing " + string(optarg)).c_str());
            }
            else {
                files.push_back(fp);
            }
            break;
        }
        case 'o':
            out = popen(create_write_pipe(optarg).c_str(), "w");
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
        FILE *fp = popen(create_read_pipe("/dev/stdin").c_str(), "r");
        files.push_back(fp);
    }
    if (!out) {
        out = popen(create_write_pipe("/dev/stdout").c_str(), "w");
    }
    uniq(files, out);
    for (TFiles::iterator i = files.begin(); i != files.end(); ++i) {
        if (*i != stdin) {
            pclose(*i);
        }
    }
    pclose(out);
}
