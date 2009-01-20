#include "SnpClusterer.hpp"
#include "SnpDataType.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unistd.h>

#include <cstring>

using namespace std;

static const char *DELIM = "//";

void AddSnps(istream &in, CSnpClusterer &clusterer)
{
    typedef CSnpClusterer::TSnps TSnps;
    CSnpClusterer::TSnpData snps;
    string line;
    while (getline(in, line)) {
        if (line == DELIM) {
            if (snps.size() < 2) {
                continue;
            }
            clusterer.AddSnps(snps);
            snps.clear();
        }
        else {
            snps.push_back(line);
        }
    }
}

void
PrintClusters(ostream &out, CSnpClusterer &clusterer)
{
    CSnpClusterer::TClusters clusters;
    clusterer.GetClusters(clusters);
    CSnpClusterer::TClusters::size_type count = 0;
    for (CSnpClusterer::TClusters::const_iterator i = clusters.begin();
         i != clusters.end(); ++i) {
        const CSnpClusterer::TSnps &snps = *i;
        for (CSnpClusterer::TSnps::const_iterator j = snps.begin();
             j != snps.end(); ++j) {
            const CSnpClusterer::TSnp &snp = *j;
            out << *snp << endl;
        }
        if (count++ < clusters.size()) {
            out << DELIM << endl;
        }
    }
}

int main(int argc, char **argv)
{
    try {
        istream *in = &cin;
        ostream *out = &cout;
        bool help = false;
        int opt;
        const char *OPTSTR = "i:o:h";
        while ((opt = getopt(argc, argv, OPTSTR)) != EOF) {
            switch (opt) {
            case 'i':
                in = new ifstream(optarg);
                break;
            case 'o':
                out = new ofstream(optarg);
                break;
            case 'h':
                help = true;
                break;
            }
        }
        if (help) {
            cerr << "usage: " << basename(argv[0])
                 << " [-i <pairwise_snps_data>] [-o <clustered_snp_output>]\n";
            cerr << "\n";
            cerr << "\ti: tab delimited pairwise data containing:\n"
                 << "\t   seq_id, from, to, orientation\n"
                 << "\t   records should be delimited by //\n"
                 << "\t   [default - stdin]\n"
                 << "\to: output for merged clusters\n"
                 << "\t   [default - stdout]"
                 << endl;
            return 1;
        }
        CSnpClusterer clusterer;
        AddSnps(*in, clusterer);
        PrintClusters(*out, clusterer);
        if (in != &cin) {
            delete in;
        }
        if (out != &cout) {
            delete out;
        }
    }
    catch (const exception &e) {
        cerr << e.what() << endl;
        return 1;
    }
    return 0;
}
