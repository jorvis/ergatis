#include "Bsml2Asn1Converter.hpp"

#include <iostream>
#include <fstream>

#include <unistd.h>

using namespace std;

int main(int argc, char **argv)
{
    try {
        const char *in = "/dev/stdin";
        const char *out = "/dev/stdout";
        ifstream contact;
        bool update = false;
        bool help = false;
        bool genProdSet = false;
        int opt;
        const char *OPTSTR = "i:o:uc:gh";
        while ((opt = getopt(argc, argv, OPTSTR)) != EOF) {
            switch (opt) {
            case 'i':
                in = optarg;
                break;
            case 'o':
                out = optarg;
                break;
            case 'u':
                update = true;
                break;
            case 'c':
                contact.open(optarg);
                break;
            case 'g':
                genProdSet = true;
                break;
            case 'h':
                help = true;
                break;
            }
        }
        if (help) {
            cerr << "usage: " << basename(argv[0]) << " [-i <bsml_file>]\n"
                 << "       [-o <asn_file>] [-u] [-c <contact_info>] [-g] [-h]"
                 << "\n" << endl;
            cerr << "\tu: update submission [default - false]\n"
                 << "\tc: submission contact info\n"
                 << "\t   if supplied, overrides data in BSML\n"
                 << "\tg: create gen-prod-set instead of nuc-prot-set "
                 << "[default - false]\n" << endl;
            cerr << "\tNOTE: -u is not currently working since data is not "
                 << "correctly captured\n"
                 << "\t      in BSML" << endl;
            return 2;
        }
        CBsml2Asn1Converter converter;
        converter.Convert(in, out, contact.is_open() ? &contact : 0, update,
                          genProdSet);
    }
#ifndef _DEBUG
    catch (const exception &e) {
        cerr << e.what() << endl;
        return 1;
    }
#else
    catch (int) { }
#endif
    return 0;
}
