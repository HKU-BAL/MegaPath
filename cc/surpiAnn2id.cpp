#include <sstream>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cassert>
#include "kxseq.h"

using namespace std;
using namespace kx;

unordered_map<string, int> name2tid;

void readName(const char *namepath) {
    gzFile fpName = string(namepath) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(namepath, "r");
    kxstream<gzFile> ksName(fpName);

    string buf;

    while (ksName.get_until(kSepLine, buf) != kEOF) {
        if (buf.find("scientific name") == string::npos) {
            continue;
        }

        istringstream is(buf);
        int id;
        char d;
        is >> id >> d;
        string name;
        string nextpart;

        is >> nextpart;

        while (nextpart != "|") {
            if (name.length() > 0) name += "_";
            name += nextpart;
            is >> nextpart;
        }

        name2tid[name] = id;
    }

    gzclose(fpName);
}

int main(int argc, const char **argv) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <names.dmp> <surpi.annotation>" << endl;
        return 1;
    }

    readName(argv[1]);

    gzFile fp = string(argv[2]) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(argv[2], "r");

    kxstream<gzFile> ks(fp);
    string buf;
    while (ks.get_until('\t', buf) != kEOF) {
        for (int i = buf.length() - 1; i >= 0; --i) {
            if (buf[i] == '#') {
                buf = buf.substr(0, i);
                break;
            }
        }
        cout << buf << '\t';

        while (ks.get_until('\t', buf) != kEOF) {
            if (buf.length() >= 9 && buf.substr(0, 9) == "species--") {
                for (int i = 0; i < buf.size(); ++i) {
                    if (isspace(buf[i])) {
                        buf[i] = '_';
                    }
                }
                cout << name2tid[buf.substr(9)] << '\n';
                break;
            }
        }
        ks.get_until(kSepLine, buf);
    }

    gzclose(fp);

    return 0;
}