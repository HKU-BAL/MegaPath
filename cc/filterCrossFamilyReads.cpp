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

struct TaxInfo {
    int parent;
    string rank;
    string name;
};

unordered_map<string, int> acc2tid;
vector<TaxInfo> taxInfo(2000000);

string getAccession(const string &ntHeader) {
    int spliterPos = ntHeader.find_first_of('|');
    if (spliterPos == string::npos) {
        // new format header, see http://www.ncbi.nlm.nih.gov/news/09-17-2014-simple-FASTA-headers-genomes-FTP/
        return ntHeader;
    } else if (ntHeader.substr(0, spliterPos) == "gi") {
        // old header, skip the GI, and then the db flag
        spliterPos = ntHeader.find_first_of('|', spliterPos + 1);
        spliterPos = ntHeader.find_first_of('|', spliterPos + 1);
        int spliterPos2 = ntHeader.find_first_of('|', spliterPos + 1);

        return ntHeader.substr(spliterPos + 1, spliterPos2 - spliterPos - 1);
    } else {
        assert(false);
    }
}

void initAccToTax(const char *filepath) {
    gzFile fp = string(filepath) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(filepath, "r");
    kxstream<gzFile> ks(fp);

    string buf, acc, tid;

    while (ks.get_until(kSepSpace, buf) != kEOF) {
        ks.get_until(kSepSpace, acc);
        ks.get_until(kSepSpace, tid);
        acc2tid[acc] = atoi(tid.c_str());
        ks.get_until(kSepLine, buf);
    }
    gzclose(fp);
}

void readNodeAndName(const char *nodepath, const char *namepath) {
    gzFile fpNode = string(nodepath) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(nodepath, "r");
    gzFile fpName = string(namepath) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(namepath, "r");

    kxstream<gzFile> ksNode(fpNode);
    kxstream<gzFile> ksName(fpName);

    string tid, buf, parent, rank;

    while (ksNode.get_until(kSepSpace, tid) != kEOF) {
        ksNode.get_until(kSepSpace, buf);
        ksNode.get_until(kSepSpace, parent);
        ksNode.get_until(kSepSpace, buf);
        ksNode.get_until(kSepSpace, rank);

        int id = atoi(tid.c_str());
        if (id >= taxInfo.size()) {
            taxInfo.resize((id + 1) * 1.5);
        }

        taxInfo[id].parent = atoi(parent.c_str());
        taxInfo[id].rank = rank;

        ksNode.get_until(kSepLine, buf);
    }

    // cerr << "Read nodes done" << endl;

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

        taxInfo[id].name = name;
    }

    gzclose(fpNode);
    gzclose(fpName);
}

void splitAcc(string &hits, vector<string> &accessions) {
    if (hits == "*" || hits.length() == 0) return;

    int start = 0;
    int next = hits.find_first_of(';');

    while (next != string::npos) {
        start = hits.find_first_of(',', start) + 1;
        accessions.push_back(getAccession(hits.substr(start, next - start)));

        // from ',' on, find next ';'
        start = next + 1;
        next = hits.find_first_of(';', start);
    }

    start = hits.find_first_of(',', start) + 1;
    accessions.push_back(getAccession(hits.substr(start, hits.size() - start)));
}

int popUpToLevel(int tid, const string &level) {
    int sp_id = tid;
    while (tid != 1 && tid != 0 && taxInfo[tid].rank != level) {
        if (taxInfo[tid].rank == "species") {
            sp_id = tid;
        }
        tid = taxInfo[tid].parent;
    }
    return tid <= 1 ? sp_id : tid;
}

int main(int argc, const char **argv) {
    if (argc < 5) {
        cerr << "Usage: " << argv[0] << " <acc2tidfile> <nodes.dmp> <names.dmp> <in.lsam> [filterLevel=family]" << endl;
        return 1;
    }

    cerr << "Reading acc2tid" << endl;
    initAccToTax(argv[1]);
    cerr << "Read acc2tid done" << endl;

    cerr << "Reading nodes and names" << endl;
    readNodeAndName(argv[2], argv[3]);
    cerr << "Read nodes and names done" << endl;

    string level = "family";
    if (argc > 5) {
        level = argv[5];
    }


    // 1-pass only
    gzFile fpLsam = string(argv[4]) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(argv[4], "r");
    kxstream<gzFile> ks(fpLsam);

    string readname, hits;
    string dump;
    string flag, score, seq, qual;

    while (ks.get_until(kSepSpace, readname) != kEOF) {
        ks.get_until(kSepSpace, flag);
        ks.get_until(kSepSpace, score);
        ks.get_until(kSepSpace, seq);
        ks.get_until(kSepSpace, qual);
        ks.get_until(kSepLine, hits);

        vector<string> accessions;
        splitAcc(hits, accessions);

        vector<int> tax;
        for (vector<string>::iterator itr = accessions.begin(); itr != accessions.end(); itr++) {
            if (acc2tid.find(*itr) != acc2tid.end()) {
                tax.push_back(popUpToLevel(acc2tid[*itr], level));
            }
        }

        sort(tax.begin(), tax.end());
        tax.resize(unique(tax.begin(), tax.end()) - tax.begin());

        if (tax.size() <= 1) {
            printf("@%s\n%s\n+\n%s\n", readname.c_str(), seq.c_str(), qual.c_str());
        }
    }

    gzclose(fpLsam);

    return 0;
}