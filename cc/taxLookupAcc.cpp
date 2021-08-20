#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <assert.h>
#include <unistd.h>
#include "kxseq.h"
#include "taxonomy.h"
#include "misc.h"

using namespace std;
using namespace kx;

TaxDB db;

vector<int> getSuperKingdom(const map<int, double> &taxid_and_score) {
    vector<int> ret;
    for (auto it = taxid_and_score.begin(); it != taxid_and_score.end(); ++it) {
        int tid = it->first;
        while (tid != 1 && tid != 0) {
            if (db.taxInfo[tid].rank == "superkingdom") {
                ret.push_back(tid);
                break;
            }
            tid = db.taxInfo[tid].parent;
        }
    }
    sort(ret.begin(), ret.end());
    ret.resize(unique(ret.begin(), ret.end()) - ret.begin());
    return ret;
}

int main(int argc, const char **argv) {
    if (argc < 5) {
        cerr << "Usage: " <<  argv[0] << " <acc2tidfile> <nodes.dmp> <names.dmp> <in.lsam>" << endl;
        return 1;
    }

    cerr << "Reading acc2tid" << endl;
    db.initAccToTax(argv[1]);
    cerr << "Read acc2tid done" << endl;

    cerr << "Reading nodes" << endl;
    db.readNode(argv[2]);
    db.readName(argv[3]);
    cerr << "Read nodes done" << endl;

    // 1-pass only

    gzFile fpLsam = string(argv[4]) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(argv[4], "r");
    kxstream<gzFile> ks(fpLsam);
    string buf;

    while (ks.get_until(kSepLine, buf) != kEOF) {
        vector<string> rec = splitBy(buf, '\t');
        assert(rec.size() >= 6);

        for (int i = 0; i < 5; ++i) {
            cout << rec[i] << '\t';
        }

        vector<pair<string, double> > accessions = splitAcc(rec[5]);
        map<int, double> taxid_and_score;

        for (auto itr = accessions.begin(); itr != accessions.end(); itr++) {
            auto dbitr = db.acc2tid.find(removeVersion(itr->first));
            if (dbitr != db.acc2tid.end()) {
                int sp = db.popUpToSpecies(dbitr->second);
                taxid_and_score[sp] = max(taxid_and_score[sp], itr->second);
            } else {
                cerr << "Error: Taxid not found for " << itr->first << endl;
            }
        }

        if (taxid_and_score.size() == 0) {
            cout << '*';
        } else {
            for (auto it = taxid_and_score.begin(); it != taxid_and_score.end(); ++it) {
                if (it != taxid_and_score.begin()) cout << ";";
                cout << it->second << ',' << it->first;
            }
        }

        for (size_t i = 6; i < rec.size(); ++i) {
            cout << '\t' << rec[i];
        }

        vector<int> sk = getSuperKingdom(taxid_and_score);
        for (size_t i = 0; i < sk.size(); ++i) {
            cout << '\t' << db.taxInfo[sk[i]].name;
        }

        cout << '\n';
    }
    gzclose(fpLsam);

    return 0;
}