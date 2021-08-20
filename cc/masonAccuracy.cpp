#include <iostream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cassert>
#include "kxseq.h"
#include "taxonomy.h"
#include "misc.h"

using namespace std;
using namespace kx;

unordered_map<string, int> genomeCnt;
unordered_map<string, int> masonAns;
TaxDB db;
int cntS = 0, cntG = 0, cntF = 0;

int singleTaxid = -1;

void getSGF(int tid, int &sp, int &g, int &f) {
    sp = g = f = -1;
    while (tid != 1 && tid != 0) {
        if (db.taxInfo[tid].rank == "species") {
           sp = tid;
        } else if (db.taxInfo[tid].rank == "genus") {
            g = tid;
        } else if (db.taxInfo[tid].rank == "family") {
            f = tid;
            break;
        }
        tid = db.taxInfo[tid].parent;
    }
}

void readMasonAns(const char *fileName) {
    gzFile fp = string(fileName) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(fileName, "r");
    kxseq<gzFile> ks(fp);

    while (ks.read() >= 0) {
        // cerr << ks.name() << ' ' << getAccession(ks.name()) << endl;
        int tid = db.acc2tid[getAccession(ks.comment())];
        if (tid == 0) {
            string name = ks.name();

            // for art reads
            for (size_t i = 0; i < name.size(); ++i) {
                if (name[i] == '-' || name[i] == '|' || name[i] == '/') {
                    name.resize(i);
                    break;
                }
            }

            // cerr << getAccession(name) << endl;
            tid = db.acc2tid[getAccession(name)];


            if (name.substr(0, 4) == "Rose") { // WARNING, hardcode for metabenchmark datasets
                tid = 267671;
            }

            genomeCnt[getAccession(name)]++;
        }

        if (singleTaxid > 0) {
            tid = singleTaxid;
        }
        int s, g, f;
        getSGF(tid, s, g, f);
        cntS += s != -1;
        cntG += g != -1;
        cntF += f != -1;
        // if (s == -1 || g == -1 || f == -1) continue;

        masonAns[trimno(ks.name())] = tid;
    }

    // cerr << cntS << ' ' << cntG << ' ' << cntF << endl;

    gzclose(fp);
}

vector<int> splitTids(string &tidLable) {
    size_t start = 0, spliter;
    vector<int> ret;
    while ((spliter = tidLable.find_first_of(';', start)) != string::npos) {
        string label = tidLable.substr(start, spliter - start);
        if (!isdigit(label[0]))
            ret.push_back(db.acc2tid[getAccession(label)]);
        else
            ret.push_back(atoi(label.c_str()));
        start = spliter + 1;
    }

    string label = tidLable.substr(start);
    if (!isdigit(label[0]))
        ret.push_back(db.acc2tid[getAccession(label)]);
    else
        ret.push_back(atoi(label.c_str()));

    return ret;
}

struct AccuracyStat {
    int sCorrect, sFalse, sUnspec;
    int gCorrect, gFalse, gUnspec;
    int fCorrect, fFalse, fUnspec;

    unordered_map<string, pair<int, int> > statPerGenomeSp;
    unordered_map<string, pair<int, int> > statPerGenomeG;
    unordered_map<string, pair<int, int> > statPerGenomeF;


    AccuracyStat() {
        sCorrect = sFalse = gCorrect = gFalse = fCorrect = fFalse = 0;
        sUnspec = gUnspec = fUnspec = 0;
    }

    bool onOrBelowLevel(int tid, const string &level) {
        while (tid != 1 && tid != 0) {
            if (db.taxInfo[tid].rank == level) {
                return true;
            }
            tid = db.taxInfo[tid].parent;
        }
        return false;
    }

    int classResult(int masonTidAtLevel, int yourTidAtLevel, int yourTid, int lca) {
        if (masonTidAtLevel == yourTidAtLevel) return 1;
        if (lca == yourTid) return 2;
        return 3;
    }

    void compare(const string &genomeName, int masonTid, int yourTid) {
        int ms, mg, mf, ys, yg, yf;
        getSGF(masonTid, ms, mg, mf);
        getSGF(yourTid, ys, yg, yf);        

        int lca = db.LCA({masonTid, yourTid});
        // cerr << lca << endl;

        if (ms != -1) {
            int cr = classResult(ms, ys, yourTid, lca);
            sCorrect += cr == 1;
            sUnspec += cr == 2;
            sFalse += cr == 3;

            statPerGenomeSp[genomeName].first += cr == 1;
            statPerGenomeSp[genomeName].second += cr == 3;
        }

        if (mg != -1) {
            int cr = classResult(mg, yg, yourTid, lca);
            gCorrect += cr == 1;
            gUnspec += cr == 2;
            gFalse += cr == 3;

            statPerGenomeG[genomeName].first += cr == 1;
            statPerGenomeG[genomeName].second += cr == 3;
        }

        if (mf != -1) {
            int cr = classResult(mf, yf, yourTid, lca);
            fCorrect += cr == 1;
            fUnspec += cr == 2;
            fFalse += cr == 3;

            statPerGenomeF[genomeName].first += cr == 1;
            statPerGenomeF[genomeName].second += cr == 3;
        }
    }
};

void calcAccuracy(const char *fileName) {
    gzFile fp = string(fileName) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(fileName, "r");
    kxstream<gzFile> ks(fp);
    string buf;
    string name;
    string label;

    AccuracyStat s;

    while (ks.get_until(kSepSpace, name) != kEOF) {
        name = trimno(name);
        ks.get_until(kSepLine, label); // label
        if (label.size() == 0 || label == "*") continue;
        int lca = db.LCA(splitTids(label));
        if (masonAns.count(name)) {
            string genome = name;

            // for art reads
            for (size_t i = 0; i < genome.size(); ++i) {
                if (genome[i] == '-') {
                    genome.resize(i);
                    break;
                }
            }

            // int sf = s.fFalse;
            s.compare(getAccession(genome), masonAns[name], lca);
            // if (sf != s.fFalse) {
            //     cerr << name << endl;
            // }
        }
    }

    cout << "Species: " << s.sCorrect << ' ' << s.sFalse << ' ' << s.sCorrect * 100.0 / cntS << ' ' << s.sFalse * 100.0 / (s.sFalse + s.sCorrect + s.sUnspec)
        << ' ' << s.sUnspec * 100.0 / (s.sFalse + s.sCorrect + s.sUnspec) << endl;
    cout << "Genus: " << s.gCorrect << ' ' << s.gFalse << ' ' << s.gCorrect * 100.0 / cntG << ' ' << s.gFalse * 100.0 / (s.gFalse + s.gCorrect + s.gUnspec)
        << ' ' << s.gUnspec * 100.0 / (s.gFalse + s.gCorrect + s.gUnspec) << endl;
    cout << "Family: " << s.fCorrect << ' ' << s.fFalse << ' ' << s.fCorrect * 100.0 / cntF << ' ' << s.fFalse * 100.0 / (s.fFalse + s.fCorrect + s.fUnspec)
        << ' ' << s.fUnspec * 100.0 / (s.fFalse + s.fCorrect + s.fUnspec) << endl;
    cout << endl;

    cout << "Stat per genome (species level):" << endl;
    for (auto it = s.statPerGenomeSp.begin(); it != s.statPerGenomeSp.end(); ++it) {
        cout << it->first << ": " << it->second.first * 100.0 / genomeCnt[it->first] << ' ' << it->second.second * 100.0 / max(1, it->second.first + it->second.second) << endl;
    }
    cout << endl;

    cout << "Stat per genome (genus level):" << endl;
    for (auto it = s.statPerGenomeG.begin(); it != s.statPerGenomeG.end(); ++it) {
        cout << it->first << ": " << it->second.first * 100.0 / genomeCnt[it->first] << ' ' << it->second.second * 100.0 / max(1, it->second.first + it->second.second) << endl;
    }
    cout << endl;

    cout << "Stat per genome (family level):" << endl;
    for (auto it = s.statPerGenomeF.begin(); it != s.statPerGenomeF.end(); ++it) {
        cout << it->first << ": " << it->second.first * 100.0 / genomeCnt[it->first] << ' ' << it->second.second * 100.0 / max(1, it->second.first + it->second.second) << endl;
    }
    cout << endl;
}

int main(int argc, const char **argv) {
    if (argc < 6) {
        cerr << "Usage: " << argv[0] << " <acc2tidfile> <nodes.dmp> <names.dmp> <classify_file> <mason.fq> [singleTaxid]" << endl;
        cerr << "classify_file format (one record per line):" << endl;
        cerr << "read_id\ttaxid;taxid;taxid..." << endl << endl;
        cerr << "OR:" << endl;
        cerr << "read_id\tref;ref;ref..." << endl << endl;
        return 1;
    }

    if (argc > 6) {
        singleTaxid = atoi(argv[6]);
    }

    cerr << "Reading acc2tid" << endl;
    db.initAccToTax(argv[1]);
    cerr << "Read acc2tid done" << endl;

    // for (auto it = db.acc2tid.begin(); it != db.acc2tid.end(); ++it) {
    //     cerr << it->first << ' ' << it->second << endl;
    // }

    cerr << "Reading nodes and names" << endl;
    db.readNode(argv[2]);
    db.readName(argv[3]);
    cerr << "Read nodes and names done" << endl;

    cerr << "Reading mason fastq" << endl;
    readMasonAns(argv[5]);
    cerr << "Reading mason fastq done" << endl;

    calcAccuracy(argv[4]);
    return 0;
}