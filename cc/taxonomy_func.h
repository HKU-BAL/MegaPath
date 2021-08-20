#ifndef TAXONOMY_H__
#define TAXONOMY_H__

#include <sstream>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cassert>
#include "kmseq.h"

// ========================== HELPER FUNCTIONS ============================

string removeVersion(const string &accesstion) {
    for (int i = accesstion.length() - 1; i >= 0; --i) {
        if (accesstion[i] == '.') {
            return accesstion.substr(0, i);
        } else if (!isdigit(accesstion[i])) {
            return accesstion;
        }
    }
    return accesstion;
}

string getAccession(const string &ntHeader) {
    int spliterPos = ntHeader.find_first_of('|');
    if (spliterPos == string::npos) {
        // new format header, see http://www.ncbi.nlm.nih.gov/news/09-17-2014-simple-FASTA-headers-genomes-FTP/
        return removeVersion(ntHeader);
    } else if (ntHeader.substr(0, spliterPos) == "gi") {
        // old header, skip the GI, and then the db flag
        spliterPos = ntHeader.find_first_of('|', spliterPos + 1);
        spliterPos = ntHeader.find_first_of('|', spliterPos + 1);
        int spliterPos2 = ntHeader.find_first_of('|', spliterPos + 1);

        return removeVersion(ntHeader.substr(spliterPos + 1, spliterPos2 - spliterPos - 1));
    } else {
        assert(false);
    }
}

// ========================== END OF HELPER FUNCTIONS ============================

struct TaxInfo {
    int parent;
    std::string rank;
    std::string name;
};

struct TaxDB {
    unordered_map<string, int> acc2tid;
    unordered_map<string, int> name2tid;
    vector<TaxInfo> taxInfo;
    static const int kInitNumTax = 2000000;

    TaxDB(): taxInfo(kInitNumTax) {}

    void initAccToTax(const char *acc2TidFile) {
        gzFile fp = string(acc2TidFile) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(acc2TidFile, "r");
        kmstream<gzFile> ks(fp);

        string buf, acc, tid;

        while (ks.get_until(kSepSpace, buf) != kEOF) {
            ks.get_until(kSepSpace, acc);
            ks.get_until(kSepSpace, tid);
            acc2tid[removeVersion(acc)] = atoi(tid.c_str());
            ks.get_until(kSepLine, buf);
        }
        gzclose(fp);
    }

    void readNode(const char *nodePath) {
        gzFile fpNode = string(nodePath) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(nodePath, "r");
        kmstream<gzFile> ksNode(fpNode);

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

        gzclose(fpNode);
    }

    void readName(const char *namePath, bool alsoAddName2tid) {
        gzFile fpName = string(namePath) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(namePath, "r");

        kmstream<gzFile> ksName(fpName);
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

            taxInfo[id].name = name;
            if (alsoAddName2tid) name2tid[name] = id;
        }

        gzclose(fpName);
    }
};

#endif