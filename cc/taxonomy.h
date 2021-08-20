#ifndef TAXONOMY_H__
#define TAXONOMY_H__

#include <sstream>
#include <unordered_map>
#include <vector>
#include <cassert>
#include <cctype>
#include <zlib.h>
#include "kxseq.h"

// ========================== HELPER FUNCTIONS ============================

inline std::string removeVersion(const std::string &accesstion) {
    for (int i = accesstion.length() - 1; i >= 0; --i) {
        if (accesstion[i] == '.') {
            return accesstion.substr(0, i);
        } else if (!isdigit(accesstion[i])) {
            return accesstion;
        }
    }
    return accesstion;
}

inline std::string getAccession(const std::string &ntHeader) {
    size_t spliterPos = ntHeader.find_first_of('|');
    if (spliterPos == std::string::npos) {
        // new format header, see http://www.ncbi.nlm.nih.gov/news/09-17-2014-simple-FASTA-headers-genomes-FTP/
        return removeVersion(ntHeader);
    } else if (ntHeader.substr(0, spliterPos) == "gi") {
        // old header, skip the GI, and then the db flag
        spliterPos = ntHeader.find_first_of('|', spliterPos + 1);
        spliterPos = ntHeader.find_first_of('|', spliterPos + 1);
        size_t spliterPos2 = ntHeader.find_first_of('|', spliterPos + 1);
        return removeVersion(ntHeader.substr(spliterPos + 1, spliterPos2 - spliterPos - 1));
    } else {
        return removeVersion(ntHeader);
    }
}

// ========================== END OF HELPER FUNCTIONS ============================

struct TaxInfo {
    int parent;
    std::string rank;
    std::string name;
};

struct TaxDB {
    static const int kInitNumTax = 2000000;

    std::unordered_map<std::string, int> acc2tid;
    std::unordered_map<std::string, int> name2tid;
    std::vector<TaxInfo> taxInfo;

    TaxDB(): taxInfo(kInitNumTax) {}

    void initAccToTax(const char *acc2TidFile) {
        gzFile fp = std::string(acc2TidFile) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(acc2TidFile, "r");
        kx::kxstream<gzFile> ks(fp);
        std::string buf, acc, tid;

        while (ks.get_until(kx::kSepSpace, buf) != kx::kEOF) {
            ks.get_until(kx::kSepSpace, acc);
            ks.get_until(kx::kSepSpace, tid);
            acc2tid[removeVersion(acc)] = atoi(tid.c_str());
            ks.get_until(kx::kSepLine, buf);
        }
        gzclose(fp);
    }

    void readNode(const char *nodePath) {
        gzFile fpNode = std::string(nodePath) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(nodePath, "r");
        kx::kxstream<gzFile> ksNode(fpNode);
        std::string tid, buf, parent, rank;

        while (ksNode.get_until(kx::kSepSpace, tid) != kx::kEOF) {
            ksNode.get_until(kx::kSepSpace, buf);
            ksNode.get_until(kx::kSepSpace, parent);
            ksNode.get_until(kx::kSepSpace, buf);
            ksNode.get_until(kx::kSepSpace, rank);

            int id = atoi(tid.c_str());
            if (id >= (int)taxInfo.size()) {
                taxInfo.resize((id + 1) * 1.5);
            }

            taxInfo[id].parent = atoi(parent.c_str());
            taxInfo[id].rank = rank;

            ksNode.get_until(kx::kSepLine, buf);
        }

        gzclose(fpNode);
    }

    void readName(const char *namePath, bool alsoAddName2tid = false) {
        gzFile fpName = std::string(namePath) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(namePath, "r");
        kx::kxstream<gzFile> ksName(fpName);
        std::string buf;

        while (ksName.get_until(kx::kSepLine, buf) != kx::kEOF) {
            if (buf.find("scientific name") == std::string::npos) {
                continue;
            }

            std::stringstream is(buf);
            int id;
            char d;
            is >> id >> d;
            std::string name;
            std::string nextpart;

            is >> nextpart;

            while (nextpart != "|") {
                if (name.length() > 0) name += " ";
                name += nextpart;
                is >> nextpart;
            }

            taxInfo[id].name = name;
            if (alsoAddName2tid) name2tid[name] = id;
        }

        gzclose(fpName);
    }

    int popUpToSpecies(int tid) {
	    while (tid != 1 && tid != 0 && taxInfo[tid].rank != "species") {
	        tid = taxInfo[tid].parent;
	    }
	    return tid;
    }

	std::string getTaxInfo(int tid) {
	    std::string s;
	    while (tid != 1 && tid != 0) {
	        if (taxInfo[tid].rank == "species" ||
	            taxInfo[tid].rank == "genus" ||
	            taxInfo[tid].rank == "family" ||
	            taxInfo[tid].rank == "superkingdom") {

	            if (s.length()) { s += "|"; }
	            s += taxInfo[tid].rank + "|" + taxInfo[tid].name;
	        }
	        tid = taxInfo[tid].parent;
	    }
	    return s;
	}

	int LCA(const std::vector<int> &tids) {
	    if (tids.size() == 1) return tids[0];
	    std::vector<std::vector<int>> lineage(tids.size());

	    for (size_t i = 0; i < tids.size(); ++i) {
	        int tid = tids[i];
	        while (tid != 0 && tid != 1) {
	            lineage[i].push_back(tid);
	            tid = taxInfo[tid].parent;
	        }
	        lineage[i].push_back(tid);
	    }

	    int lca = 0;
	    for (size_t k = 0; ; ++k) {
	        int lcaCandidate = *(lineage[0].rbegin() + k);
	        for (size_t i = 1; i < tids.size(); ++i) {
	            if (lineage[i].size() < k + 1 || *(lineage[i].rbegin() + k) != lcaCandidate) {
	                return lca;
	            }
	        }
	        lca = lcaCandidate;
	    }

	    return lca;
	}
};

#endif