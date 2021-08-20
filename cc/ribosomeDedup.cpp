#include <string>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <bitset>
#include <assert.h>
#include <unordered_map>
#include <algorithm>
#include <stdint.h>
#include <zlib.h>
#include "kxseq.h"

using namespace kx;
using namespace std;

const size_t kMinLenToKeep = 100;
const size_t kKmerSize = 32;
const uint64_t kMask = kKmerSize == 32 ? 0xFFFFFFFFFFFFFFFFULL : ((1ULL << 2 * kKmerSize) - 1);

uint nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

unordered_map<uint64_t, vector<uint64_t> > kmer2pos;
vector<string> seqs;
vector<string> names;

void addToHashTable(const string &s, size_t sid) {	
	if (s.length() < kMinLenToKeep) {
		return;
	}

	for (size_t i = 0; i + kKmerSize <= s.length(); i += kKmerSize) {
		uint64_t h = 0;
		for (size_t j = 0; j < kKmerSize; ++j) {
			uint b = min(3U, nst_nt4_table[(int)s[i + j]]);
			h = (h << 2) | b;
		}

		kmer2pos[h].push_back((sid << 32) | i);
	}
}

bool filtered(const string &s, size_t sid) {
	if (s.length() < kMinLenToKeep) {
		return true;
	}

	uint64_t h = 0;

	for (size_t i = 0; i < s.length() && i < kKmerSize * 2 - 1; i++) {
		uint b = min(3U, nst_nt4_table[(int)s[i]]);
		h = (h << 2) | b;
		h &= kMask;

		if (i < kKmerSize - 1) continue;

		auto p = kmer2pos.find(h);
		if (p != kmer2pos.end()) {
			for (auto it = p->second.begin(); it != p->second.end(); ++it) {
				size_t tPos = *it & 0XFFFFFFFFULL;
				if ((*it >> 32) == sid || tPos < (i + 1 - kKmerSize)) continue;
				size_t id = *it >> 32;
				size_t offset = tPos - (i + 1 - kKmerSize);

				if (offset + s.length() > seqs[id].length()) { continue; }
				if (seqs[id].length() == s.length() && id > sid) { continue; }
				if (seqs[id].substr(offset, s.length()) == s) { return true; }
			}
		}
	}

	return false;
}

char Comp(char c) {
    switch (c) {
    case 'A':
    case 'a':
        return 'T';

    case 'C':
    case 'c':
        return 'G';

    case 'G':
    case 'g':
        return 'C';

    case 'T':
    case 't':
        return 'A';

    default:
    	return c;
    }
}

string RevComp(const string &s) {
    string ret;

    for (unsigned i = 0; i < s.length(); ++i) {
        ret.push_back(Comp(s[s.length() - 1 - i]));
    }

    return ret;
}

bool filteredRev(const string &s, size_t sid) {
	return filtered(RevComp(s), sid);
}

int main(int argc, char **argv) {
	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " <ribosome.fa>" << endl;
		return 1;
	}

	gzFile fp = gzopen(argv[1], "r");
	kxseq<gzFile> ks(fp);

	size_t num0 = 0, bp0 = 0;
	while (ks.read() >= 0) {
		names.push_back(ks.name());
		seqs.push_back(ks.seq());
		addToHashTable(ks.seq(), num0);
		num0++;
		bp0 += ks.seq().length();
	}

	cerr << "Read " << num0 << " sequences, " << bp0 << "bp" << endl;
	gzclose(fp);

	size_t num1 = 0, bp1 = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+: num1, bp1)
	for (size_t i = 0; i < seqs.size(); ++i) {
		if (!filtered(seqs[i], i) && !filteredRev(seqs[i], i)) {
			num1++;
			bp1 += seqs[i].length();
			printf(">%s\n%s\n", names[i].c_str(), seqs[i].c_str());
		}
	}

	cerr << "Keep " << num1 << " sequences, " << bp1 << "bp" << endl;

	return 0;
}