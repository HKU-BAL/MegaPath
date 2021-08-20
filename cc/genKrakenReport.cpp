#include <cstdio>
#include <string>
#include <vector>
#include <set>
#include <assert.h>
#include <unistd.h>
#include <math.h>
#include <unordered_map>
#include "kxseq.h"
#include "misc.h"
#include "taxonomy.h"

using namespace std;
using namespace kx;

const static int kRoot = -2;
const static vector<string> Ranks = {"domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"};
//const static string AllowedChar = "[],. _()";
TaxDB db;

struct TaxCounter {

	int totalReads;
	unordered_map<int, int> tidCount;
	unordered_map<int, int> tidAccCount; // accumulative count 
	unordered_map<int, set<int> > sons;

	TaxCounter() : totalReads(0) {
	}

	string foo(string a) {
// !!! pavian will fail if the name contains single quote
		string b = "";
		for (int i=0;i<a.length();++i) 
			//if (a[i] >='a' && a[i] <='z' || a[i] >='A' && a[i] <='Z' || find(AllowedChar.begin(), AllowedChar.end(), a[i]) != AllowedChar.end())
			if (a[i] != '\'')
				b+=a[i];
		return b;
	}

	static char getRankCode(const string &rank) {
		if (rank=="superkingdom") 
			return 'D';
		else if (find(Ranks.begin(),Ranks.end(), rank) != Ranks.end()) 
			return rank[0]+'A'-'a';
		else 
			return '-';
	}

	void createLineage_(const vector<int> &tids) {
		for (int i=0;i<tids.size()-1;i++) sons[tids[i+1]].insert(tids[i]);
		sons[kRoot].insert(tids.back());
	}

	void add1read(const string &s) {
		++totalReads;
		if (s.length() == 0 || s == "*") {
			tidCount[0]++;
			return;
		}

		vector<int> tids;
		vector<string> tid_and_score = splitBy(s, ';');
		for (unsigned i = 0; i < tid_and_score.size(); ++i) {
			auto s_tid = splitBy(tid_and_score[i], ',');
			tids.push_back(atoi(s_tid[1].c_str()));
		}

		int lca = db.LCA(tids);
		if (lca==0) lca = 1; /// !!!
		tidCount[lca]++;
		tids.clear();
		while (lca != 1 && lca != 0) {
			tids.push_back(lca);
			tidAccCount[lca]++;
			lca = db.taxInfo[lca].parent;
		}
		tids.push_back(lca);
		tidAccCount[lca]++;
		createLineage_(tids);

	}

        struct cmp_ {
                unordered_map<int, int> *counter;
                cmp_(unordered_map<int, int> *counter) : counter(counter) {}
                bool operator() (int a, int b) {
                        return (*counter)[a]> (*counter)[b];
                }
        };

	void printTable_(int tid, int depth) {
		if (tid >= 0 && (tid & 0xC0000000) == 0) {
			printf("%6.2f\t", tidAccCount[tid] * 100.0 / totalReads);
			printf("%d\t", tidAccCount[tid]);
			printf("%d\t", tidCount[tid]);
			printf("%c\t", getRankCode(db.taxInfo[tid].rank));
			printf("%d\t", tid);
			printf("%d\t", depth);
			for (int i=0;i<depth;i++) printf("  ");
			printf("%s\n", foo(db.taxInfo[tid].name).c_str());
		}

		vector<int> mySons(sons[tid].begin(), sons[tid].end());
		sort(mySons.begin(), mySons.end(), cmp_(&tidAccCount));

		for (auto it = mySons.begin(); it != mySons.end(); ++it) printTable_(*it,depth+1);
	
	}

	void printTable() {
		printf("perc\tn-clade\tn-stay\tlevel\ttaxonid\tdepth\tname\n");

		if (true || tidCount[0]) { 
		/// !!! must have this even if unclassiied is 0; otherwise Pavian will fail
			printf("%6.2f\t", tidCount[0] * 100.0 / totalReads);
			printf("%d\t", tidCount[0]);
			printf("%d\t", tidCount[0]);
			printf("U\t"); 
			printf("0\t");
			printf("0\t");
			printf("unclassified\n");
		}
		printTable_(1,0);
	}
};

int main(int argc, char **argv) {
	if (argc < 4) {
		fprintf(stderr, "Usage: %s <node.dmp> <names.dmp> <xxx.lsam.id> [scoreThreshold = 40]\n", argv[0]);
		exit(1);
	}

	int scoreT = 40;
	if (argc >= 5) {
		scoreT = atoi(argv[4]);
	}

	db.readNode(argv[1]);
	db.readName(argv[2]);

	gzFile fp = string(argv[3]) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(argv[3], "r");
	kxstream<gzFile> in(fp);
	string buf;

	TaxCounter taxCounter = {}; // need init

	while (in.get_until(kSepLine, buf) != kEOF) {
		// name, flag, score, seq, qual, taxLabels
		// fprintf(stderr, "%s\n", buf.c_str());
		vector<string> rec = splitBy(buf, '\t');
		int score = atoi(rec[2].c_str());
		if (score < scoreT) rec[5] = "*";

		taxCounter.add1read(rec[5]);
	}

	gzclose(fp);
	taxCounter.printTable();
	return 0;
}
