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
TaxDB db;

struct TaxCounter {

	struct Counter {
		int uniqCount;
		int nonUniqCount;
		Counter() { uniqCount = 0; nonUniqCount = 0; }
	};

	unordered_map<int, Counter> tidCount;
	unordered_map<int, set<int> > sons;

	void createLineage_(int sp, int g, int f, int sk) {
		sons[g].insert(sp);
		sons[f].insert(g);
		sons[sk].insert(f);
		sons[kRoot].insert(sk);
	}

	void addToTidCount_(const set<int> &st) {
		if (st.size() == 1) {
			tidCount[*st.begin()].uniqCount++;
		} else {
			for (auto it = st.begin(); it != st.end(); ++it) {
				tidCount[*it].nonUniqCount++;
			}
		}
	}

	void add1read(const string &s) {
		if (s.length() == 0 || s == "*") return;
		vector<string> tid_and_score = splitBy(s, ';');
		vector<string> tids;
		for (unsigned i = 0; i < tid_and_score.size(); ++i) {
			auto s_tid = splitBy(tid_and_score[i], ',');
			tids.push_back(s_tid[1]);
		}

		set<int> st_sp, st_g, st_f, st_sk;

		for (size_t i = 0; i < tids.size(); ++i) {
			int sp = -1, g = -1, f = -1, sk = -1;
			int tid = atoi(tids[i].c_str());
			while (tid != 1 && tid != 0) {
				if (db.taxInfo[tid].rank == "species") {
					sp = tid;
				} else if (db.taxInfo[tid].rank == "genus") {
					g = tid;
				} else if (db.taxInfo[tid].rank == "family") {
					f = tid;
				} else if (db.taxInfo[tid].rank == "superkingdom") {
					sk = tid;
				}
				tid = db.taxInfo[tid].parent;
			}

			if (sp == -1) { continue; }
			else { st_sp.insert(sp); }

			g = g == -1 ? (sp | 0x80000000) : g;
			f = f == -1 ? (g | 0x40000000) : f;

			st_g.insert(g);
			st_f.insert(f);
			st_sk.insert(sk);

			createLineage_(sp, g, f, sk);
		}

		addToTidCount_(st_sp);
		addToTidCount_(st_g);
		addToTidCount_(st_f);
		addToTidCount_(st_sk);
	}

	struct cmp_ {
		unordered_map<int, Counter> *counter;
		cmp_(unordered_map<int, Counter> *counter) : counter(counter) {}
		bool operator() (int a, int b) {
			return (*counter)[a].uniqCount > (*counter)[b].uniqCount;
		}
	};

	void printTable_(int tid) {
		if (tid >= 0 && (tid & 0xC0000000) == 0) {
			printf("%s", db.taxInfo[tid].rank.c_str());

			// get lineage
			int tmpTid = tid;
			vector<string> lineage(4, "-");

			while (tmpTid != 1 && tmpTid != 0) {
				if (db.taxInfo[tmpTid].rank == "species") {
					lineage[3] = db.taxInfo[tmpTid].name;
				} else if (db.taxInfo[tmpTid].rank == "genus") {
					lineage[2] = db.taxInfo[tmpTid].name;
				} else if (db.taxInfo[tmpTid].rank == "family") {
					lineage[1] = db.taxInfo[tmpTid].name;
				} else if (db.taxInfo[tmpTid].rank == "superkingdom") {
					lineage[0] = db.taxInfo[tmpTid].name;
				}
				tmpTid = db.taxInfo[tmpTid].parent;
			}

			for (size_t i = 0; i < lineage.size(); ++i) {
				printf("\t%s", lineage[i].c_str());
			}

			printf("\t%d\t%d\n", tidCount[tid].uniqCount, tidCount[tid].nonUniqCount);
		}

		vector<int> mySons(sons[tid].begin(), sons[tid].end());
		sort(mySons.begin(), mySons.end(), cmp_(&tidCount));

		for (auto it = mySons.begin(); it != mySons.end(); ++it) {
			printTable_(*it);
		}
	}

	void printTable() {
		printTable_(kRoot);
	}
};

void printLabel(const vector<string> &rec, FILE *fpOutLabel) {
	for (int i = 0; i < 5; ++i) {
		fprintf(fpOutLabel, "%s\t", rec[i].c_str());
	}

	if (rec[5].length() == 0 || rec[5] == "*") return;
	vector<string> tids = splitBy(rec[5], ';');
	set<int> st_sp;

	for (size_t i = 0; i < tids.size(); ++i) {
		int tid = atoi(tids[i].c_str());
		while (tid != 1 && tid != 0) {
			if (db.taxInfo[tid].rank == "species") {
				st_sp.insert(tid);
				break;
			}
			tid = db.taxInfo[tid].parent;
		}
	}

	for (auto it = st_sp.begin(); it != st_sp.end(); ++it) {
		if (it == st_sp.begin())
			fprintf(fpOutLabel, "%s", db.getTaxInfo(*it).c_str());
		else
			fprintf(fpOutLabel, ";%s", db.getTaxInfo(*it).c_str());
	}
	fprintf(fpOutLabel, "\n");
}

int main(int argc, char **argv) {
	if (argc < 4) {
		fprintf(stderr, "Usage: %s <node.dmp> <names.dmp> <xxx.lsam.id> [scoreThreshold = 40] [out.label]\n", argv[0]);
		exit(1);
	}

	int scoreT = 40;
	if (argc >= 5) {
		scoreT = atoi(argv[4]);
	}

	FILE *fpOutLabel = NULL;
	if (argc >= 6) {
		fpOutLabel = fopen(argv[5], "w");
	}

	db.readNode(argv[1]);
	db.readName(argv[2]);

	gzFile fp = string(argv[3]) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(argv[3], "r");
	kxstream<gzFile> in(fp);
	string buf;

	TaxCounter taxCounter;

	while (in.get_until(kSepLine, buf) != kEOF) {
		// name, flag, score, seq, qual, taxLabels
		// fprintf(stderr, "%s\n", buf.c_str());
		vector<string> rec = splitBy(buf, '\t');
		int score = atoi(rec[2].c_str());

		if (score >= scoreT) {
			taxCounter.add1read(rec[5]);
		}

		if (fpOutLabel) {
			printLabel(rec, fpOutLabel);
		}
	}

	if (fpOutLabel) fclose(fpOutLabel);
	gzclose(fp);
	taxCounter.printTable();

	return 0;
}