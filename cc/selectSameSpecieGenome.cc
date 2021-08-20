#include <iostream>
#include <map>
#include <vector>
#include <cstdlib>
#include <zlib.h>

#include "taxonomy.h"
#include "kxseq.h"

using namespace std;
using namespace kx;

int target_tids[] = { 337042,12461,11908,12092,333760,11292,11277,1714621,11269,12475,11709,11676,10298,1868658,122928,85106,10798,11161,11041,186538,10407,11103,162145,46839,46839,46839,46839,46839,46839,46839,46839,46839,46839,46839,46839,11620,290028,10245,335341,335341,335341,335341,335341,335341,335341,335341,37296,493803,28875,28875,28875,28875,28875,28875,28875,28875,28875,28875,28875,651580,11588,68887 };

TaxDB db;

int main(int argc, char **argv) {
	srand(10086);
	if (argc < 5) {
		cerr << "Usage: " << argv[0] << " <acc2tid> <nodes.dmp> <names.dmp> <nt.fa[.gz]>" << endl;
		return 1;
	}

	cerr << "Reading acc2tid" << endl;
	db.initAccToTax(argv[1]);
	cerr << "Reading nodes" << endl;
	db.readNode(argv[2]);
	cerr << "Reading names" << endl;
	db.readName(argv[3]);

	map<int, string> sp_seq;
	map<int, int> seen;
	for (int i = 0; i < sizeof(target_tids) / sizeof(int); ++i) {
		int id = db.popUpToSpecies(target_tids[i]);
		if (id != 0 && !sp_seq.count(id)) {
			cerr << "Added: " << id << db.taxInfo[id].name << endl;
			sp_seq[id] = "";
			seen[id] = 0;
		}
	}

	gzFile fp = string(argv[4]) == "" ? gzdopen(fileno(stdin), "r") : gzopen(argv[4], "r");
	kxseq<gzFile> ks(fp);

	while (ks.read() >= 0) {
		int tid = db.popUpToSpecies(db.acc2tid[getAccession(ks.name())]);
		if (sp_seq.count(tid) && ks.comment().find("complete genome") != string::npos) {
			cerr << ks.name() << ' ' << ks.comment() << ' ' << ks.seq().size() << endl;
			seen[tid]++;
			if (rand() % seen[tid] == 0) {
				sp_seq[tid] = ">" + ks.name() + " " + ks.comment() + "\n" + ks.seq();
			}
		}
	}

	for (auto it = sp_seq.begin(); it != sp_seq.end(); ++it) {
		cerr << it->first << " " << db.taxInfo[it->first].name << " seen by " << seen[it->first] << " times\n";
		if (seen[it->first] > 0) {
			cout << it->second << '\n';
		}
	}

	return 0;
}