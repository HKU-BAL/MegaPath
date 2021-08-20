#include <unordered_map>
#include <string>
#include <fstream>
#include <iostream>
#include <zlib.h>
#include "kxseq.h"

using namespace std;
using namespace kx;

static const int kUndef = -100861008;

struct option {
	string map_file;
	int taxid;

	option() {
		map_file = "";
		taxid = kUndef;
	}
} opt;

unordered_map<string, int> seqid2taxid;

int main(int argc, char **argv) {
	int c;
	while ((c = getopt(argc, argv, "m:i:")) >= 0) {
		if (c == 'm') opt.map_file = optarg;
		else if (c == 'i') opt.taxid = atoi(optarg);
	}

	if (argc == optind || (opt.taxid == kUndef && opt.map_file == "") ||
		(opt.taxid != kUndef && opt.map_file != "")) {
		cerr << "Usage: " << argv[0] << " -m seqid2taxid.map ref.fa[.gz]" << endl;
		cerr << "   or: " << argv[0] << " -i taxid ref.fa[.gz]" << endl;
		return 1;
	}

	if (opt.map_file != "") {
		ifstream f_map(opt.map_file);
		string seqid;
		int taxid;
		while (f_map >> seqid >> taxid) {
			seqid2taxid[seqid] = taxid;
		}
	}

	gzFile fp = string(argv[optind]) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(argv[optind], "r");
	kxseq<gzFile> ks(fp);

	while (ks.read() >= 0) {
		if (opt.taxid >= 0) {
			cout << '>' << ks.name() << "|kraken:taxid|" << opt.taxid << '\n';
			cout << ks.seq() << '\n';
		} else {
			auto iter = seqid2taxid.find(ks.name());
			if (iter == seqid2taxid.end()) {
				cerr << "Error: cannot find taxid for " << ks.name() << endl;
			} else {
				cout << '>' << ks.name() << "|kraken:taxid|" << iter->second << '\n';
				cout << ks.seq() << '\n';
			}
		}
	}

	gzclose(fp);

	return 0;
}