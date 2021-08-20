#include <map>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <zlib.h>
#include "kxseq.h"
#include "misc.h"

using namespace std;
using namespace kx;

namespace option {
	vector<string> str_contids = {"9606", "32630"};
	double top = 0.94;
	double frac = 0.95;
}

map<int, int> counts;
map<pair<int, int>, int> contaminantOverlapCounts;

int main(int argc, char **argv) {
	int c;
	while ((c = getopt(argc, argv, "c:p:t:")) >= 0) {
		if (c == 'c') option::str_contids = splitBy(optarg, ',');
		else if (c == 't') option::top = atof(optarg);
		else if (c == 'f') option::frac = atof(optarg);
	}

	if (argc == optind) {
		cerr << "Usage: " << argv[0] << " [options] <in.lsam.id>" << endl;
		cerr << "options:" << endl;
		cerr << "    -c <INT[,INT...]>      tax id of contaminants [9606,32630]" << endl;
		cerr << "    -t <FLOAT>             ratio to determine a read is homologous to contaminants [0.94]" << endl;
		cerr << "    -f <FLOAT>             min fraction of reads to consider a species as contaminants [0.95]" << endl;
		return 1;
	}

	vector<int> contids;
	for (auto it = option::str_contids.begin(); it != option::str_contids.end(); ++it) {
		contids.push_back(atoi(it->c_str()));
		cerr << contids.back() << endl;
	}

	{ // 1-pass
		gzFile fp = string(argv[optind]) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(argv[optind], "r");
		kxstream<gzFile> ks(fp);
		string buf;

		while (ks.get_until(kSepLine, buf) != kEOF) {
			vector<string> rec = splitBy(buf, '\t');
	        assert(rec.size() >= 6);
	        vector<pair<string, double> > str_hits = splitAcc(rec[5]);
	        vector<pair<int, double> > hits;

	        for (auto it = str_hits.begin(); it != str_hits.end(); ++it) {
	        	hits.push_back(make_pair(atoi(it->first.c_str()), it->second));
	        }

	        map<int, double> contaminantScores;

	       	for (auto it = hits.begin(); it != hits.end(); ++it) {
	       		if (find(contids.begin(), contids.end(), it->first) != contids.end()) {
	       			contaminantScores[it->first] = it->second;
	       		}
	       	}

	       	for (auto it = hits.begin(); it != hits.end(); ++it) {
	   			for (size_t ci = 0; ci < contids.size(); ++ci) {
	   				if (it->first != contids[ci] &&
	   					it->second * option::top <= contaminantScores[contids[ci]] &&
	   					!contaminantScores.count(it->first)) {
	   					++contaminantOverlapCounts[make_pair(it->first, contids[ci])];
	   				}
	   			}

	   			++counts[it->first];
	       	}
		}

		gzclose(fp);
	}

	for (auto it = contaminantOverlapCounts.begin(); it != contaminantOverlapCounts.end(); ++it) {
		if (it->second > counts[it->first.first] * option::frac) {
			cerr << it->first.first << " is explained by " << it->first.second << ":"
				<< it->second
				<< "/" << counts[it->first.first] << endl;
		} else {
			cerr << it->first.first << " is not explained by " << it->first.second << ":"
				<< it->second
				<< "/" << counts[it->first.first] << endl;
		}
	}

	{ // 2-pass
		gzFile fp = string(argv[optind]) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(argv[optind], "r");
		kxstream<gzFile> ks(fp);
		string buf;

		while (ks.get_until(kSepLine, buf) != kEOF) {
			vector<string> rec = splitBy(buf, '\t');
	        assert(rec.size() >= 6);

	        vector<pair<string, double> > str_hits = splitAcc(rec[5]);
	        vector<pair<int, double> > hits;

	        for (auto it = str_hits.begin(); it != str_hits.end(); ++it) {
	        	hits.push_back(make_pair(atoi(it->first.c_str()), it->second));
	        }

	        map<int, double> contaminantScores;

	        for (auto it = hits.begin(); it != hits.end(); ++it) {
	       		if (find(contids.begin(), contids.end(), it->first) != contids.end()) {
	       			contaminantScores[it->first] = it->second;
	       		}
	       	}

	       	size_t numContaminants = contaminantScores.size();

	        vector<int> explainedBy;

   			for (size_t ci = 0; ci < contids.size(); ++ci) {
       			for (auto it = hits.begin(); it != hits.end(); ++it) {
	   				if (it->first != contids[ci] &&
	   					it->second * option::top <= contaminantScores[contids[ci]] &&
	   					contaminantOverlapCounts[make_pair(it->first, contids[ci])] > counts[it->first] * option::frac) {
	   					explainedBy.push_back(contids[ci]);
	   					break;
	   				}
	   			}
	       	}

	       	for (size_t i = 0; i < 5; ++i) {
	       		cout << rec[i] << '\t';
	       	}

	       	if (explainedBy.size() > 0) {
	       		for (size_t i = 0; i < explainedBy.size(); ++i) {
	       			if (i > 0) cout << ";";
	       			cout << contaminantScores[explainedBy[i]] << ',' << explainedBy[i];
	       		}
	       	} else if (hits.size() > 0) {
	       		for (size_t i = 0, first = 1; i < hits.size(); ++i) {
	       			if (hits.size() > numContaminants &&
	       				contaminantScores.count(hits[i].first)) 
	       			{
	       				continue;
	       			}

	       			if (!first) cout << ";";
	       			cout << hits[i].second << ',' << hits[i].first;
	       			first = 0;
	       		}
	       	} else {
	       		cout << "*";
	       	}

	       	for (size_t i = 6; i < rec.size(); ++i) {
	       		cout << '\t' << rec[i];
	       	}

	       	cout << '\n';
		}

		gzclose(fp);
	}

	return 0;
}