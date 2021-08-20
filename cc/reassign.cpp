#include <map>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <zlib.h>
#include <omp.h>
#include <unordered_map>
#include <unordered_set>
#include "kxseq.h"
#include "misc.h"

using namespace std;
using namespace kx;

namespace option {
	double u = 20;
	double v = 0.05;
	double t = 40;
	int n_threads = 1;
	bool output_seq = false;
}

unordered_map<int, int> counts;
unordered_map<int, int> uniq_counts;
unordered_map<uint64_t, int> intersect;
unordered_set<uint64_t> explains;
unordered_set<int> weakly_explained;

uint64_t pairup(int tax1, int tax2) {
	if (tax1 < tax2) swap(tax1, tax2);
	return (tax1 * 1ULL << 32) | tax2;
}

bool weakly_explain(int tax1, int tax2) {
	if (uniq_counts[tax1] <= option::u * uniq_counts[tax2]) return false;
	if (counts[tax1] - intersect[pairup(tax1, tax2)] <= option::v * counts[tax1]) return false;
	return true;
}

int main(int argc, char **argv) {
	int c;
	while ((c = getopt(argc, argv, "u:v:p:t:s")) >= 0) {
		if (c == 'u') option::u = atof(optarg);
		else if (c == 'v') option::v = atof(optarg);
		else if (c == 't') option::t = atof(optarg);
		else if (c == 'p') option::n_threads = atoi(optarg);
		else if (c == 's') option::output_seq = true;
	}

	if (argc == optind) {
		cerr << "Usage: " << argv[0] << " [options] <in.lsam.id>" << endl;
		cerr << "options:" << endl;
		cerr << "    -u <FLOAT>      [20]" << endl;
		cerr << "    -v <FLOAT>      [0.05]" << endl;
		cerr << "    -t <INT>        Alignment threshold [40]" << endl;
		cerr << "    -p <INT>        Number of CPU threads [1]" << endl;
		cerr << "    -s              Output sequences and Q-Values" << endl;
		return 1;
	}

	omp_set_num_threads(option::n_threads);

	{ // 1-pass

		gzFile fp = string(argv[optind]) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(argv[optind], "r");
		kxstream<gzFile> ks(fp);
		vector<string> buf(1 << 20);
		size_t n = 0;

		while (true) {
			while (n < buf.size() && ks.get_until(kSepLine, buf[n]) != kEOF) {
				++n;
			}

			vector<unordered_map<int, int> > counts(option::n_threads);
			vector<unordered_map<int, int> > uniq_counts(option::n_threads);
			vector<unordered_map<uint64_t, int> > intersect(option::n_threads);

#pragma omp parallel for schedule(dynamic)
			for (size_t i = 0; i < n; ++i) {
				vector<string> rec = splitBy(buf[i], '\t');
		        assert(rec.size() >= 6);
		        vector<pair<string, double> > str_hits = splitAcc(rec[5]);
		        vector<pair<int, double> > hits;

		        if (atof(rec[2].c_str()) < option::t || str_hits.size() == 0) {
		        	continue;
		        }

		        int thread_id = omp_get_thread_num();

		        for (auto it = str_hits.begin(); it != str_hits.end(); ++it) {
		        	int taxid = atoi(it->first.c_str());
		        	counts[thread_id][taxid]++;
		        	if (str_hits.size() == 1) {
		        		uniq_counts[thread_id][taxid]++;
		        	} else {
		        		for (auto jt = hits.begin(); jt != hits.end(); ++jt) {
		        			intersect[thread_id][pairup(jt->first, taxid)]++;
		        		}
		        	}

		        	hits.push_back(make_pair(taxid, it->second));
		        }
		    }

			// join all results
			for (int ti = 0; ti < (option::n_threads); ++ti) {
				for (auto it = counts[ti].begin(); it != counts[ti].end(); ++it) {
					::counts[it->first] += it->second;
				}


				for (auto it = uniq_counts[ti].begin(); it != uniq_counts[ti].end(); ++it) {
					::uniq_counts[it->first] += it->second;
				}


				for (auto it = intersect[ti].begin(); it != intersect[ti].end(); ++it) {
					::intersect[it->first] += it->second;
				}
			}

			if (n == 0) { break; }
			n = 0;
		}

		// calculate explain

		for (auto it = intersect.begin(); it != intersect.end(); ++it) {
			int t1 = it->first >> 32, t2 = it->first & 0xFFFFFFFF;
			if (weakly_explain(t1, t2)) {
				weakly_explained.insert(t2);
			} else if (weakly_explain(t2, t1)) {
				weakly_explained.insert(t1);
			}
		}

		for (auto it = intersect.begin(); it != intersect.end(); ++it) {
			int t1 = it->first >> 32, t2 = it->first & 0xFFFFFFFF;
			if (weakly_explain(t1, t2)) {
				if (!weakly_explained.count(t1)) {
					explains.insert((t1 * 1ULL << 32) | t2);
					cerr << t1 << " explains " << t2 << "\n";
				}
			} else if (weakly_explain(t2, t1)) {
				if (!weakly_explained.count(t2)) {
					explains.insert((t2 * 1ULL << 32) | t1);
					cerr << t2 << " explains " << t1 << "\n";
				}

			}
		}

		gzclose(fp);
	}


	{ // 2-pass
		gzFile fp = string(argv[optind]) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(argv[optind], "r");
		kxstream<gzFile> ks(fp);
		vector<string> buf(1 << 20);
		size_t n = 0;

		while (true) {
			while (n < buf.size() && ks.get_until(kSepLine, buf[n]) != kEOF) {
				++n;
			}

#pragma omp parallel for schedule(dynamic) // for pairing up
			for (size_t i = 0; i < n; i++) {
				vector<string> rec = splitBy(buf[i], '\t');

		        assert(rec.size() >= 6);
		        vector<pair<string, double> > str_hits = splitAcc(rec[5]);
		        vector<pair<int, double> > hits;

		        for (auto it = str_hits.begin(); it != str_hits.end(); ++it) {
		        	hits.push_back(make_pair(atoi(it->first.c_str()), it->second));
		        }

		        if (!option::output_seq) {
		        	rec[3] = rec[4] = "*";
		        }

		        buf[i] = rec[0] + '\t' + rec[1] + '\t' + rec[2] + '\t' + rec[3] + '\t' + rec[4] + '\t';
		        bool first_taxid = true;

		        for (auto it = hits.begin(); it != hits.end(); ++it) {
		        	bool remain = true;
		        	for (auto jt = hits.begin(); jt != hits.end(); ++jt) {
		        		if (explains.count((jt->first * 1ULL << 32) | it->first)) {
		        			remain = false;
		        			break;
		        		}
		        	}
		        	if (remain) {
		        		if (!first_taxid) buf[i] += ";";
		        		else first_taxid = false;
		        		buf[i] += to_string((long long)(it->second)) + "," + to_string((long long)(it->first));
		        	}
		        }

		        if (hits.size() == 0) {
		        	buf[i] += "*";
		        }

		        for (size_t j = 6; j < rec.size(); ++j) {
		        	buf[i] += '\t' + rec[j];
		        }
		        buf[i] += "\n";
			}

			for (size_t i = 0; i < n; i++) {
				fputs(buf[i].c_str(), stdout);
			}

			if (n == 0) { break; }
			n = 0;
		}

		gzclose(fp);
	}

	return 0;
}