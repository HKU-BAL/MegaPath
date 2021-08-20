#include <map>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <zlib.h>
#include "kxseq.h"

using namespace std;
using namespace kx;

typedef long long LL;

int getCoverage(vector<pair<LL, LL> > &v, int n) {
	if (n == 0) return 0;

	sort(v.begin(), v.begin() + n);
	auto intv = v[0];
	LL cov = 0;

	for (int i = 1; i < n; ++i) {
		if (v[i].first <= intv.second) {
			intv.second = max(v[i].second, intv.second);
		} else {
			cov += intv.second - intv.first + 1;
			intv = v[i];
		}
	}

	cov += intv.second - intv.first + 1;
	return cov;
}

int main(int argc, char **argv) {
	if (argc != 1 && argc != 3) {
		cerr << "Usage: " << argv[0] << " [ <ref.fa> <contig.fa> ] < in.m8" << endl;
		return 2;
	}

	map<string, int> QLen, TLen;
	bool calcAvg = false;

	if (argc >= 3) {
		calcAvg = true;

		gzFile fp = gzopen(argv[1], "r");
		kxseq<gzFile> seqT(fp);
    	while (seqT.read() >= 0) {
    		TLen[seqT.name()] = seqT.seq().length();
		}
		gzclose(fp);

		fp = gzopen(argv[2], "r");
		kxseq<gzFile> seqQ(fp);
    	while (seqQ.read() >= 0) {
    		QLen[seqQ.name()] = seqQ.seq().length();
		}
		gzclose(fp);
	}

	map<string, vector<pair<LL, LL> > > intervals;
	map<string, vector<pair<LL, int> > > queryAlgnmtLen; // first: len, second: index
	map<string, vector<string> > queryID;
	string qid, sid, lastq;
	double identity;
	LL algn_len, mm, gaps, qs, qe, ss, se;
	double e_value, bit_score;

	while (cin >> qid >> sid
		       >> identity >> algn_len
		       >> mm >> gaps >> qs >> qe
		       >> ss >> se >> e_value >> bit_score) {
		if (lastq != qid) {
			if (ss > se) { swap(ss, se); }
			intervals[sid].push_back(make_pair(ss, se));
			queryAlgnmtLen[sid].push_back(make_pair(abs(qe - qs) + 1, (int)queryAlgnmtLen[sid].size()));
			queryID[sid].push_back(qid);
			lastq = qid;
		}
	}

	for (auto it = queryAlgnmtLen.begin(); it != queryAlgnmtLen.end(); ++it) {
		sort(it->second.rbegin(), it->second.rend());
		cout << "Target: " << it->first;
		if (calcAvg) { cout << "\t" << TLen[it->first]; }
		cout << endl;

		vector<pair<LL, LL> > v;
		double totalMappingLen = 0, totalLength = 0;
		int NC50 = 0;
		for (int i = 0; i < it->second.size(); ++i) {
			v.push_back(intervals[it->first][it->second[i].second]);
			totalLength += QLen[queryID[it->first][it->second[i].second]];
			totalMappingLen += it->second[i].first;

			cout << it->second[i].first << '\t'
				<< QLen[queryID[it->first][it->second[i].second]] << '\t' 
				<< double(it->second[i].first) / QLen[queryID[it->first][it->second[i].second]] << '\t'
				<< getCoverage(v, i+1) << endl;

			if (NC50 == 0 && totalMappingLen >= .5 * TLen[it->first]) {
				NC50 = it->second[i].first;
			}
		}
		cout << "Mapping Ratio: " << totalMappingLen / totalLength
			<< "\tAvg Mapping Length: " << totalMappingLen / v.size()
			<< "\tNC50: " << NC50 << endl;
	}

	return 0;
}