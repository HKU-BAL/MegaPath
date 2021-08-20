#include <map>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

typedef long long LL;

int main(int argc, char **argv) {
	if (argc > 1) {
		cerr << "Usage: " << argv[0] << " < in.m8" << endl;
		return 2;
	}


	map<string, vector<pair<LL, LL> > > intervals;
	string qid, sid;
	double identity;
	LL algn_len, mm, gaps, qs, qe, ss, se;
	double e_value, bit_score;

	while (cin >> qid >> sid
		       >> identity >> algn_len
		       >> mm >> gaps >> qs >> qe
		       >> ss >> se >> e_value >> bit_score) {
		if (ss > se) { swap(ss, se); }
		intervals[sid].push_back(make_pair(ss, se));
	}

	for (auto it = intervals.begin(); it != intervals.end(); ++it) {
		sort(it->second.begin(), it->second.end());
		auto intv = it->second[0];
		LL cov = 0;

		cout << it->first << '\t';
		for (int i = 1; i < it->second.size(); ++i) {
			if (it->second[i].first <= intv.second) {
				intv.second = max(it->second[i].second, intv.second);
			} else {
				cov += intv.second - intv.first + 1;
				cout << intv.first << "," << intv.second << ";";
				intv = it->second[i];
			}
		}

		cov += intv.second - intv.first + 1;
		cout << intv.first << "," << intv.second << ";";
		cout << '\t' << cov << endl;
		for (int i = 0; i < it->second.size(); ++i) {
			cout << it->second[i].first << ' ' << it->second[i].second << endl;
		}
	}

	return 0;
}