#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <assert.h>
#include <zlib.h>
#include "kxseq.h"
#include "taxonomy.h"

using namespace std;
using namespace kx;

unordered_map<string, int> acc2tid;
vector<TaxInfo> taxInfos(2000000);
unordered_set<string> taxNamesToFilter;

void printUsage(char *programName) {
	cerr << programName << ": filter out sequences belongs to specified taxonomy" << endl;
	cerr << "Usage: " << programName << " <nt.fasta> <accession2taxid> <nodes.dmp> <names.dmp> [name1] [name2] ... >out.fasta" << endl;
	cerr << "All files can be gzip'ed." << endl;
}

vector<string> header2acc(const string &header) {
	vector<string> ret;

	if (header.length() > 7 && header.substr(0, 7) == "gnl|uv|") {
		size_t end = header.find_first_of(":");
		ret.push_back(removeVersion(header.substr(7, end - 7)));
		return ret;
	}

	size_t start = 0, end = header.find_first_of('|');
	while (start != string::npos) {
		if (end != string::npos) {
			assert(header.substr(start, end - start) == "gi");

			// old header
			size_t spliter = header.find_first_of('|', end + 1);
			size_t spliter2 = header.find_first_of('|', spliter + 1);
			size_t spliter3 = header.find_first_of('|', spliter2 + 1);
			ret.push_back(removeVersion(header.substr(spliter2 + 1, spliter3 - spliter2 - 1)));
			start = header.find_first_of('\1', spliter3 + 1);
		} else {
			// new header
			end = start;
			while (end < header.size() && !isspace(header[end]) && header[end] != '\1' && header[end] != '|') {
				++end;
			}
			ret.push_back(removeVersion(header.substr(start, end - start + 1)));
			if (end == header.size()) {
				start = string::npos;
			} else {
				start = header.find_first_of('\1', end);
			}
		}

		if (start == string::npos) {
			return ret;
		} else {
			start++;
			end = header.find_first_of('|', start);
		}
	}

	// cannot come here
	return ret;
}

struct ParseHeader {
	void operator() (const kxseq<gzFile> &ks) {
		vector<string> acc = header2acc(ks.name() + " " + ks.comment());
		for (auto it = acc.begin(); it != acc.end(); ++it) {
			acc2tid[*it] = -1;
		}
	}
};

struct Filter {
	static bool belongsTo(const string &acc, const unordered_set<string> &taxNames) {
		int tid = acc2tid[acc];
	    while (tid > 1) {
	    	if (taxNames.count(taxInfos[tid].name)) return true;
	        tid = taxInfos[tid].parent;
	    }

	    return false;
	}

	void operator() (const kxseq<gzFile> &ks) {
		vector<string> acc = header2acc(ks.name() + " " + ks.comment());
		bool keep = true;

		for (auto it = acc.begin(); it != acc.end(); ++it) {
			if (belongsTo(*it, taxNamesToFilter)) {
				keep = false;
				cerr << "Discarding " << ks.name() << ' ' << ks.comment() << endl;
				break;
			}
		}

		if (keep) {
			cout << ">" << ks.name() << ' ' << ks.comment() << '\n';
			cout << ks.seq() << '\n';
		}
	}
};

template<typename NtFuncType>
void processNT(const string &ntFile, NtFuncType ntFunc) {
	gzFile fpNT = ntFile == "-" ? gzdopen(fileno(stdin), "r") : gzopen(ntFile.c_str(), "r");
	kxseq<gzFile> ksNT(fpNT);

	while (ksNT.read() >= 0) {
		ntFunc(ksNT);
	}

	gzclose(fpNT);
}

void initAccToTax(const string &acc2tidFile, bool prefilter = false) {
	gzFile fp = acc2tidFile == "-" ? gzdopen(fileno(stdin), "r") : gzopen(acc2tidFile.c_str(), "r");
	kxstream<gzFile> ks(fp);

	string buf, acc, tid;

	while (ks.get_until(kSepSpace, buf) != kEOF) {
		ks.get_until(kSepSpace, acc);
		ks.get_until(kSepSpace, tid);
		acc = removeVersion(acc);

		if (acc2tid.count(acc)) {
			acc2tid[acc] = atoi(tid.c_str());

			if (prefilter) {
				cout << buf << '\t' << acc << '\t' << tid << "\t*\n";
			}
		}
		ks.get_until(kSepLine, buf);
	}
	gzclose(fp);
}

void readNodeAndName(const string &nodeFile, const string &nameFile) {
	gzFile fpNode = nodeFile == "-" ? gzdopen(fileno(stdin), "r") : gzopen(nodeFile.c_str(), "r");
	gzFile fpName = nameFile == "-" ? gzdopen(fileno(stdin), "r") : gzopen(nameFile.c_str(), "r");

	kxstream<gzFile> ksNode(fpNode);
	kxstream<gzFile> ksName(fpName);

	string tid, buf, parent, rank;

    while (ksNode.get_until(kSepSpace, tid) != kEOF) {
    	ksNode.get_until(kSepSpace, buf);
    	ksNode.get_until(kSepSpace, parent);
    	ksNode.get_until(kSepSpace, buf);
    	ksNode.get_until(kSepSpace, rank);

    	int id = atoi(tid.c_str());
    	if (id >= (int)taxInfos.size()) {
    		taxInfos.resize((id + 1) * 1.5);
    	}

    	taxInfos[id].parent = atoi(parent.c_str());
    	taxInfos[id].rank = rank;

    	ksNode.get_until(kSepLine, buf);
    }

    while (ksName.get_until(kSepLine, buf) != kEOF) {
    	if (buf.find("scientific name") == string::npos) {
    		continue;
    	}

    	istringstream is(buf);
    	int id;
    	char d;
    	is >> id >> d;
    	string name;
    	string nextpart;

    	is >> nextpart;

    	while (nextpart != "|") {
    		if (name.length() > 0) name += "_";
    		name += nextpart;
    		is >> nextpart;
    	}

    	taxInfos[id].name = name;
    }

    gzclose(fpNode);
    gzclose(fpName);
}

int main(int argc, char **argv) {
	if (argc < 5) {
		printUsage(argv[0]);
		return 1;
	}

	string ntFile = argv[1];
	string acc2tidFile = argv[2];
	string nodeFile = argv[3];
	string nameFile = argv[4];

	for (int i = 5; i < argc; ++i) {
		string s = argv[i];
		for (size_t j = 0; j < s.size(); ++j) {
			if (isspace(s[j])) s[j] = '_';
		}
		taxNamesToFilter.insert(s);
	}

	cerr << "Names to be filtered:" << endl;
	for (auto it = taxNamesToFilter.begin(); it != taxNamesToFilter.end(); ++it) {
		cerr << *it << endl;
	}

	// read nt and uv, keep headers in acc2tid
	cerr << "Reading NT" << endl;
	processNT(ntFile, ParseHeader());
	cerr << "Number of ACCESSIONs: " << acc2tid.size() << endl;

	// read accession2taxid
	cerr << "Reading accession2taxid" << endl;
	initAccToTax(acc2tidFile);

	// read nodes and names
	cerr << "Reading nodes and names" << endl;
	readNodeAndName(nodeFile, nameFile);

	// read nt and uv, keep non-artifical part of nt, and all uv, and output to stdout
	processNT(ntFile, Filter());

	return 0;
}