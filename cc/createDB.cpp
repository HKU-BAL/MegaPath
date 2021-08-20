#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <assert.h>
#include <zlib.h>
#include "kxseq.h"
#include "taxonomy.h"

using namespace std;
using namespace kx;

unordered_map<string, int> acc2tid;
vector<TaxInfo> taxInfos(2000000);

void printUsage(char *programName) {
	cerr << programName << ": create a FASTA file from NCBI-nt and NCBI-UniVec database." << endl;
	cerr << "- remove all artificial sequences in NCBI-nt" << endl;
	cerr << "- add UniVec to the DB" << endl;
	cerr << "- remove all sequences that cannot matched any record in NCBI Taxonmy database" << endl;
	cerr << "- reformat the FASTA header to >ACCESSION.VERSION,ACCESSION.VERSION..." << endl;
	cerr << endl;
	cerr << "Usage: " << programName << " <nt> <UniVec> <human_seq> <accession2taxid> <nodes.dmp> <names.dmp> [prefilter_mode=0] >out.fasta" << endl;
	cerr << "All files can be gzip'ed." << endl;
	cerr << "prefilter_mode will output accession2taxid that matches ACCESSION.VERSION in nt and UniVec to stdout" << endl;
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
		if (end != string::npos && header.substr(start, end - start) == "gi") {
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
			ret.push_back(removeVersion(header.substr(start, end - start)));
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

struct OutputAcc2Tid {
	int tid;
	bool output;

	OutputAcc2Tid(int tid, bool output): tid(tid), output(output) {}

	void operator() (const kxseq<gzFile> &ks) {
		if (output) {
			vector<string> acc = header2acc(ks.name() + " " + ks.comment());
			cout << acc[0] << '\t' << acc[0] << '\t' << tid << "\t*\n";
		}
	}
};

struct ReformatNT {
	static bool belongsTo(const string &acc, const string &taxName) {
		int tid = acc2tid[acc];
	    while (tid != 1 && tid != 0) {
	    	if (taxInfos[tid].name == taxName) return true;
	        tid = taxInfos[tid].parent;
	    }

	    return false;
	}

	void operator() (const kxseq<gzFile> &ks) {
		vector<string> acc = header2acc(ks.name() + " " + ks.comment());
		string newHeader;

		for (auto it = acc.begin(); it != acc.end(); ++it) {
			if (acc2tid[*it] != -1 && !belongsTo(*it, "artificial sequences")) {
				if (newHeader.length() > 0) newHeader.push_back(',');
				newHeader.append(*it);
			}
		}

		if (newHeader.length() > 0) {
			cout << ">" << newHeader << '\n';
			cout << ks.seq() << '\n';
		}
	}
};

struct ReformatUVHG {
	void operator() (const kxseq<gzFile> &ks) {
		vector<string> acc = header2acc(ks.name() + " " + ks.comment());
		string newHeader;

		for (auto it = acc.begin(); it != acc.end(); ++it) {
			if (newHeader.length() > 0) newHeader.push_back(',');
			newHeader.append(*it);
		}

		if (newHeader.length() > 0) {
			cout << ">" << newHeader << '\n';
			cout << ks.seq() << '\n';
		}
	}
};

template<typename NtFuncType, typename UvFuncType, typename HgFuncType>
void processNtAndUv(const string &ntFile, const string &uvFile, const string &hgFile, NtFuncType ntFunc, UvFuncType uvFunc, HgFuncType hgFunc) {
	gzFile fpNT = ntFile == "-" ? gzdopen(fileno(stdin), "r") : gzopen(ntFile.c_str(), "r");
	kxseq<gzFile> ksNT(fpNT);

	while (ksNT.read() >= 0) {
		ntFunc(ksNT);
	}

	gzclose(fpNT);

	gzFile fpUV = gzopen(uvFile.c_str(), "r");
	kxseq<gzFile> ksUV(fpUV);

	while (ksUV.read() >= 0) {
		uvFunc(ksUV);
	}

	gzclose(fpUV);

	gzFile fpHG = gzopen(hgFile.c_str(), "r");
	kxseq<gzFile> ksHG(fpHG);

	while (ksHG.read() >= 0) {
		hgFunc(ksHG);
	}

	gzclose(fpHG);
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
	if (argc < 7) {
		printUsage(argv[0]);
		return 1;
	}

	string ntFile = argv[1];
	string uvFile = argv[2];
	string hgFile = argv[3];
	string acc2tidFile = argv[4];
	string nodeFile = argv[5];
	string nameFile = argv[6];

	bool prefilter = argc > 7;

	// read nt and uv, keep headers in acc2tid
	cerr << "Reading NT and UV" << endl;
	processNtAndUv(ntFile, uvFile, hgFile, ParseHeader(), OutputAcc2Tid(32630, prefilter), OutputAcc2Tid(9606, prefilter));
	cerr << "Number of ACCESSIONs: " << acc2tid.size() << endl;

	// read accession2taxid
	cerr << "Reading accession2taxid" << endl;
	initAccToTax(acc2tidFile, prefilter);

	if (!prefilter) {
		// read nodes and names
		cerr << "Reading nodes and names" << endl;
		readNodeAndName(nodeFile, nameFile);

		// read nt and uv, keep non-artifical part of nt, and all uv, and output to stdout
		processNtAndUv(ntFile, uvFile, hgFile, ReformatNT(), ReformatUVHG(), ReformatUVHG());
	}

	cerr << "Unfound accessions (may be none):" << endl;
	for (auto it = acc2tid.begin(); it != acc2tid.end(); ++it) {
		if (it->second == -1) {
			cerr << it->first << endl;
		}
	}

	return 0;
}