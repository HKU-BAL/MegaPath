#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "kxseq.h"

using namespace std;
using namespace kx;

int A = 1, B = 2, O = 3, E = 1;
double d = 0.95;
bool lsam = false;

int calc_score(const char* cigar, int ed) {
	char *next_;
	int mp_len = 0;
	int indel = 0;
	int indel_p = 0;
	while (isdigit(*cigar)) {
		int len = strtol(cigar, &next_, 10);
		if (*next_ == 'M') {
			mp_len += len;
		} else if (*next_ == 'I' || *next_ == 'D') {
			indel += len;
			indel_p += O + E * len;
		}
		cigar = next_ + 1;
	}

	return mp_len * A - (ed - indel) * (A + B) - indel_p;
}

int main(int argc, char **argv) {
	char c;
	while ((c = getopt(argc, argv, "A:B:O:E:d:l")) >= 0) {
		if (c == 'A') A = atoi(optarg);
		else if (c == 'B') B = atoi(optarg);
		else if (c == 'O') O = atoi(optarg);
		else if (c == 'E') E = atoi(optarg);
		else if (c == 'd') d = atof(optarg);
		else if (c == 'l') lsam = true;
		else return 1;
	}

	char comp[256];
	for (int i = 0; i < 10; ++i) {
		comp["ACGTacgtNn"[i]] = "TGCAtgcaNn"[i];
	}

	if (optind != argc) {
		fprintf(stderr, "Usage %s [options] <in.sam >out.cfq\n", argv[0]);
		fprintf(stderr, "Options:\n"
			            "	-A <int>    match score [1]\n"
			            "	-B <int>    mismatch penalty [4]\n"
						"	-O <int>    gap open penalty [6]\n"
						"	-E <int>    gap extension penalty [1]\n"
						"	-d <float>  dropout ratio [1]\n"
						"	-l          output lsam\n");
		return 1;
	}

	string query, target, pos, seq, qual, buf;
	const char* cigar = NULL;
	int flag = 0;

	kxstream<int> in(fileno(stdin));

	while (in.get_until(kSepLine, buf) != kEOF) {
		if (buf[0] == '@') continue;
		size_t pos = 0;
		size_t next = 0;
		int ed = 0;
		bool NM_done = false, XA_done = false, XC_done = false;
		int max_score = 0;
		bool secondary = false;
		vector<pair<string, int> > v_aln;

		for (int i = 0; pos < buf.length(); ++i, pos = next + 1) {
			next = buf.find_first_of('\t', pos);
			if (next == string::npos) next = buf.length();
			if (i == 0) query = buf.substr(pos, next - pos);
			else if (i == 1) {
				flag = atoi(buf.c_str() + pos);
				if (flag & 0x100) {
					secondary = true;
					break;
				}
			}
			else if (i == 2) {
				target = buf.substr(pos, next - pos);
				size_t taxid_pos;
				if ((taxid_pos = target.find("kraken:taxid|")) != string::npos) {
					target = target.substr(taxid_pos + 13);
				}
			}
			else if (i == 5) cigar = buf.c_str() + pos;
			else if (i == 9) seq = buf.substr(pos, next - pos);
			else if (i == 10) qual = buf.substr(pos, next - pos);

			else if (!NM_done && buf.substr(pos, 5) == "NM:i:") {
				ed = atoi(buf.c_str() + pos + 5);
				NM_done = true;
			} else if (!XA_done && buf.substr(pos, 5) == "XA:Z:") {
				pos += 5;

				while (pos < buf.length() && !isspace(buf[pos])) {
					size_t next_ = buf.find_first_of(',', pos);
					string t = buf.substr(pos, next_ - pos);
					size_t taxid_pos;
					if ((taxid_pos = t.find("kraken:taxid|")) != string::npos) {
						t = t.substr(taxid_pos + 13);
					}
					pos = next_ + 1;

					pos = buf.find_first_of(',', pos) + 1;
					const char* cigar_ = buf.c_str() + pos;

					pos = buf.find_first_of(',', pos) + 1;
					char* next_p_;
					int ed_ = strtol(buf.c_str() + pos, &next_p_, 10);

					int score = calc_score(cigar_, ed_);
					if (score >= max_score * d - 1e-6) {
						v_aln.push_back(make_pair(t, score));
						if (score > max_score) max_score = score;
					}

					pos = next_p_ - buf.c_str() + 1;
				}
				XA_done = true;
			} else if (!XC_done && buf.substr(pos, 5) == "XC:Z:") {
				pos += 5;
				while (pos < buf.length() && !isspace(buf[pos])) {
					size_t next_ = buf.find_first_of(',', pos);
					string t = buf.substr(pos, next_ - pos);
					pos = next_ + 1;

					char* next_p_;
					int score = strtol(buf.c_str() + pos, &next_p_, 10);
					if (score >= max_score * d - 1e-6) {
						v_aln.push_back(make_pair(t, score));
						if (score > max_score) max_score = score;
					}

					pos = next_p_ - buf.c_str() + 1;
				}
				XC_done = true;
			}

			if (XC_done && XA_done && NM_done) {
				break;
			}
		}

		if (secondary) continue;

		if (!(flag & 0x4)) {
			int score = calc_score(cigar, ed);
			if (score >= max_score * d - 1e-6) { 
				v_aln.push_back(make_pair(target, score));
				if (score > max_score) max_score = score;
			}
		}

		bool is_fq = qual.length() == seq.length();

		if (flag & 0x10) {
			reverse(qual.begin(), qual.end());
			reverse(seq.begin(), seq.end());
			for (unsigned i = 0; i < seq.length(); ++i) {
				seq[i] = comp[seq[i]];
			}
		}

		if (lsam) {
			fputs(query.c_str(), stdout);
			printf("\t%d\t%d\t", flag, max_score);
			fputs(seq.c_str(), stdout);
			putchar('\t');
			fputs(qual.c_str(), stdout);
			putchar('\t');

			bool first = true;
			for (unsigned i = 0; i < v_aln.size(); ++i) {
				if (v_aln[i].second >= max_score * d - 1e-6) {
					if (!first) putchar(';');
					printf("%d,", v_aln[i].second);
					fputs(v_aln[i].first.c_str(), stdout);
					first = false;
				}
			}
			puts(first ? "*" : "");
		} else {
			putchar(is_fq ? '@' : '>');
			fputs(query.c_str(), stdout);

			if (!v_aln.empty()) {
				fputs("\tXC:Z:", stdout);
				for (unsigned i = 0; i < v_aln.size(); ++i) {
					if (v_aln[i].second >= max_score * d - 1e-6) {
						fputs(v_aln[i].first.c_str(), stdout);
						printf(",%d;", v_aln[i].second);
					}
				}
			}
			puts("");

			fputs(seq.c_str(), stdout);
			putchar('\n');
			if (is_fq) {
				fputs("+\n", stdout);
				fputs(qual.c_str(), stdout);
				puts("");
			}
		}
	}

	return 0;
}