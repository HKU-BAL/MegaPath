#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string>
#include "kxseq.h"

using namespace kx;
using namespace std;

static inline void trim_readno(string &s) {
    if (s.length() > 2 && s[s.length()-2] == '/' && isdigit(s[s.length()-1]))
        s.resize(s.length() - 2);
}

string name1, comm1, seq1, qual1;

inline void update1(kxseq<int> &seq) {
    name1 = seq.name();
    comm1 = seq.comment();
    seq1 = seq.seq();
    qual1 = seq.qual();
}

void print_seq(const string &name, const string &comm, const string &seq, const string &qual, int whichEnd, FILE* fp = stdout) {
	fprintf(fp, "%c", qual.length() ? '@' : '>');
	fprintf(fp, "%s", name.c_str());
    if (whichEnd > 0) {
        fprintf(fp, "/%d", whichEnd);
    }
	if (comm.length()) {
		fprintf(fp, " %s", comm.c_str());
	}
	fprintf(fp, "\n%s\n", seq.c_str());
	if (qual.length()) {
		fprintf(fp, "+\n%s\n", qual.c_str());
	}
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "cat a.fq | %s <output_prefix>\n", argv[0]);
        exit(1);
    }

    FILE *out_se = fopen((string(argv[1]) + ".se.fq").c_str(), "w");
    FILE *out_pe1 = fopen((string(argv[1]) + ".pe_1.fq").c_str(), "w");
    FILE *out_pe2 = fopen((string(argv[1]) + ".pe_2.fq").c_str(), "w");

    kxseq<int> seq(fileno(stdin));
    bool has_last = false;

    while (seq.read() >= 0) {
        trim_readno(seq.name());
        if (has_last) {
             if (name1 == seq.name()) {
             	print_seq(name1, comm1, seq1, qual1, 1, out_pe1);
                print_seq(seq.name(), seq.comment(), seq.seq(), seq.qual(), 2, out_pe2);
                has_last = false;
             } else {
                print_seq(name1, comm1, seq1, qual1, 0, out_se);
                update1(seq);
             }
        } else {
            has_last = true;
            update1(seq);
        }
    }

    if (has_last) {
        print_seq(name1, comm1, seq1, qual1, 0, out_se);
    }

    fclose(out_se);
    fclose(out_pe1);
    fclose(out_pe2);

    return 0;
}