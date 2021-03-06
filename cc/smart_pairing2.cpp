#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string>
#include "kmseq.h"

using namespace kmlib;

static inline void trim_readno(kmstring &s) {
    if (s.length() > 2 && s[s.length()-2] == '/' && isdigit(s[s.length()-1]))
        s.resize(s.length() - 2);
}

string name1, comm1, seq1, qual1;

inline void update1(kmseq<gzFile> &seq) {
    name1 = seq.name();
    comm1 = seq.comment();
    seq1 = seq.seq();
    qual1 = seq.qual();
}

void print_seq(const string &name, const string &comm, const string &seq, const string &qual, FILE* fp = stdout) {
	fprintf(fp, "%c", qual.length() ? '@' : '>');
	fprintf(fp, "%s", name.c_str());
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
        fprintf(stderr, "cat a.fq | %s <SE_file_name>\n", argv[0]);
        exit(1);
    }

    FILE *out_se = fopen(argv[1], "w");
    kmseq<int> seq(fileno(stdin));
    bool has_last = false;

    while (seq.read() >= 0) {
        trim_readno(seq.name());
        if (has_last) {
             if (name1 == seq.name()) {
             	print_seq(name1, comm1, seq1, qual1);
                print_seq(seq.name(), seq.comment(), seq.seq(), seq.qual());
                has_last = false;
             } else {
                print_seq(name1, comm1, seq1, qual1, out_se);
                update1(seq);
             }
        } else {
            has_last = true;
            update1(seq);
        }
    }

    if (has_last) {
        print_seq(name1, comm1, seq1, qual1, out_se);
    }

    return 0;
}