#include "kxseq.h"
#include <stdio.h>
#include <stdlib.h>

using namespace std;
using namespace kx;

static inline char Complement(char c) {
    if (c >= 0 && c < 4) {
        return 3 - c;
    }

    switch (c) {
    case 'A': {
        return 'T';
    }

    case 'C': {
        return 'G';
    }

    case 'G': {
        return 'C';
    }

    case 'T': {
        return 'A';
    }

    default: {
        return c;
    }
    }
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
	if (argc != 1) {
		fprintf(stderr, "Usage: %s < reads.fq > out.fq\n", argv[0]);
		return 1;
	}

    kxseq<int> seq(fileno(stdin));

    while (seq.read() >= 0) {
        string s = seq.seq();
        string q = seq.qual();
        reverse(s.begin(), s.end());
        reverse(q.begin(), q.end());
        for (int i = 0; i < s.size(); ++i) { s[i] = Complement(s[i]); }
        print_seq(seq.name(), seq.comment(), s, q);
    }

	return 0;
}