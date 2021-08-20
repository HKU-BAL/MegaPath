#include "kxseq.h"
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <ctype.h>
#include <zlib.h>

using namespace std;
using namespace kx;

int main(int argc, char **argv) {
	if (argc < 2) {
		fprintf(stderr, "Usage: %s <in.fa> [min_mask_len=15]\n", argv[0]);
		return 1;
	}

	int min_mask_len = 15;
	if (argc >= 3) {
		min_mask_len = atoi(argv[2]);
	}

	gzFile fp = string(argv[1]) == "-" ? gzdopen(fileno(stdin), "r") : gzopen(argv[1], "r");
    kxseq<gzFile> seq(fp);

    while (seq.read() >= 0) {
    	string s = seq.seq();
    	for (size_t i = 0, j = 0; i < s.size(); i = j) {
    		j = i + 1;
    		if (!islower(s[i])) { continue; }

    		while (j < s.size() && islower(s[j])) {
    			++j;
    		}

    		if (j - i >= min_mask_len) {
    			for (size_t k = i; k < j; ++k) {
    				s[k] = 'N';
    			}
    		}
    	}
    	printf(">%s%s\n%s\n", seq.name().c_str(), seq.comment().size() ? (" " + seq.comment()).c_str(): s.c_str());
    }

    gzclose(fp);

	return 0;
}
