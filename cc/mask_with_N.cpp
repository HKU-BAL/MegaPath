#include "kxseq.h"
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <unordered_map>

using namespace std;
using namespace kx;

int main(int argc, char **argv) {
	if (argc < 2) {
		fprintf(stderr, "Usage: %s <mapped.region> < in.fa > out.fa\n", argv[0]);
		return 1;
	}

    kxseq<int> seq(fileno(stdin));
    unordered_map<string, string> ref;

    while (seq.read() >= 0) {
        ref[seq.name()] = seq.seq();
    }

    ifstream mr(argv[1]);
    string ref_id;
    int pos, len;

    while (mr >> ref_id >> pos >> len) {
        if (pos + len > ref[ref_id].length()) {
            pos = ref[ref_id].length() - len;
        }
        fill(ref[ref_id].begin() + pos, ref[ref_id].begin() + pos + len, 'N');
    }

    for (auto it = ref.begin(); it != ref.end(); ++it) {
        printf(">%s\n%s\n", it->first.c_str(), it->second.c_str());
    }

	return 0;
}