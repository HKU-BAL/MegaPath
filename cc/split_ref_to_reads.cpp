#include "kxseq.h"
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

using namespace std;
using namespace kx;

int main(int argc, char **argv) {
	if (argc < 3) {
		fprintf(stderr, "Usage: %s <read_len> <overlap> < ref.fa > reads.fa\n", argv[0]);
		return 1;
	}

    kxseq<int> seq(fileno(stdin));
    int read_len = atoi(argv[1]);
    int overlap = atoi(argv[2]);

    while (seq.read() >= 0) {
    	for (unsigned i = 0; ; i += overlap) {
    		if (i + read_len > seq.seq().length()) {
    			i = max(0, (int)seq.seq().length() - read_len);
    		}
    		printf(">%s_%d\n%s\n", seq.name().c_str(), i, seq.seq().substr(i, min(read_len, (int)seq.seq().length())).c_str());
    		if (i + read_len >= seq.seq().length()) break;
    	}
    }
    
	return 0;
}