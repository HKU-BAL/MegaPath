#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string>
#include "kxseq.h"
#include "misc.h"

using namespace kx;
using namespace std;

const int kScoreIgnore = -1;
bool outputSeq = true;

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

void print_lsam_line(const string &name, const string &comm, const string &seq, const string &qual, int whichEnd, FILE* fp = stdout) {
    fputs(name.c_str(), fp);
    fputs("\t", fp);
    if (whichEnd == 1) {
        fprintf(fp, "%d", 0x40);
    } else if (whichEnd == 2) {
        fprintf(fp, "%d", 0x80);
    } else {
        fprintf(fp, "%d", 0x0);
    }
    fputs("\t", fp);

    int score = (comm == "IGNORE") ? kScoreIgnore :  atoi(comm.c_str() + 6);
    fprintf(fp, "%d\t", score);

    if (outputSeq) {
        fputs(seq.c_str(), fp);
        fputs("\t", fp);
        fputs(qual.c_str(), fp);
        fputs("\t", fp);
    } else {
        fputs("*\t*\t", fp);
    }

    if (score <= 0) {
        fputs("*", fp);
    } else {
        vector<string> score_acc = splitBy(comm, ';');
        bool isFirst = true;
        for (size_t i = 1; i < score_acc.size(); ++i) {
            // the first is score, so skip 0
            vector<string> sub = splitBy(score_acc[i], ',');
            for (size_t j = 1; j < sub.size(); ++j) {
                if (!isFirst) {
                    fputs(";", fp);
                } else {
                    isFirst = false;
                }

                fputs(sub[0].c_str(), fp);
                fputs(",", fp);
                fputs(sub[j].c_str(), fp);
            }
        }
    }

    if (score == kScoreIgnore)  {
        fputs("\tIGNORE", fp);
    }
    fputs("\n", fp);
}

int main(int argc, char **argv) {
    if (argc > 2) {
        fprintf(stderr, "cat a.fq | %s [outputSeq=1] >out.lsam\n", argv[0]);
        exit(1);
    }

    if (argc == 2) {
        outputSeq = atoi(argv[1]);
    }

    kxseq<int> seq(fileno(stdin));
    bool has_last = false;

    while (seq.read() >= 0) {
        trim_readno(seq.name());
        if (has_last) {
             if (name1 == seq.name()) {
             	print_lsam_line(name1, comm1, seq1, qual1, 1);
                print_lsam_line(seq.name(), seq.comment(), seq.seq(), seq.qual(), 2);
                has_last = false;
             } else {
                print_lsam_line(name1, comm1, seq1, qual1, 0);
                update1(seq);
             }
        } else {
            has_last = true;
            update1(seq);
        }
    }

    if (has_last) {
        print_lsam_line(name1, comm1, seq1, qual1, 0);
    }

    return 0;
}