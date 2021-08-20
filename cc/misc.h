#ifndef MISC_H__
#define MISC_H__

#include <string>
#include <vector>
#include <algorithm>
#include "taxonomy.h"

inline std::string trimno(const std::string &s) {
    if (s.size() > 2U && s[s.size() - 2] == '/' && (s[s.size() - 1] == '1' || s[s.size() - 1] == '2')) {
        return s.substr(0, s.size() - 2);
    }
    return s;
}

std::vector<std::string> splitBy(const std::string &s, char delimiter) {
    std::vector<std::string> ret;
    for (size_t i = 0, j; i < s.size(); i = j + 1) {
        j = i;
        while (j < s.size() && s[j] != delimiter) {
            ++j;
        }
        ret.push_back(s.substr(i, j - i));
    }
    return ret;
}

std::string getCorrectAcc(const std::string &ntHeader) {
    size_t spliterPos = ntHeader.find_first_of('|');
    if (spliterPos == std::string::npos) {
        // new format header, see http://www.ncbi.nlm.nih.gov/news/09-17-2014-simple-FASTA-headers-genomes-FTP/
        return ntHeader;
    } else if (ntHeader.substr(0, spliterPos) == "gi") {
        // old header, skip the GI, and then the db flag
        spliterPos = ntHeader.find_first_of('|', spliterPos + 1);
        spliterPos = ntHeader.find_first_of('|', spliterPos + 1);
        size_t spliterPos2 = ntHeader.find_first_of('|', spliterPos + 1);

        return ntHeader.substr(spliterPos + 1, spliterPos2 - spliterPos - 1);
    } else {
        return ntHeader;
    }
}


std::vector<std::pair<std::string, double> > splitAcc(const std::string &hits) {
    std::vector<std::pair<std::string, double> > ret;
    if (hits.length() == 0 || hits == "*") { return ret; }

    auto rec = splitBy(hits, ';');
    for (size_t i = 0; i < rec.size(); ++i) {
        if (rec[i].size() == 0) continue;
        auto score_and_acc = splitBy(rec[i], ',');
        score_and_acc[1] = getCorrectAcc(score_and_acc[1]);
        ret.push_back(std::make_pair(score_and_acc[1], atof(score_and_acc[0].c_str())));
    }
    return ret;
}

#endif