/* The MIT License
   Copyright (c) 2016 Dinghua Li <voutcn@gmail.com>
   
   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:
   
   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.
   
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

   This is a C++ reinplement of kseq.h:
      https://github.com/attractivechaos/klib/blob/master/kseq.h
*/

#ifndef V_KX_SEQ_H
#define V_KX_SEQ_H

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <algorithm>

#define HAVE_ZLIB
#ifdef HAVE_ZLIB
#include <zlib.h>
#endif

namespace kx {

#ifndef V_KX_STREAM_T
#define V_KX_STREAM_T

enum Seperator {
    kSepSpace = 0,
    kSepTab = 1,
    kSepLine = 2
};

static const int kEOF = -1;

template <typename file_t> struct kxreadfunc {};

template <> struct kxreadfunc<int> {
    size_t operator() (int fp, void *buf, size_t size) {
        return read(fp, buf, size);
    }   
};

#ifdef HAVE_ZLIB
template <> struct kxreadfunc<gzFile> {
    int operator() (gzFile fp, void *buf, int size) {
        return gzread(fp, buf, size);
    }
};
#endif

template <typename file_t, typename read_t = kxreadfunc<file_t> >
class kxstream {
private:
    static const int kBufferSize = 4096;
    file_t file_;
    read_t read_;
    char *buf_;
    int begin_, end_;
    bool eof_;

public:

    kxstream(file_t f)
        : file_(f), read_(read_t()), buf_(NULL),
          begin_(0), end_(0), eof_(false) {
        buf_ = new char[kBufferSize];
    }
    ~kxstream() { delete[] buf_; }

    bool eof() const {
        return eof_;
    }

    int getc() {
        if (eof_) { return kEOF; }
        if (begin_ >= end_) {
            begin_ = 0;
            end_ = read_(file_, buf_, kBufferSize);
            if (end_ == 0) { eof_ = true; return kEOF; }
        }
        return static_cast<int>(buf_[begin_++]);
    }

    int get_until(int delimiter, std::string &str, bool append = false) {
        bool got_any = false;
        int i = 0;
        int dret = 0;
        if (!append) str.clear();

        while (true) {
            if (begin_ >= end_) {
                if (!eof_) {
                    begin_ = 0;
                    end_ = read_(file_, buf_, kBufferSize);
                    if (end_ == 0) { eof_ = true; break; }
                } else { break; }
            }

            if (delimiter == kSepLine) {
                for (i = begin_; i < end_; ++i) {
                    if (buf_[i] == '\n') break;
                }
            } else if (delimiter > kSepLine) {
                for (i = begin_; i < end_; ++i) {
                    if (buf_[i] == delimiter) break;
                }
            } else if (delimiter == kSepSpace) {
                for (i = begin_; i < end_; ++i) {
                    if (isspace(buf_[i])) break;
                }
            } else if (delimiter == kSepTab) {
                for (i = begin_; i < end_; ++i) {
                    if (isspace(buf_[i]) && buf_[i] != ' ') break;
                }
            }

            got_any = true;
            str.append(buf_ + begin_, i - begin_);
            begin_ = i + 1;

            if (i < end_) {
                dret = buf_[i];
                break;
            }
        }

        if (!got_any && eof_) return kEOF;
        if (delimiter == kSepLine && str.length() > 1 && str[str.length() - 1] == '\r') str.resize(str.length() - 1);
        return dret;
    }
};

#endif

static const int kTruncatedQual = -2;

template <typename file_t, typename read_t = kxreadfunc<file_t> >
class kxseq {
private:
    kxstream<file_t, read_t> stream_;
    std::string name_, comment_, seq_, qual_;
    int last_char_;

public:
    kxseq(file_t f): stream_(f), last_char_(0) {}
    ~kxseq() {}
    const std::string &name() const { return name_; }
    const std::string &comment() const { return comment_; }
    const std::string &seq() const { return seq_; }
    const std::string &qual() const { return qual_; }
    std::string &name() { return name_; }
    std::string &comment() { return comment_; }
    std::string &seq() { return seq_; }
    std::string &qual() { return qual_; }

    int read() {
        int c;
        if (last_char_ == 0) {
            while ((c = stream_.getc()) != kEOF && c != '>' && c != '@') { continue; }
            if (c == kEOF) return kEOF; // end of the file
            last_char_ = c;
        }

        name_.clear(); comment_.clear(); seq_.clear(); qual_.clear();

        if ((c = stream_.get_until(kSepSpace, name_)) == kEOF) return kEOF; // got nothing
        if (c != '\n') stream_.get_until(kSepLine, comment_);

        while ((c = stream_.getc()) != kEOF && c != '>' && c != '+' && c != '@') {
            if (c == '\n') continue;
            seq_.push_back(c);
            stream_.get_until(kSepLine, seq_, true);
        }

        if (c == '>' || c == '@') last_char_ = c;
        if (c != '+') return seq_.length();

        while ((c = stream_.getc()) != kEOF && c != '\n') { continue; }
        if (c == kEOF) return kTruncatedQual;

        while (stream_.get_until(kSepLine, qual_, true) != kEOF && qual_.length() < seq_.length()) {
            continue;
        }

        last_char_ = 0;
        if (qual_.length() != seq_.length()) return kTruncatedQual;
        return seq_.length();
    }
};

}

#endif