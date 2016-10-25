/*
 * align.h
 *
 *  Created on: Oct 12, 2016
 *      Author: isovic
 */

#ifndef ALIGN_IS_H_
#define ALIGN_IS_H_

#include <stdio.h>
#include <stdint.h>
#include <string>
#include <vector>
#include <algorithm>

namespace is {

typedef enum {
  kGlobal,
  kLocal
} AlignType;


// const int ALN_MOVE_EQ = 0;
// const int ALN_MOVE_D = 2;
// const int ALN_MOVE_I = 1;
// const int ALN_MOVE_X = 3;

const int8_t ALN_OP_EQ = 0;
const int8_t ALN_OP_X = 3;
const int8_t ALN_OP_D = 2;
const int8_t ALN_OP_I = 1;
const int8_t ALN_OP_S = 4;
const int8_t ALN_OP_H = 5;
const int8_t ALN_OP_NOP = 6;
const int8_t ALN_OP_M = 7;
const char ALN_OP_TO_CHAR[] = "=IDXSH";
const char ALN_OP_TO_MATCH[] = "|  X-~";
const char ALN_OP_TO_BASIC_CHAR[] = "MIDMSH";
const int8_t ALN_OP_TO_EDLIB[] = {ALN_OP_EQ, ALN_OP_I, ALN_OP_D, ALN_OP_X, ALN_OP_S, ALN_OP_H, ALN_OP_NOP};

const int32_t MINUS_INF = std::numeric_limits<int32_t>::min() + 1000000; // Allow a margin to skip overflow.

class CigarOp {
public:
  CigarOp() : op(ALN_OP_NOP), count(0) { };
  CigarOp(int8_t _op, int32_t _count) : op(_op), count(_count) { };

  int8_t op;
  int32_t count;
};

class Penalties {
 public:
  Penalties(int32_t _match, int32_t _mismatch, int32_t _gopen, int32_t _gext)
      : match(_match),
        mismatch(_mismatch),
        gopen(_gopen),
        gext(_gext) {
  }
  ;
  int32_t match, mismatch, gopen, gext;
};

// If any global margin is true, then the corresponding will be penalized.
// Concretely, if top/left are true, then the first row/column will be initialized
// to the multiple of the gap extend penalty in global alignment.
// If bottom is false, the maximum of last row will be found instead of taking
// the bottom right corner for global alignment.
// If right is false, the maximum of last column will be found instead of taking
// the bottom right corner for global alignment.
class GlobalMargins {
 public:
  GlobalMargins()
      : top(true),
        left(true),
        bottom(true),
        right(true) {
  }
  GlobalMargins(bool _top, bool _left, bool _bottom, bool _right)
      : top(_top),
        left(_left),
        bottom(_bottom),
        right(_right) {
  }
  bool top, left, bottom, right;
};

class Align {
 public:
  Align(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p,
        AlignType aln_type, GlobalMargins gm);
  Align(const std::string &q, const std::string &t, Penalties p,
        AlignType aln_type, GlobalMargins gm);
  ~Align();

  int GetAlignment(int32_t &q_start, int32_t &q_end, int32_t &t_start, int32_t &t_end, int32_t &score, std::vector<CigarOp> &cigar);
  // int GetCigar(std::vector<CigarOp> &cigar);
  // int GetCigarAsString(std::string &cigar_string);
  void FormatAlignment(const std::string &q, const std::string &t, std::string &alnq, std::string &alnt, std::string &alnm);
  void Verbose(const std::string &q, const std::string &t, std::ostream &os);
  
 private:
  Align(const Align& op) = delete;
  Align& operator=(const Align& op) = delete;

  int Align_(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p, GlobalMargins gm, AlignType aln_type,
                   int32_t &q_start, int32_t &q_end, int32_t &t_start, int32_t &t_end, int32_t &score , std::vector<is::CigarOp> &cigar);
  int AlignGlobal_(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p, GlobalMargins gm,
                   int32_t &q_start, int32_t &q_end, int32_t &t_start, int32_t &t_end, int32_t &score, std::vector<is::CigarOp> &cigar);
  int AlignLocal_(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p);
  int Traceback_(const char* q, int64_t ql, const char* t, int64_t tl,
		  std::vector<std::vector<int32_t> > &dir, int32_t row, int32_t col, std::vector<CigarOp> &cigar);

  Penalties p_;
  GlobalMargins gm_;
  AlignType aln_type_;
  std::vector<is::CigarOp> cigar_;
  int32_t q_start_, q_end_, t_start_, t_end_;
  int32_t score_;
};

std::string CigarToString(const std::vector<CigarOp> &cigar);
std::string CigarToBasicString(const std::vector<CigarOp> &cigar);
void CigarToEdlibAln(const std::vector<CigarOp> &cigar, std::vector<int8_t>& alignment);
void CigarToAlignment(const char* q, int64_t ql, int64_t q_start, int64_t q_end,
            const char* t, int64_t tl, int64_t t_start, int64_t t_end, AlignType aln_type,
            const std::vector<CigarOp> &cigar, std::string &alnq, std::string &alnt, std::string &alnm);

}

#endif /* ALIGN_H_ */
