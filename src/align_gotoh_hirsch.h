/*
 * align.h
 *
 *  Created on: Oct 12, 2016
 *      Author: isovic
 */

#ifndef ALIGN_GOTOH2_IS_H_
#define ALIGN_GOTOH2_IS_H_

#include <stdio.h>
#include <stdint.h>
#include <string>
#include <vector>
#include <algorithm>

#include "align_utility.h"

#include "align.h"

namespace is {

class AlignGotohHirsch : public Align {
 public:
  AlignGotohHirsch(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p,
        AlignType aln_type, GlobalMargins gm);
  AlignGotohHirsch(const std::string &q, const std::string &t, Penalties p,
        AlignType aln_type, GlobalMargins gm);
  ~AlignGotohHirsch();

 // private:
 protected:
  AlignGotohHirsch(const AlignGotohHirsch& op) = delete;
  AlignGotohHirsch& operator=(const AlignGotohHirsch& op) = delete;

  int Align_(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p, GlobalMargins gm, AlignType aln_type,
                   int32_t &q_start, int32_t &q_end, int32_t &t_start, int32_t &t_end, int32_t &score , std::vector<is::CigarOp> &cigar);
  int AlignGlobal_(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p, GlobalMargins gm,
                   int32_t &q_start, int32_t &q_end, int32_t &t_start, int32_t &t_end, int32_t &score, std::vector<is::CigarOp> &cigar);
  int AlignLocal_(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p);

  // M is a two row matrix for calculating the alignment (main matrix).
  // V is also two row, the vertical Gotoh matrix.
  // H0 is the initial row penalty (zeroth element of the row).
  int AlignBlock_(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p, GlobalMargins gm,
                  MatrixType &M, MatrixType &V, cell_t H0);

};

}

#endif /* ALIGN_H_ */
