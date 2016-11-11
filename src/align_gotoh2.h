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

class AlignGotoh2 : public Align {
 public:
  AlignGotoh2(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p,
        AlignType aln_type, GlobalMargins gm);
  AlignGotoh2(const std::string &q, const std::string &t, Penalties p,
        AlignType aln_type, GlobalMargins gm);
  ~AlignGotoh2();

 // private:
 protected:
  AlignGotoh2(const AlignGotoh2& op) = delete;
  AlignGotoh2& operator=(const AlignGotoh2& op) = delete;

  int Align_(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p, GlobalMargins gm, AlignType aln_type,
                   int32_t &q_start, int32_t &q_end, int32_t &t_start, int32_t &t_end, int32_t &score , std::vector<is::CigarOp> &cigar);
  int AlignGlobal_(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p, GlobalMargins gm,
                   int32_t &q_start, int32_t &q_end, int32_t &t_start, int32_t &t_end, int32_t &score, std::vector<is::CigarOp> &cigar);
  int AlignLocal_(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p);

};

}

#endif /* ALIGN_H_ */
