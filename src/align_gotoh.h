/*
 * align.h
 *
 *  Created on: Oct 12, 2016
 *      Author: isovic
 */

#ifndef ALIGN_GOTOH_IS_H_
#define ALIGN_GOTOH_IS_H_

#include <stdio.h>
#include <stdint.h>
#include <string>
#include <vector>
#include <algorithm>

#include "align_utility.h"

namespace is {

class AlignGotoh {
 public:
  AlignGotoh(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p,
        AlignType aln_type, GlobalMargins gm);
  AlignGotoh(const std::string &q, const std::string &t, Penalties p,
        AlignType aln_type, GlobalMargins gm);
  ~AlignGotoh();

  int GetAlignment(int32_t &q_start, int32_t &q_end, int32_t &t_start, int32_t &t_end, int32_t &score, std::vector<CigarOp> &cigar);
  // int GetCigar(std::vector<CigarOp> &cigar);
  // int GetCigarAsString(std::string &cigar_string);
  void FormatAlignment(const std::string &q, const std::string &t, std::string &alnq, std::string &alnt, std::string &alnm);
  void Verbose(const std::string &q, const std::string &t, std::ostream &os);
  
 private:
  AlignGotoh(const AlignGotoh& op) = delete;
  AlignGotoh& operator=(const AlignGotoh& op) = delete;

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

}

#endif /* ALIGN_H_ */
