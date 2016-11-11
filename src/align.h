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

#include "align_utility.h"

namespace is {

typedef int32_t cell_t;
typedef std::vector<std::vector<cell_t> > MatrixType;

class Align {
 public:
  // This constructor simply forwards the parameters to the constructor which takes C-style strings (char arrays).
  Align(const std::string &q, const std::string &t, Penalties p,
        AlignType aln_type, GlobalMargins gm);
  // Sets parameters and aligns the sequences.
  Align(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p,
        AlignType aln_type, GlobalMargins gm);
  virtual ~Align();

  // Gets all the results at the same time.
  int GetAlignment(int32_t &q_start, int32_t &q_end, int32_t &t_start, int32_t &t_end, int32_t &score, std::vector<CigarOp> &cigar);

  // Formats the alignment in a BLAST-like string fashion.
  void FormatAlignment(const std::string &q, const std::string &t, std::string &alnq, std::string &alnt, std::string &alnm);

  // Verbose the alignment results to an output stream.
  void Verbose(const std::string &q, const std::string &t, std::ostream &os);
  
 // private:
 protected:
  // Function Align_ automatically calls AlignGlobal_ or AlignLocal_. Can be reimplemented.
  virtual int Align_(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p, GlobalMargins gm, AlignType aln_type,
                   int32_t &q_start, int32_t &q_end, int32_t &t_start, int32_t &t_end, int32_t &score , std::vector<is::CigarOp> &cigar);

  // Virtual methods to be implemented in the derived classes.
  virtual int AlignGlobal_(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p, GlobalMargins gm,
                   int32_t &q_start, int32_t &q_end, int32_t &t_start, int32_t &t_end, int32_t &score, std::vector<is::CigarOp> &cigar);
  virtual int AlignLocal_(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p);

  // Forbid the copy constructor.
  Align(const Align& op) = delete;
  // Forbit assignment.
  Align& operator=(const Align& op) = delete;

 // protected:

  // A generic traceback which works for both the global and local alignment modes.
  int Traceback_(const char* q, int64_t ql, const char* t, int64_t tl,
		  std::vector<std::vector<int32_t> > &dir, int32_t row, int32_t col, std::vector<CigarOp> &cigar);

  // Alignment parameters.
  Penalties p_;
  GlobalMargins gm_;
  AlignType aln_type_;

  // Alignment results.
  std::vector<is::CigarOp> cigar_;
  int32_t q_start_, q_end_, t_start_, t_end_;
  int32_t score_;
};

}

#endif /* ALIGN_H_ */
