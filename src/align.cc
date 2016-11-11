/*
 * align.cpp
 *
 *  Created on: Oct 12, 2016
 *      Author: isovic
 */

#include "align.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <cmath>
#include "timer.hpp"

namespace is {

Align::Align(const std::string& q, const std::string& t, Penalties p,
             AlignType aln_type, GlobalMargins gm)
    : Align(q.c_str(), q.size(), t.c_str(), t.size(), p, aln_type, gm) {
}

Align::Align(const char* q, int64_t ql, const char* t, int64_t tl, Penalties p,
             AlignType aln_type, GlobalMargins gm) :  p_(p), gm_(gm), aln_type_(aln_type),
		     q_start_(0), q_end_(0), t_start_(0), t_end_(0), score_(0)  {
  // Align_(q, ql, t, tl, p, gm, aln_type, q_start_, q_end_, t_start_, t_end_, score_, cigar_);
}

Align::~Align() {

}

int Align::GetAlignment(int32_t &q_start, int32_t &q_end, int32_t &t_start, int32_t &t_end, int32_t &score, std::vector<CigarOp> &cigar) {
	cigar = cigar_;
	q_start = q_start_;
	q_end = q_end_;
	t_start = t_start_;
	t_end = t_end_;
	score = score_;
	return 0;
}

void Align::FormatAlignment(const std::string &q, const std::string &t, std::string &alnq, std::string &alnt, std::string &alnm) {
	CigarToAlignment(q.c_str(), q.size(), q_start_, q_end_, t.c_str(), t.size(), t_start_, t_end_, aln_type_, cigar_, alnq, alnt, alnm);
}

void Align::Verbose(const std::string &q, const std::string &t, std::ostream &os) {
  // Debug output and conversion.
  std::string alnq, alnt, alnm;
  CigarToAlignment(q.c_str(), q.size(), q_start_, q_end_, t.c_str(), t.size(), t_start_, t_end_, aln_type_, cigar_, alnq, alnt, alnm);
  os << alnt << std::endl;
  os << alnm << std::endl;
  os << alnq << std::endl;
  os << CigarToString(cigar_) << std::endl;
  os << CigarToBasicString(cigar_) << std::endl;
  os << "q = [" << q_start_ << ", " << q_end_ << "]" << std::endl;
  os << "t = [" << t_start_ << ", " << t_end_ << "]" << std::endl;
  os << "score = " << score_ << std::endl;
}

int Align::Align_(const char* q, int64_t ql, const char* t, int64_t tl, Penalties p, GlobalMargins gm, AlignType aln_type,
                       int32_t &q_start, int32_t &q_end, int32_t &t_start, int32_t &t_end, int32_t &score, std::vector<is::CigarOp> &cigar) {
  fprintf (stderr, "Not yet implemented!\n");
  return 1;
}

int Align::AlignGlobal_(const char* q, int64_t ql, const char* t, int64_t tl, Penalties p, GlobalMargins gm,
                       int32_t &q_start, int32_t &q_end, int32_t &t_start, int32_t &t_end, int32_t &score, std::vector<is::CigarOp> &cigar) {
  fprintf (stderr, "Not yet implemented!\n");
  return 1;
}

int Align::AlignLocal_(const char* q, int64_t ql, const char* t, int64_t tl,
                      Penalties p) {
  fprintf (stderr, "Not yet implemented!\n");
  return 1;
}

int Align::Traceback_(const char* q, int64_t ql, const char* t, int64_t tl,
					  std::vector<std::vector<int32_t> > &dir, int32_t row, int32_t col, std::vector<is::CigarOp> &cigar) {

	cigar.clear();
	cigar.reserve((int32_t) std::sqrt(row*row + col*col));	// An approximation, to save from reallocating every move.

	q_start_ = 0;
	q_end_ = ql;
	t_start_ = 0;
	t_end_= tl;

	// if ((ql - row) > 0) { cigar.push_back(is::CigarOp(ALN_OP_I, (ql - row))); q_end_ = (ql - row); }
	// if ((tl - col) > 0) { cigar.push_back(is::CigarOp(ALN_OP_D, (tl - col))); t_end_ = (tl - col); }
	if ((ql - row) > 0) { q_end_ = (ql - row); }
	if ((tl - col) > 0) { t_end_ = (tl - col); }

  	int8_t op = dir[row][col];
	if (op == ALN_OP_M || op == ALN_OP_EQ || op == ALN_OP_X) {
		op = (q[row-1] == t[col-1]) ? ALN_OP_EQ : ALN_OP_X;
	}
	  is::CigarOp cigarop(op, 0);

	  while (row > 0 && col > 0) {
	  	op = dir[row][col];
		if (op == ALN_OP_M || op == ALN_OP_EQ || op == ALN_OP_X) {
			op = (q[row-1] == t[col-1]) ? ALN_OP_EQ : ALN_OP_X;
		}

		if (op == cigarop.op) {
		  cigarop.count += 1;
		} else {
			cigar.push_back(cigarop);
			cigarop.op = op;
			cigarop.count = 1;
		}

		if (op == ALN_OP_EQ || op == ALN_OP_X) {
			row -= 1; col -= 1;
		} else if (op == ALN_OP_I || op == ALN_OP_S) {
			row -= 1;
		} else if (op == ALN_OP_D) {
			col -= 1;
		}
	  }

	  if (cigarop.count > 0) {
	  	cigar.push_back(cigarop);
	  }

	  // if (row > 0) { cigar.push_back(CigarOp(ALN_OP_I, row)); q_start_ = row; }
	  // if (col > 0) { cigar.push_back(CigarOp(ALN_OP_D, col)); t_start_ = col; }
	  if (row > 0) { q_start_ = row; }
	  if (col > 0) { t_start_ = col; }

	  std::reverse(cigar.begin(), cigar.end());

	  return 0;
}

}
