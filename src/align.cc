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

namespace is {

Align::Align(const std::string& q, const std::string& t, Penalties p,
             AlignType aln_type, GlobalMargins gm)
    : Align(q.c_str(), q.size(), t.c_str(), t.size(), p, aln_type, gm) {
}

Align::Align(const char* q, int64_t ql, const char* t, int64_t tl, Penalties p,
             AlignType aln_type, GlobalMargins gm) :  p_(p), gm_(gm), aln_type_(aln_type),
		     q_start_(0), q_end_(0), t_start_(0), t_end_(0), score_(0)  {
  Align_(q, ql, t, tl, p, gm, aln_type, q_start_, q_end_, t_start_, t_end_, score_, cigar_);
//  if (aln_type == kGlobal) {
//    AlignGlobal_(q, ql, t, tl, p, gm, q_start_, q_end_, t_start_, t_end_, cigar_);
//  } else if (aln_type == kLocal) {
//    AlignLocal_(q, ql, t, tl, p);
//  }
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

template <class T>
void PrintMatrix(std::vector<std::vector<T> > m) {
  for (int32_t i=0; i<m.size(); i++) {
    for (int32_t j=0; j<m[i].size(); j++) {
      std::cout << m[i][j] << "\t";
    }
    std::cout << std::endl;
  }
}

int Align::Align_(const char* q, int64_t ql, const char* t, int64_t tl, Penalties p, GlobalMargins gm, AlignType aln_type,
                       int32_t &q_start, int32_t &q_end, int32_t &t_start, int32_t &t_end, int32_t &score, std::vector<is::CigarOp> &cigar) {

  // Storing only two lines for the main (M) and the vertical matrix (V).
  std::vector<std::vector<int32_t> > M(2, std::vector<int32_t>(tl+1, 0));
  std::vector<std::vector<int32_t> > V(2, std::vector<int32_t>(tl+1, 0));
  // Only current and previous element of the horizontal matrix are needed.
  std::vector<int32_t> H(2, 0);
  // The direction backtrack matrix is stored entirely.
  std::vector<std::vector<int32_t> > dir(ql+1, std::vector<int32_t>(tl+1, 0));

//  int32_t w[] = {p.match, p.mismatch};    // Match score and mismatch penalty, for easier lookup using a profile.
//  int32_t w_op[] = {ALN_OP_EQ, ALN_OP_X}; // Operation lookup (for profile).

  // Penalize the first row.
  if (aln_type == kGlobal && gm.top) {
    for (int32_t i=0; i<(tl+1); i++) {
      M[0][i] = i * p.gext;
      V[0][i] = MINUS_INF;
    }
  }

  // Variables for tracking the maximum (for SW).
  int32_t bt_row = 0;  // Backtrace start column.
  int32_t bt_col = 0;  // Backtrace start row.
  int32_t bt_val = 0;  // Value of M at the bt coordinates.

  for (int32_t i=1; i<(ql+1); i++) {
    if (aln_type == kGlobal && gm.top) { H[0] = MINUS_INF; }
    // Penalize the first column.
    if (aln_type == kGlobal && gm.left) {	M[1][0] = i * p.gext; }

    for (int32_t j=1; j<(tl+1); j++) {
      V[1][j] = std::max((V[0][j] + p.gext), M[0][j] + p.gopen);
      H[1] = std::max(H[0] + p.gext, M[1][j-1] + p.gopen);

      int32_t diag = M[0][j-1] + ((q[i-1] == t[j-1]) ? p.match : p.mismatch);
      dir[i][j] = ALN_OP_M;  // Match/Mismatch. Needs to be replaced in the traceback by a correct op (faster this way).
      M[1][j] = diag;
      if (V[1][j] > M[1][j]) { M[1][j] = V[1][j]; dir[i][j] = ALN_OP_I; }  // Insertion.
      if (H[1] > M[1][j]) { M[1][j] = H[1]; dir[i][j] = ALN_OP_D; }  // Deletion.

      // Local alignment - reset scores smaller than 0.
      if (aln_type == kLocal && M[1][j] < 0) { M[1][j] = 0; }
      // Local alignment - find maximum.
      if (M[1][j] > bt_val) { bt_val = M[1][j]; bt_row = i; bt_col = j; }
    }
    V[0] = V[1];
    M[0] = M[1];
    H[0] = H[1];
  }

  // General case for global alignment. Start from the bottom right corner.
  if (aln_type == kGlobal) {
    bt_row = ql;
    bt_col = tl;
    bt_val = M[1][tl];
  }

  // Handle semiglobal right margin. Allow the start position at any row of the last column.
  if (aln_type == kGlobal && !gm.right) {
    for (int32_t i=0; i<(ql+1); i++) {
      if (M[i][tl] > bt_val) {
        bt_row = i;
        bt_val = M[i][tl];
      }
    }
  }

  // Handle semiglobal bottom margin. Allow the start position at any column of the last row.
  if (aln_type == kGlobal && !gm.bottom) {
    for (int32_t i=0; i<(tl+1); i++) {
      if (M[ql][i] > bt_val) {
        bt_col = i;
        bt_val = M[ql][i];
      }
    }
  }

  score = bt_val;

  // Perform traceback.
  Traceback_(q, ql, t, tl, dir, bt_row, bt_col, cigar);

  return 0;
}

int Align::AlignGlobal_(const char* q, int64_t ql, const char* t, int64_t tl, Penalties p, GlobalMargins gm,
                       int32_t &q_start, int32_t &q_end, int32_t &t_start, int32_t &t_end, int32_t &score, std::vector<is::CigarOp> &cigar) {
  return 0;
}

int Align::AlignLocal_(const char* q, int64_t ql, const char* t, int64_t tl,
                      Penalties p) {
  return 0;
}

int Align::Traceback_(const char* q, int64_t ql, const char* t, int64_t tl,
					  std::vector<std::vector<int32_t> > &dir, int32_t row, int32_t col, std::vector<is::CigarOp> &cigar) {

	cigar.clear();
	cigar.reserve((int32_t) sqrt(row*row + col*col));	// An approximation, to save from reallocating every move.

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

std::string CigarToString(const std::vector<CigarOp> &cigar) {
	std::stringstream ss;
	for (uint32_t i=0; i<cigar.size(); i++) {
		ss << cigar[i].count << ALN_OP_TO_CHAR[cigar[i].op]	;
	}
	return ss.str();
}

std::string CigarToBasicString(const std::vector<CigarOp> &cigar) {
	std::stringstream ss;
	std::vector<CigarOp> cigar_basic(cigar);

	int32_t offset = 0;
	for (int32_t i=1; i<(int32_t) cigar_basic.size(); i++) {
		if (ALN_OP_TO_BASIC_CHAR[cigar_basic[i].op] == ALN_OP_TO_BASIC_CHAR[cigar_basic[i-offset-1].op]) {
			cigar_basic[i-offset-1].count += cigar_basic[i].count;
			offset += 1;
			continue;
		}
		cigar_basic[i-offset] = cigar_basic[i];
	}
	for (int32_t i=0; i<(((int32_t) cigar_basic.size())-offset); i++) {
		ss << cigar_basic[i].count << ALN_OP_TO_BASIC_CHAR[cigar_basic[i].op]	;
	}
	return ss.str();
}

void CigarToEdlibAln(const std::vector<CigarOp> &cigar, std::vector<int8_t>& alignment) {
	alignment.clear();
	int32_t aln_len = 0;
	for (uint32_t i=0; i<cigar.size(); i++) {
		aln_len += cigar[i].count;
	}
	alignment.reserve(aln_len);
	for (uint32_t i=0; i<cigar.size(); i++) {
		std::vector<int8_t> temp(cigar[i].count, ALN_OP_TO_EDLIB[cigar[i].op]);
		alignment.insert(alignment.end(), temp.begin(), temp.end());
	}
}

void CigarToAlignment(const char* q, int64_t ql, int64_t q_start, int64_t q_end,
					  const char* t, int64_t tl, int64_t t_start, int64_t t_end, AlignType aln_type,
					  const std::vector<CigarOp> &cigar, std::string &alnq, std::string &alnt, std::string &alnm) {
	std::stringstream ssq, sst, ssm;

	int32_t posq = 0, post = 0;

  if (aln_type == kGlobal) {
    for (int32_t i=0; i<t_start; i++) {
      ssq << "-"; sst << t[post++]; ssm << " ";
    }
    for (int32_t i=0; i<q_start; i++) {
      ssq << q[posq++]; sst << "-"; ssm << " ";
    }
  }
  posq = q_start;
  post = t_start;

	for (uint32_t i=0; i<cigar.size(); i++) {
		if (cigar[i].op == ALN_OP_EQ || cigar[i].op == ALN_OP_X) {
			for (int32_t j=0; j<cigar[i].count; j++) {
				ssq << q[posq++];	sst << t[post++];	ssm << ALN_OP_TO_MATCH[cigar[i].op];
			}
		}
		if (cigar[i].op == ALN_OP_I) {
			for (int32_t j=0; j<cigar[i].count; j++) {
				ssq << q[posq++];	sst << "-";	ssm << " ";
			}
		}
		if (cigar[i].op == ALN_OP_D) {
			for (int32_t j=0; j<cigar[i].count; j++) {
				ssq << "-";	sst << t[post++];	ssm << " ";
			}
		}
		if (cigar[i].op == ALN_OP_S) {
			for (int32_t j=0; j<cigar[i].count; j++) {
				ssq << q[posq++];	sst << "-";	ssm << "-";
			}
		}

	}

  if (aln_type == kGlobal) {
    for (; post < tl; post++) {
      ssq << "-"; sst << t[post]; ssm << " ";
    }
    for (; posq < ql; posq++) {
      ssq << q[posq]; sst << "-"; ssm << " ";
    }
  }

	alnq = ssq.str();
	alnt = sst.str();
	alnm = ssm.str();
}

}
