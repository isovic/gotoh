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
             AlignType aln_type, GlobalMargins gm) :  q_(q), ql_(ql), t_(t), tl_(tl), p_(p), gm_(gm)  {
  if (aln_type == kGlobal) {
    AlignGlobal_(q, ql, t, tl, p, gm);
  } else if (aln_type == kLocal) {
    AlignLocal_(q, ql, t, tl, p);
  }
}

Align::~Align() {

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

int Align::AlignGlobal_(const char* q, int64_t ql, const char* t, int64_t tl,
                       Penalties p, GlobalMargins gm) {

  std::vector<std::vector<int32_t> > M(ql+1, std::vector<int32_t>(tl+1, 0));
  std::vector<std::vector<int32_t> > dir(ql+1, std::vector<int32_t>(tl+1, 0));
  std::vector<std::vector<int32_t> > V(2, std::vector<int32_t>(tl+1, 0));
  std::vector<int32_t> H(tl+1, 0);

  int32_t w[] = {p.match, p.mismatch};    // Match score and mismatch penalty, for easier lookup using a profile.
  int32_t w_op[] = {ALN_OP_EQ, ALN_OP_X}; // Operation lookup (for profile).



  // Penalize the first column.
  if (gm.left) {
    for (int32_t i=0; i<(ql+1); i++) { M[i][0] = i * p.gext; }
  }
  // Penalize the first row.
  if (gm.top) {
    for (int32_t i=0; i<(tl+1); i++) {
      M[0][i] = i * p.gext;
      V[0][i] = MINUS_INF;
    }
  }

  std::abs(-5);

  for (int32_t i=1; i<(ql+1); i++) {
    if (gm.top) { H[0] = MINUS_INF; }

    for (int32_t j=1; j<(tl+1); j++) {
      V[1][j] = std::max((V[0][j] + p.gext), M[i-1][j] + p.gopen);
      H[j] = std::max(H[j-1] + p.gext, M[i][j-1] + p.gopen);

      int32_t diag = M[i-1][j-1] + ((q[i-1] == t[j-1]) ? p.match : p.mismatch);
      dir[i][j] = (q[i-1] == t[j-1]) ? ALN_OP_EQ : ALN_OP_X;  // Match/Mismatch.
      M[i][j] = diag;
      if (V[1][j] > M[i][j]) { M[i][j] = V[1][j]; dir[i][j] = ALN_OP_I; }  // Insertion.
      if (H[j] > M[i][j]) { M[i][j] = H[j]; dir[i][j] = ALN_OP_D; }  // Deletion.
    }
    V[0] = V[1];
  }

  int32_t bt_row = ql;  // Backtrace start column.
  int32_t bt_col = tl;  // Backtrace start row.
  int32_t bt_val = M[ql][tl]; // Value of M at the bt coordinates.

  // Handle semiglobal right margin.
  if (gm.right) {
    for (int32_t i=0; i<(ql+1); i++) {
      if (M[i][tl] > bt_val) {
        bt_row = i;
        bt_val = M[i][tl];
      }
    }
  }

  // Handle semiglobal bottom margin.
  if (gm.bottom) {
    for (int32_t i=0; i<(tl+1); i++) {
      if (M[ql][i] > bt_val) {
        bt_col = i;
        bt_val = M[ql][i];
      }
    }
  }

  // Perform traceback.
  std::vector<is::CigarOp> cigar;
  Traceback_(q, ql, t, tl, M, dir, bt_row, bt_col, cigar);

  // Debug output and conversion.
  std::string alnq, alnt, alnm;
  CigarToAlignment(q, ql, t, tl, cigar, alnq, alnt, alnm);
  printf ("%s\n", alnt.c_str());
  printf ("%s\n", alnm.c_str());
  printf ("%s\n", alnq.c_str());
  printf ("%s\n", CigarToString(cigar).c_str());
  printf ("%s\n", CigarToBasicString(cigar).c_str());

  return 0;
}

int Align::AlignLocal_(const char* q, int64_t ql, const char* t, int64_t tl,
                      Penalties p) {
  return 0;
}

int Align::Traceback_(const char* q, int64_t ql, const char* t, int64_t tl,
		std::vector<std::vector<int32_t> > &M, std::vector<std::vector<int32_t> > &dir, int32_t row, int32_t col, std::vector<is::CigarOp> &cigar) {

	  cigar.clear();
	  cigar.reserve((int32_t) sqrt(row*row + col*col));	// An approximation, to save from reallocating every move.

	  if ((ql - row) > 0) { cigar.push_back(is::CigarOp(ALN_OP_I, (ql - row))); }
	  if ((tl - col) > 0) { cigar.push_back(is::CigarOp(ALN_OP_D, (tl - col))); }

	  is::CigarOp cigarop(dir[row][col], 0);

	  while (row > 0 && col > 0) {
		  if (dir[row][col] == cigarop.op) {
			  cigarop.count += 1;
		  } else {
		  	cigar.push_back(cigarop);
		  	cigarop.op = dir[row][col];
		  	cigarop.count = 1;
		  }

		  if (dir[row][col] == ALN_OP_EQ || dir[row][col] == ALN_OP_X) {
		    row -= 1; col -= 1;
		  } else if (dir[row][col] == ALN_OP_I || dir[row][col] == ALN_OP_S) {
		    row -= 1;
      } else if (dir[row][col] == ALN_OP_D) {
        col -= 1;
		  }
	  }

	  if (cigarop.count > 0) {
	  	cigar.push_back(cigarop);
	  }

	  if (row > 0) { cigar.push_back(CigarOp(ALN_OP_I, row)); }
	  if (col > 0) { cigar.push_back(CigarOp(ALN_OP_D, col)); }

	  std::reverse(cigar.begin(), cigar.end());

	  return 0;
}

std::string CigarToString(const std::vector<CigarOp> &cigar) {
	std::stringstream ss;
	for (int32_t i=0; i<cigar.size(); i++) {
		ss << cigar[i].count << ALN_OP_TO_CHAR[cigar[i].op]	;
	}
	return ss.str();
}

std::string CigarToBasicString(const std::vector<CigarOp> &cigar) {
	std::stringstream ss;
	std::vector<CigarOp> cigar_basic(cigar);

	int32_t offset = 0;
	for (int32_t i=1; i<cigar_basic.size(); i++) {
		if (ALN_OP_TO_BASIC_CHAR[cigar_basic[i].op] == ALN_OP_TO_BASIC_CHAR[cigar_basic[i-offset-1].op]) {
			cigar_basic[i-offset-1].count += cigar_basic[i].count;
			offset += 1;
			continue;
		}
		cigar_basic[i-offset] = cigar_basic[i];
	}
	for (int32_t i=0; i<(cigar_basic.size()-offset); i++) {
		ss << cigar_basic[i].count << ALN_OP_TO_BASIC_CHAR[cigar_basic[i].op]	;
	}
	return ss.str();
}

void CigarToEdlibAln(const std::vector<CigarOp> &cigar, std::vector<int8_t>& alignment) {
	alignment.clear();
	int32_t aln_len = 0;
	for (int32_t i=0; i<cigar.size(); i++) {
		aln_len += cigar[i].count;
	}
	alignment.reserve(aln_len);
	for (int32_t i=0; i<cigar.size(); i++) {
		std::vector<int8_t> temp(cigar[i].count, ALN_OP_TO_EDLIB[cigar[i].op]);
		alignment.insert(alignment.end(), temp.begin(), temp.end());
	}
}

void CigarToAlignment(const char* q, int64_t ql, const char* t, int64_t tl,
							 const std::vector<CigarOp> &cigar, std::string &alnq, std::string &alnt, std::string &alnm) {
	std::stringstream ssq, sst, ssm;

	int32_t posq = 0, post = 0;

	for (int32_t i=0; i<cigar.size(); i++) {
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

	alnq = ssq.str();
	alnt = sst.str();
	alnm = ssm.str();
}

}
