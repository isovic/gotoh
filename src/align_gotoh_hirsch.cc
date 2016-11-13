/*
 * align.cpp
 *
 *  Created on: Oct 12, 2016
 *      Author: isovic
 */

#include "align_gotoh_hirsch.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <cmath>
#include "timer.hpp"

namespace is {

AlignGotohHirsch::AlignGotohHirsch(const std::string& q, const std::string& t, Penalties p,
             AlignType aln_type, GlobalMargins gm)
    : AlignGotohHirsch(q.c_str(), q.size(), t.c_str(), t.size(), p, aln_type, gm) {
}

AlignGotohHirsch::AlignGotohHirsch(const char* q, int64_t ql, const char* t, int64_t tl, Penalties p,
             AlignType aln_type, GlobalMargins gm)
              : Align(q, ql, t, tl, p, aln_type, gm) {
  Align_(q, ql, t, tl, p, gm, aln_type, q_start_, q_end_, t_start_, t_end_, score_, cigar_);
}

AlignGotohHirsch::~AlignGotohHirsch() {

}

int AlignGotohHirsch::Align_(const char* q, int64_t ql, const char* t, int64_t tl, Penalties p, GlobalMargins gm, AlignType aln_type,
                       int32_t &q_start, int32_t &q_end, int32_t &t_start, int32_t &t_end, int32_t &score, std::vector<is::CigarOp> &cigar) {

  if (aln_type == kGlobal) {
   return AlignGlobal_(q, ql, t, tl, p, gm, q_start_, q_end_, t_start_, t_end_, score_, cigar_);
  } else if (aln_type == kLocal) {
   return AlignLocal_(q, ql, t, tl, p);
  }

  // Unknown alignment type.
  return 1;
}

int AlignGotohHirsch::AlignGlobal_(const char* q, int64_t ql, const char* t, int64_t tl, Penalties p, GlobalMargins gm,
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
  if (gm.top) {
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
    if (gm.top) { H[0] = MINUS_INF; }
    // Penalize the first column.
    if (gm.left) { M[1][0] = i * p.gext; }

    for (int32_t j=1; j<(tl+1); j++) {
      V[1][j] = std::max((V[0][j] + p.gext), M[0][j] + p.gopen);
      H[1] = std::max(H[0] + p.gext, M[1][j-1] + p.gopen);

      int32_t diag = M[0][j-1] + ((q[i-1] == t[j-1]) ? p.match : p.mismatch);
      dir[i][j] = ALN_OP_M;  // Match/Mismatch. Needs to be replaced in the traceback by a correct op (faster this way).
      M[1][j] = diag;
      if (V[1][j] > M[1][j]) { M[1][j] = V[1][j]; dir[i][j] = ALN_OP_I; }  // Insertion.
      if (H[1] > M[1][j]) { M[1][j] = H[1]; dir[i][j] = ALN_OP_D; }  // Deletion.
    }
    V[0] = V[1];
    M[0] = M[1];
    H[0] = H[1];
  }

  // General case for global alignment. Start from the bottom right corner.
  bt_row = ql;
  bt_col = tl;
  bt_val = M[1][tl];

  // Handle semiglobal right margin. Allow the start position at any row of the last column.
  if (!gm.right) {
    for (int32_t i=0; i<(ql+1); i++) {
      if (M[i][tl] > bt_val) {
        bt_row = i;
        bt_val = M[i][tl];
      }
    }
  }

  // Handle semiglobal bottom margin. Allow the start position at any column of the last row.
  if (!gm.bottom) {
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

int AlignGotohHirsch::AlignLocal_(const char* q, int64_t ql, const char* t, int64_t tl,
                      Penalties p) {
  return 0;
}

int AlignGotohHirsch::Hirschberg_(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p, GlobalMargins gm) {

  // Storing only two lines for the main (M) and the vertical matrix (V).
  MatrixType M(2, std::vector<cell_t>(tl+1, 0));
  MatrixType V(2, std::vector<cell_t>(tl+1, 0));



  return 0;
}

// int AlignGotohHirsch::AlignBlock_(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p, GlobalMargins gm,
//                 MatrixType &M, MatrixType &V, const std::vector<cell_t> &H0_column) {
int AlignGotohHirsch::AlignBlock_(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p, GlobalMargins gm,
                std::vector<cell_t> &last_row) {

  // Storing only two lines for the main (M) and the vertical matrix (V).
  std::vector<std::vector<cell_t> > M(2, std::vector<cell_t>(tl+1, 0));
  std::vector<std::vector<cell_t> > V(2, std::vector<cell_t>(tl+1, MINUS_INF));
  // Only current and previous element of the horizontal matrix are needed.
  std::vector<int32_t> H(2, 0);

  for (int32_t i=0; i<(tl+1); i++) {
    M[0][i] = i * p.gext;
  }

  for (int32_t i=1; i<(ql+1); i++) {
//    cell_t[2] H = {H0, 0};
    M[1][0] = i * p.gext;
    V[1][0] = MINUS_INF;
    H[0] = MINUS_INF;

    for (int32_t j=1; j<(tl+1); j++) {
      V[1][j] = std::max((V[0][j] + p.gext), M[0][j] + p.gopen);
      H[1] = std::max(H[0] + p.gext, M[1][j-1] + p.gopen);

      int32_t diag = M[0][j-1] + ((q[i-1] == t[j-1]) ? p.match : p.mismatch);
      M[1][j] = diag; // Match/mismatch.
      if (V[1][j] > M[1][j]) { M[1][j] = V[1][j]; }  // Insertion.
      if (H[1] > M[1][j]) { M[1][j] = H[1]; }  // Deletion.

      H[0] = H[1];
    }
    V[0] = V[1];
    M[0] = M[1];
  }


//  cell_t H1 = 0;
//  cell_t[2] H = {H0, 0};

  // for (int32_t i=1; i<(ql+1); i++) {
  //   cell_t[2] H = {H0, 0};

  //   for (int32_t j=1; j<(tl+1); j++) {
  //     V[1][j] = std::max((V[0][j] + p.gext), M[0][j] + p.gopen);
  //     H[1] = std::max(H[0] + p.gext, M[1][j-1] + p.gopen);

  //     int32_t diag = M[0][j-1] + ((q[i-1] == t[j-1]) ? p.match : p.mismatch);
  //     M[1][j] = diag; // Match/mismatch.
  //     if (V[1][j] > M[1][j]) { M[1][j] = V[1][j]; }  // Insertion.
  //     if (H[1] > M[1][j]) { M[1][j] = H[1]; }  // Deletion.

  //     H[0] = H[1];
  //   }
  //   V[0] = V[1];
  //   M[0] = M[1];
  // }

  return 0;
}

}
