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

const int32_t MINUS_INF = std::numeric_limits<int32_t>::min() + 1000000; // Allow a margin to skip overflow.

Align::Align(const std::string& q, const std::string& t, Penalties p,
             AlignType aln_type, GlobalMargins gm)
    : Align(q.c_str(), q.size(), t.c_str(), t.size(), p, aln_type, gm) {
}

Align::Align(const char* q, int64_t ql, const char* t, int64_t tl, Penalties p,
             AlignType aln_type, GlobalMargins gm) {
  if (aln_type == kGlobal) {
    AlignGlobal(q, ql, t, tl, p, gm);
  } else if (aln_type == kLocal) {
    AlignLocal(q, ql, t, tl, p);
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

int Align::AlignGlobal(const char* q, int64_t ql, const char* t, int64_t tl,
                       Penalties p, GlobalMargins gm) {

  std::vector<std::vector<int32_t> > M(ql+1, std::vector<int32_t>(tl+1, 0));
  std::vector<std::vector<int32_t> > dir(ql+1, std::vector<int32_t>(tl+1, 0));
  std::vector<std::vector<int32_t> > V(2, std::vector<int32_t>(tl+1, 0));
  std::vector<int32_t> H(tl+1, 0);

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

  for (int32_t i=1; i<(ql+1); i++) {
    if (gm.top) { H[0] = MINUS_INF; }

    for (int32_t j=1; j<(tl+1); j++) {
      V[1][j] = std::max((V[0][j] + p.gext), M[i-1][j] + p.gopen);
      H[j] = std::max(H[j-1] + p.gext, M[i][j-1] + p.gopen);

      int32_t diag = M[i-1][j-1] + ((q[i-1] == t[j-1]) ? p.match : p.mismatch);
      M[i][j] = diag; dir[i][j] = 0;  // Match/Mismatch.
      if (V[1][j] > M[i][j]) { M[i][j] = V[1][j]; dir[i][j] = 1; }  // Deletion.
      if (H[j] > M[i][j]) { M[i][j] = H[j]; dir[i][j] = 2; }  // Insertion.
    }
    V[0] = V[1];
  }

//  PrintMatrix(M);
//  printf ("\n");
//  PrintMatrix(dir);
//  printf ("Performing traceback.\n");

  // Traceback.
  std::stringstream ssq, sst, ssm;
  int32_t row = ql;
  int32_t col = tl;
  while (row > 0 && col > 0) {
    if (dir[row][col] == 0) {
      ssq << q[row-1];
      sst << t[col-1];
      ssm << ((q[row-1] == t[col-1]) ? "|" : "X");
      row -= 1; col -= 1;
    } else if (dir[row][col] == 2) {
      ssq << "-";
      sst << t[col-1];
      ssm << " ";
      col -= 1;
    } else if (dir[row][col] == 1) {
      ssq << q[row-1];
      sst << "-";
      ssm << " ";
      row -= 1;
    } else {
      fprintf (stderr, "ERROR: Unknown alignment move! Exiting.\n");
      exit(1);
    }
  }

  while (row > 0) {
    ssq << q[row-1];
    sst << "-";
    ssm << " ";
    row -= 1;
  }
  while (col > 0) {
    ssq << "-";
    sst << t[col-1];
    ssm << " ";
    col -= 1;
  }

  std::string alnq = ssq.str();
  std::string alnt = sst.str();
  std::string alnm = ssm.str();

  std::reverse(alnq.begin(), alnq.end());
  std::reverse(alnt.begin(), alnt.end());
  std::reverse(alnm.begin(), alnm.end());

  printf ("%s\n", alnt.c_str());
  printf ("%s\n", alnm.c_str());
  printf ("%s\n", alnq.c_str());

  return 0;
}

int Align::AlignLocal(const char* q, int64_t ql, const char* t, int64_t tl,
                      Penalties p) {
  return 0;
}

}
