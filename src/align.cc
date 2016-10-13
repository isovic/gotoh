/*
 * align.cpp
 *
 *  Created on: Oct 12, 2016
 *      Author: isovic
 */

#include "align.h"

namespace is {

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

int Align::AlignGlobal(const char* q, int64_t ql, const char* t, int64_t tl,
                       Penalties p, GlobalMargins gm) {

  std::vector<std::vector<int32_t> > M(ql+1, std::vector<int32_t>(tl+1, 0));
  std::vector<std::vector<int32_t> > V(2, std::vector<int32_t>(tl+1, 0));
  std::vector<std::vector<int32_t> > H(2, std::vector<int32_t>(tl+1, 0));

  // Penalize the first column.
  if (gm.left) {
    for (int32_t i=0; i<ql; i++) { M[i][0] = i * p.gext; }
  }
  // Penalize the first row.
  if (gm.top) {
    for (int32_t i=0; i<tl; i++) {
      M[0][i] = i * p.gext;
      V[0][i] = i * p.gext;
    }
  }

  for (int32_t i=1; i<(ql+1); i++) {
    if (gm.top) { H[0][0] = i * p.gext; }

    for (int32_t j=1; j<(tl+1); j++) {
      V[1][j] = std::max((V[0][j] + p.gext), M[i-1][j] + p.gopen);
      H[1][j] = std::max(H[1][j-1] + p.gext, M[i][j-1] + p.gopen);
    }

    V[0] = V[1];
    H[0] = H[1];
  }

  return 0;
}

int Align::AlignLocal(const char* q, int64_t ql, const char* t, int64_t tl,
                      Penalties p) {
  return 0;
}

}
