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

namespace is {

typedef enum {
  kGlobal,
  kLocal
} AlignType;

class Penalties {
 public:
  Penalties(int32_t _match, int32_t _mismatch, int32_t _gopen, int32_t _gext)
      : match(_match),
        mismatch(_mismatch),
        gopen(_gopen),
        gext(_gext) {
  }
  ;
  int32_t match, mismatch, gopen, gext;
};

class GlobalMargins {
 public:
  GlobalMargins()
      : top(false),
        left(false),
        bottom(false),
        right(false) {
  }
  GlobalMargins(bool _top, bool _left, bool _bottom, bool _right)
      : top(_top),
        left(_left),
        bottom(_bottom),
        right(_right) {
  }
  bool top, left, bottom, right;
};

class Align {
 public:
  Align(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p,
        AlignType aln_type, GlobalMargins gm);
  Align(const std::string &q, const std::string &t, Penalties p,
        AlignType aln_type, GlobalMargins gm);
  ~Align();

 private:
  int AlignGlobal(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p, GlobalMargins gm);
  int AlignLocal(const char *q, int64_t ql, const char *t, int64_t tl, Penalties p);

};

}

#endif /* ALIGN_H_ */
