/*
 * align_test.cc
 *
 *  Created on: Oct 12, 2016
 *      Author: isovic
 */

#include "align.h"
#include <string>
#include "../lib/gtest/gtest.h"

// TEST(AlignTest1, BasicTest1) {
void AlignTest1() {
  is::Align aln1(std::string("GGATCGA"), std::string("GAATTCAGTTA"), is::Penalties(1, -1, -1, -1), is::kGlobal, is::GlobalMargins());
  printf ("\n");
//  EXPECT_EQ(
//      "TTTACAGGATAGTGCCGCCAATCTTCCAGTGATACCCCGTGCCGCCAATCTTCCAGTATATACAGCACGA"
//      "GTAGC",
//      pc->Sequence);

  is::Align aln2(std::string("ACGA"), std::string("ACCCGA"), is::Penalties(1, -1, -1, -1), is::kGlobal, is::GlobalMargins());
  printf ("\n");
  is::Align aln3(std::string("A"), std::string("A"), is::Penalties(1, -1, -1, -1), is::kGlobal, is::GlobalMargins());
  printf ("\n");
  is::Align aln4(std::string("A"), std::string("T"), is::Penalties(1, -1, -1, -1), is::kGlobal, is::GlobalMargins());
  printf ("\n");

//  EXPECT_EQ(1, 1);
}
