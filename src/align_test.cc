/*
 * align_test.cc
 *
 *  Created on: Oct 12, 2016
 *      Author: isovic
 */

#include "align.h"
#include <iostream>
#include <string>
#include "../lib/gtest/gtest.h"

// TEST(AlignTest1, BasicTest1) {
void AlignTestGlobal1() {
  is::AlignType aln_type = is::kGlobal;

  std::string q1("GGATCGA"), t1("GAATTCAGTTA");
  is::Align aln1(q1, t1, is::Penalties(1, -1, -1, -1), aln_type, is::GlobalMargins());
  aln1.Verbose(q1, t1, std::cout);
  printf ("\n");
//  EXPECT_EQ(
//      "TTTACAGGATAGTGCCGCCAATCTTCCAGTGATACCCCGTGCCGCCAATCTTCCAGTATATACAGCACGA"
//      "GTAGC",
//      pc->Sequence);

  std::string q2("ACGA"), t2("ACCCGA");
  is::Align aln2(q2, t2, is::Penalties(1, -1, -1, -1), aln_type, is::GlobalMargins());
  aln2.Verbose(q2, t2, std::cout);
  printf ("\n");

  std::string q3("A"), t3("A");
  is::Align aln3(q3, t3, is::Penalties(1, -1, -1, -1), aln_type, is::GlobalMargins());
  aln3.Verbose(q3, t3, std::cout);
  printf ("\n");

  std::string q4("A"), t4("T");
  is::Align aln4(q4, t4, is::Penalties(1, -1, -1, -1), aln_type, is::GlobalMargins());
  aln4.Verbose(q4, t4, std::cout);
  printf ("\n");

  std::string q5("CCGA"), t5("ACCCGA");
  is::Align aln5(q5, t5, is::Penalties(1, -1, -1, -1), aln_type, is::GlobalMargins());
  aln5.Verbose(q5, t5, std::cout);
  printf ("\n");

//  EXPECT_EQ(1, 1);
}

void AlignTestLocal1() {
  is::AlignType aln_type = is::kLocal;

  std::string q1("GGATCGA"), t1("GAATTCAGTTA");
  is::Align aln1(q1, t1, is::Penalties(1, -1, -1, -1), aln_type, is::GlobalMargins());
  aln1.Verbose(q1, t1, std::cout);
  printf ("\n");
//  EXPECT_EQ(
//      "TTTACAGGATAGTGCCGCCAATCTTCCAGTGATACCCCGTGCCGCCAATCTTCCAGTATATACAGCACGA"
//      "GTAGC",
//      pc->Sequence);

  std::string q2("ACGA"), t2("ACCCGA");
  is::Align aln2(q2, t2, is::Penalties(1, -1, -1, -1), aln_type, is::GlobalMargins());
  aln2.Verbose(q2, t2, std::cout);
  printf ("\n");

  std::string q3("A"), t3("A");
  is::Align aln3(q3, t3, is::Penalties(1, -1, -1, -1), aln_type, is::GlobalMargins());
  aln3.Verbose(q3, t3, std::cout);
  printf ("\n");

  std::string q4("A"), t4("T");
  is::Align aln4(q4, t4, is::Penalties(1, -1, -1, -1), aln_type, is::GlobalMargins());
  aln4.Verbose(q4, t4, std::cout);
  printf ("\n");

  std::string q5("CCGA"), t5("ACCCGA");
  is::Align aln5(q5, t5, is::Penalties(1, -1, -1, -1), aln_type, is::GlobalMargins());
  aln5.Verbose(q5, t5, std::cout);
  printf ("\n");

//  EXPECT_EQ(1, 1);
}
