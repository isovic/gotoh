#ifndef TEST1_
#define TEST1_

#include "gtest/gtest.h"
#include "align.h"

TEST(AlignGlobal, CigarComparison1) {
  std::string q("GGATCGA");
  std::string t("GAATTCAGTTA");
  is::Penalties p(1, -1, -1, -1);
  is::AlignType aln_type = is::kGlobal;

  is::Align aln(q, t, p, aln_type, is::GlobalMargins());

  int32_t q_start = 0, q_end = 0, t_start = 0, t_end = 0, score = 0;
  std::vector<is::CigarOp> cigar;
  aln.GetAlignment(q_start, q_end, t_start, t_end, score, cigar);
  std::string cigar_string = CigarToString(cigar);

  EXPECT_EQ(cigar_string, std::string("1=1X1=1D2=1D1=2D1="));
}

TEST(AlignGlobal, CigarComparison2) {
  std::string q("ACGA");
  std::string t("ACCCGA");
  is::Penalties p(1, -1, -1, -1);
  is::AlignType aln_type = is::kGlobal;

  is::Align aln(q, t, p, aln_type, is::GlobalMargins());

  int32_t q_start = 0, q_end = 0, t_start = 0, t_end = 0, score = 0;
  std::vector<is::CigarOp> cigar;
  aln.GetAlignment(q_start, q_end, t_start, t_end, score, cigar);
  std::string cigar_string = CigarToString(cigar);

  EXPECT_EQ(cigar_string, std::string("1=2D3="));
}

TEST(AlignGlobal, CigarComparison3) {
  std::string q("A");
  std::string t("A");
  is::Penalties p(1, -1, -1, -1);
  is::AlignType aln_type = is::kGlobal;

  is::Align aln(q, t, p, aln_type, is::GlobalMargins());

  int32_t q_start = 0, q_end = 0, t_start = 0, t_end = 0, score = 0;
  std::vector<is::CigarOp> cigar;
  aln.GetAlignment(q_start, q_end, t_start, t_end, score, cigar);
  std::string cigar_string = CigarToString(cigar);

  EXPECT_EQ(cigar_string, std::string("1="));
}

TEST(AlignGlobal, CigarComparison4) {
  std::string q("A");
  std::string t("T");
  is::Penalties p(1, -1, -1, -1);
  is::AlignType aln_type = is::kGlobal;

  is::Align aln(q, t, p, aln_type, is::GlobalMargins());

  int32_t q_start = 0, q_end = 0, t_start = 0, t_end = 0, score = 0;
  std::vector<is::CigarOp> cigar;
  aln.GetAlignment(q_start, q_end, t_start, t_end, score, cigar);
  std::string cigar_string = CigarToString(cigar);

  EXPECT_EQ(cigar_string, std::string("1X"));
}

TEST(AlignGlobal, CigarComparison5) {
  std::string q("CCGA");
  std::string t("ACCCGA");
  is::Penalties p(1, -1, -1, -1);
  is::AlignType aln_type = is::kGlobal;

  is::Align aln(q, t, p, aln_type, is::GlobalMargins());

  int32_t q_start = 0, q_end = 0, t_start = 0, t_end = 0, score = 0;
  std::vector<is::CigarOp> cigar;
  aln.GetAlignment(q_start, q_end, t_start, t_end, score, cigar);
  std::string cigar_string = CigarToString(cigar);

  EXPECT_EQ(cigar_string, std::string("4="));
}

#endif
