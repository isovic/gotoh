//============================================================================
// Name        : main.cc
// Author      : Ivan Sovic
// Version     :
// Copyright   : Copyright Ivan Sovic, 2016
// Description :
//============================================================================

#include <stdio.h>
#include <stdint.h>
#include <iostream>
#include <string>
#include "../lib/gtest/gtest.h"
#include "align.h"
#include "align_test.h"
using namespace std;

//G A A T T C A G T T A
//|   |   | |   |     |
//G G A _ T C _ G _ _ A
//
//G _ A A T T C A G T T A
//|     |   | |   |     |
//G G _ A _ T C _ G _ _ A

int main(int argc, char **argv) {
//  std::string target = "GAATTCAGTTA";
//  std::string query = "GGATCGA";
//
////  std::string target = "ACCCGA";
////  std::string query = "ACGA";
//
//  printf ("Aligning:\n  query: '%s'\n  target: '%s'\n", query.c_str(), target.c_str());
//
//  is::Align aln(query, target, is::Penalties(1, -1, -1, -1), is::kGlobal, is::GlobalMargins());

//  AlignTest1();
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
