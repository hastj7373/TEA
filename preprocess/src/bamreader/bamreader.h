#include <iostream>
#include <climits>
#include "api/BamAlignment.h"

void writeFQ(std::ofstream &os, BamTools::BamAlignment &al);
void writeFQpair(std::ofstream &f1, BamTools::BamAlignment &a,
                 std::ofstream &f2, BamTools::BamAlignment &b, int n_cutoff);
void split_read(BamTools::BamAlignment &al, int frag_size, int n_cutoff);
int getMateNumber(BamTools::BamAlignment &al);
bool isBigS(BamTools::BamAlignment &al, BamTools::CigarOp &cop, int n);
