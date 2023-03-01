/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __DRAW_H__
#define __DRAW_H__

#include <fstream>
#include "vcf_data.h"

using namespace std;

int draw_header(ofstream & fout);
int draw_footer(ofstream & fout);
int graphviz_header(ofstream & fout);
int graphviz_footer(ofstream & fout);
string graphviz_gt_color(genotype gt);

#endif
