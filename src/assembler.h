/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __ASSEMBLER_H__
#define __ASSEMBLER_H__

#include <fstream>
#include <string>
#include "bundle_base.h"
#include "bundle.h"
#include "transcript.h"
#include "splice_graph.h"
#include "gurobi_c++.h"

using namespace std;

class assembler
{
public:
	assembler();
	~assembler();

private:
	GRBEnv *env;
	samFile *sfn;
	bam_hdr_t *hdr;
	bam1_t *b1t;
	bundle_base bb1;		// +
	bundle_base bb2;		// -
	vector<bundle_base> pool;
	map<string, int> choose_map;

	int index;
	bool terminate;
	int qcnt;
	double qlen;
	vector<transcript> trsts;

public:
	int assemble();

private:
	int process(int n);
	int assemble(const splice_graph &gr, const hyper_set &hs, bool is_allelic);
	int assign_RPKM();
	int write();
	bool determine_regional_graph(splice_graph &gr);
	int read_prediction();
};

#endif
