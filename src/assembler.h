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
// #include "bundle.h"
#include "transcript.h"
#include "splice_graph.h"
#include "hyper_set.h"
#include "transcript_set.h"
#include "phaser.h"

using namespace std;

class assembler
{
public:
	assembler();
	~assembler();

private:
	samFile *sfn;
	bam_hdr_t *hdr;
	bam1_t *b1t;
	bundle_base bb1;		// +
	bundle_base bb2;		// -
	vector<bundle_base> pool;

	int hid;
	int index;
	bool terminate;
	int qcnt;
	double qlen;
	vector< vector<transcript> > trsts;
	vector< vector<transcript> > nonfull_trsts;
	vector< vector<transcript> > specific_full_trsts; // 0: nonspecific; 1: ALLELE1 spec; 2: ALLELE2 spec
	vector<transcript> recovered_allele1;
	vector<transcript> recovered_allele2;

public:
	int assemble();

private:
	int process(int n);
	int assemble(const splice_graph &gr, const hyper_set &hs, bool is_allelic, vector<transcript_set> &ts1, vector<transcript_set> &ts2);
	int assign_RPKM();
	int write();
	bool determine_regional_graph(splice_graph &gr);
};

#endif
