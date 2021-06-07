/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __GTF_TRANSCRIPT_H__
#define __GTF_TRANSCRIPT_H__

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "item.h"
#include "htslib/faidx.h"
#include "../../src/as_pos32.hpp"

using namespace std;

typedef pair<as_pos32, as_pos32> PI32;

class transcript
{
public:
	transcript(const item &ie);
	transcript();
	~transcript();

public:
	bool operator< (const transcript &t) const;

public:
	string seqname;
	string source;
	string feature;
	string gene_id;
	string transcript_id;
	string gene_type;
	string transcript_type;
	as_pos32 start;
	as_pos32 end;
	double score;
	char strand;
	int frame;
	double coverage;
	double covratio;
	double RPKM;
	double FPKM;
	double TPM;

	vector<PI32> exons;

public:
	int add_exon(as_pos32 s, as_pos32 t);
	int add_exon(const item &e);
	int assign_RPKM(double factor);
	int sort();
	int clear();
	int shrink();
	int assign(const item &e);
	int length() const;
	PI32 get_bounds() const;
	PI32 get_first_intron() const;
	vector<PI32> get_intron_chain() const;
	bool intron_chain_match(const transcript &t) const;
	string label() const;
	int write(ostream &fout) const;
	int write_gvf(ostream &fout) const;
	int write_fasta(ostream &fout, int line_len, faidx_t *fai) const;
	int write_fasta_AS_only(ostream &fout, int line_len, faidx_t *fai) const;

	static int reverse_complement_DNA(string &, const string);
};

#endif
