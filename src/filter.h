/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __FILTER_H__
#define __FILTER_H__

#include "gene.h"

class filter
{
public:
	filter(const vector<transcript> &v);

public:
	vector<transcript> trs;

public:
	// int simple_phase_set_by_coverage();
	// int simple_phase_set_by_variant_number();
	// int keep_as_transcripts_only();

	int join_single_exon_transcripts();
	int filter_length_coverage();
	int remove_nested_transcripts();
	int merge_single_exon_transcripts();
	int merge_single_exon_transcripts(vector<transcript> &trs0);
	int print();

private:
	bool join_transcripts();
	int locate_next_transcript(int t);
};

bool transcript_cmp(const transcript &x, const transcript &y);

#endif
