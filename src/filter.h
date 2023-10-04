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

namespace specific_trsts {
    vector<transcript> intersection_of(const vector<transcript>& v1, const vector<transcript>& v2);
    vector<transcript> exclusive_of_1 (const vector<transcript>& v1, const vector<transcript>& v2);
	vector<transcript> union_of		  (const vector<transcript>& v1, const vector<transcript>& v2);
	vector<transcript> recover_full_from_partial_transcripts
		(const vector<transcript>& full_txs, const vector<transcript>& part_txs, double min_chain_overlap_ratio, bool will_change_gt);
}

class filter
{
public:
	filter(const vector<transcript> &v);

public:
	vector<transcript> trs;

public:
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
