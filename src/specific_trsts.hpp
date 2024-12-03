/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Altai
(c) 2023 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef _SPECIFIC_TRSTS_
#define _SPECIFIC_TRSTS_

#include "gene.h"
#include "transcript.h"

namespace specific_trsts {
	int get_allele_spec_trsts(vector<transcript>& trsts1, vector<transcript>& trsts2,  double min_valid_cov);
	int hash_multi_exon_transcripts(const vector<transcript>& ts, vector<int>& hash_values);
	int get_multi_allele_spec_trsts(vector<transcript>& ts, const vector<int>& hasing_self, const set<int>& hasing_against, genotype gtspec);
	int get_single_exon_allele_spec_trsts(vector<transcript>& v1_single_exon, vector<transcript>& v2_single_exon);
    vector<transcript> intersection_of(const vector<transcript>& v1, const vector<transcript>& v2);
    vector<transcript> exclusive_of_1 (const vector<transcript>& v1, const vector<transcript>& v2);
	vector<transcript> union_of		  (const vector<transcript>& v1, const vector<transcript>& v2);
	vector<transcript> recover_full_from_partial_transcripts
		(const vector<transcript>& full_txs, const vector<transcript>& part_txs, double min_chain_overlap_ratio, bool will_change_gt);
}

#endif
