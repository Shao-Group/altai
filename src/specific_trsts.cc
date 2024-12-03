/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Altai
(c) 2023 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include "config.h"
#include "as_pos32.hpp"
#include "specific_trsts.hpp"
#include <cassert>
#include <algorithm>

/*
** assign spec in place
** does not add new spec trsts, only assign NONSPECIFICity to ALLELE_ transcripts
*/
int specific_trsts::get_allele_spec_trsts(vector<transcript>& v1, vector<transcript>& v2, double min_valid_cov)
{
    // get hash values, multi trsts, single-exon transcripts
    vector<int> v1_hasing;
    vector<transcript> v1_multi;
    vector<transcript> v1_single_exon;
    int v1_nonvalid_counter = 0;
	for(const transcript& t1: v1)
	{
        assert (t1.coverage > 0);
        if (t1.coverage < min_valid_cov) {v1_nonvalid_counter ++ ; continue;}
		if (t1.exons.size() <= 1) v1_single_exon.push_back(t1);
        else v1_multi.push_back(t1);
	}
    specific_trsts::hash_multi_exon_transcripts(v1_multi, v1_hasing);
    assert(v1_multi.size() == v1_hasing.size());

    vector<transcript> v2_multi;
    vector<transcript> v2_single_exon;
    vector<int> v2_hasing;
    int v2_nonvalid_counter = 0;
    for(const transcript& t2: v2)
	{
        assert (t2.coverage > 0);
        if (t2.coverage < min_valid_cov) {v2_nonvalid_counter ++ ; continue;}
		if (t2.exons.size() <= 1) v2_single_exon.push_back(t2);
        else v2_multi.push_back(t2);
	}
    specific_trsts::hash_multi_exon_transcripts(v2_multi, v2_hasing);
    assert(v2_multi.size() == v2_hasing.size());

    // assign NONSPECIFICity to multi-exon transcripts in-place
    // FIXME:  get spec single exons
    specific_trsts::get_multi_allele_spec_trsts(v1_multi, v1_hasing, set<int>(v2_hasing.begin(), v2_hasing.end()), ALLELE1);
    specific_trsts::get_multi_allele_spec_trsts(v2_multi, v2_hasing, set<int>(v1_hasing.begin(), v1_hasing.end()), ALLELE2);
    specific_trsts::get_single_exon_allele_spec_trsts(v1_single_exon, v2_single_exon);
    

    // return in place
    int v1_before_size = v1.size();
    v1.clear();
    for(transcript t1: v1_multi) v1.push_back(t1);
    for(transcript t1: v1_single_exon) v1.push_back(t1);
    assert(v1_before_size == v1_nonvalid_counter + v1.size());

    int v2_before_size = v2.size();
    v2.clear();
    for(transcript t2: v2_multi) v2.push_back(t2);
    for(transcript t2: v2_single_exon) v2.push_back(t2);
    assert(v2_before_size == v2_nonvalid_counter + v2.size());
    
    return 0;
}

/*
** assign multi-exon spec in place
** does not add new spec trsts, only assign NONSPECIFICity to ALLELE_ transcripts
*/
int specific_trsts::get_multi_allele_spec_trsts(vector<transcript>& ts, const vector<int>& hasing_self, const set<int>& hasing_against, genotype gtspec)
{
    assert(ts.size() == hasing_self.size());
    assert(gtspec == ALLELE1 || gtspec == ALLELE2);
	for(int i = 0; i < ts.size(); i++)
	{    
        transcript& t = ts[i];
		assert(t.exons.size() >= 2);
        assert(!gt_conflict(t.gt, gtspec));
        if (t.gt == NONSPECIFIC) continue;
		if (hasing_against.find(hasing_self[i]) == hasing_against.end()) t.assign_gt(gtspec);
        else t.assign_gt(NONSPECIFIC);
	}
    return 0;
}   

int specific_trsts::get_single_exon_allele_spec_trsts(vector<transcript>& v1_single_exon, vector<transcript>& v2_single_exon)
{
    //FIXME:
    cerr << "specific_trsts::get_single_exon_allele_spec_trsts not implemented yet" << endl;
    return 0;
}

/*
** calculate hash values of `ts` and populate `hash_values`
*/
int specific_trsts::hash_multi_exon_transcripts(const vector<transcript>& ts, vector<int>& hash_values)
{
	for(const transcript& t: ts)
	{
		if(t.exons.size() <= 1) throw runtime_error("hash_multi_exon_transcripts cannot hash single-exon transcripts");

		size_t hashing = t.get_intron_chain_hashing();
		hash_values.push_back(hashing);
	}
    assert(hash_values.size() == ts.size());
    return 0;
}

// return MULTI-exon trsts in intersection of v1, v2
vector<transcript> specific_trsts::intersection_of(const vector<transcript>& v1, const vector<transcript>& v2)
{
	vector<transcript> v0;
	if(v1.size() == 0 || v2.size() == 0) return v0;

	set<int> v2_hasing;
	for(const transcript& t2: v2)
	{
		if (t2.exons.size() <= 1) continue;
		size_t t2_hashing = t2.get_intron_chain_hashing();
		v2_hasing.insert(t2_hashing);
	}

	for(const transcript& t1: v1)
	{
		if (t1.exons.size() <= 1) continue;
		size_t t1_hashing = t1.get_intron_chain_hashing();
		if (v2_hasing.find(t1_hashing) == v2_hasing.end()) continue;
		v0.push_back(t1);
	}

	return v0;
}

// return MULTI-exon trsts present in only v1 
vector<transcript> specific_trsts::exclusive_of_1(const vector<transcript>& v1, const vector<transcript>& v2)
{
	vector<transcript> v0;
	if(v1.size() == 0) return v0;

	set<int> v2_hasing;
	for(const transcript& t2: v2)
	{
		if (t2.exons.size() <= 1) continue;
		size_t t2_hashing = t2.get_intron_chain_hashing();
		v2_hasing.insert(t2_hashing);
	}

	for(const transcript& t1: v1)
	{
		if (t1.exons.size() <= 1) continue;
		size_t t1_hashing = t1.get_intron_chain_hashing();
		if (v2_hasing.find(t1_hashing) != v2_hasing.end()) continue;
		v0.push_back(t1);
	}

	return v0;
}

// return MULTI-exon trsts in unions of v1, v2
// TODO: mark which allele the trst is from
vector<transcript> specific_trsts::union_of(const vector<transcript>& v1, const vector<transcript>& v2)
{
	vector<transcript> v0;
	if(v1.size() == 0 && v2.size() == 0) return v0;

	set<int> hashing;
	for(const transcript& t1: v1)
	{
		if (t1.exons.size() <= 1) continue;
		size_t t1_hashing = t1.get_intron_chain_hashing();
		if (hashing.find(t1_hashing) != hashing.end()) continue;
		
		v0.push_back(t1);
		hashing.insert(t1_hashing);
	}


	for(const transcript& t2: v2)
	{
		if (t2.exons.size() <= 1) continue;
		size_t t2_hashing = t2.get_intron_chain_hashing();
		if (hashing.find(t2_hashing) != hashing.end()) continue;

		v0.push_back(t2);
		hashing.insert(t2_hashing);
	}

	return v0;
}

/*
*	recovers full transcripts from partial transcripts, if >50% match
*	only works for multi-exon transcripts
*	@return a copy of recovered transcripts
*/

/*
vector<transcript> specific_trsts::recover_full_from_partial_transcripts
	(const vector<transcript>& full_txs, const vector<transcript>& part_txs, double min_chain_overlap_ratio, bool will_change_gt)
{
	vector<transcript> recovered;
	set<int> recovered_hasing;

	if (min_chain_overlap_ratio <= 0 || full_txs.size() == 0 || part_txs.size() == 0) return recovered;
	assert(min_chain_overlap_ratio > 0 && min_chain_overlap_ratio <= 1);

	for(const transcript& fullt: full_txs)
	{
		if (fullt.exons.size() <= 1) continue;
		size_t fullt_hashing = fullt.get_intron_chain_hashing();
		if (recovered_hasing.find(fullt_hashing) != recovered_hasing.end()) continue;
		
		const vector<PI32>& fullt_chain = fullt.get_intron_chain();

		for(const transcript& partt: part_txs)
		{
			if (recovered_hasing.find(fullt_hashing) != recovered_hasing.end()) continue;
			const vector<PI32>& part_chain = partt.get_intron_chain();

			vector<PI32> intersected;
			set_intersection(fullt_chain.begin(), fullt_chain.end(), part_chain.begin(), part_chain.end(), back_inserter(intersected));
			if (intersected.size() >= fullt_chain.size() * min_chain_overlap_ratio)
			{
				transcript t(fullt);
				t.coverage = partt.coverage;

				if (will_change_gt) t.transform_gt(partt.gt);
				else assert(!gt_conflict(fullt.gt, partt.gt));
				
				recovered.push_back(t);
				recovered_hasing.insert(fullt_hashing);
				break;
			}
		}
	}

	assert(recovered.size() == recovered_hasing.size());
	return recovered;
}
*/