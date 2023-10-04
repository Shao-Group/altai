/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include "filter.h"
#include "config.h"
#include <cassert>
#include <algorithm>
#include "as_pos32.hpp"

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

filter::filter(const vector<transcript> &v)
	:trs(v)
{}

int filter::filter_length_coverage()
{
	// calculate relative ratio (only for multi-exon transcripts)
	double sum_cov = 0;
	vector<double> ratio;
	for(int i = 0; i < trs.size(); i++)
	{
		int e = trs[i].exons.size();
		if(e <= 1) continue;
		sum_cov += trs[i].coverage;
	}
	for(int i = 0; i < trs.size(); i++)
	{
		int e = trs[i].exons.size();
		if(e <= 1) ratio.push_back(1.0);
		else ratio.push_back(trs[i].coverage / sum_cov);
	}
	
	vector<transcript> v;
	for(int i = 0; i < trs.size(); i++)
	{
		int e = trs[i].exons.size();
		int minl = min_transcript_length_base + e * min_transcript_length_increase;
		if(trs[i].length() < minl) continue;
		if(e == 1 && trs[i].coverage < min_single_exon_coverage) continue;
		if(e >= 2 && trs[i].coverage < min_transcript_coverage) continue;
		v.push_back(trs[i]);
	}
	trs = v;
	return 0;
}

int filter::remove_nested_transcripts()
{
	set<int> s;
	for(int i = 0; i < trs.size(); i++)
	{
		vector<PI32> v = trs[i].exons;
		if(v.size() <= 1) continue;
		double w1 = trs[i].coverage;
		bool b = false;
		for(int k = 1; k < v.size(); k++)
		{
			as_pos32 p = v[k - 1].second;
			as_pos32 q = v[k - 0].first;

			for(int j = 0; j < trs.size(); j++)
			{
				if(trs[j].exons.size() <= 1) continue;
				PI32 pq = trs[j].get_bounds();
				double w2 = trs[j].coverage;

				if(w2 >= w1 && pq.first.rightto(p) && pq.second.leftto(q))
				{
					b = true;
					break;
				}
			}
			if(b == true) break;
		}
		if(b == true) s.insert(i);
	}

	vector<transcript> v;
	for(int i = 0; i < trs.size(); i++)
	{
		if(s.find(i) != s.end()) continue;
		v.push_back(trs[i]);
	}

	trs = v;
	return 0;
}

int filter::join_single_exon_transcripts()
{
	while(true)
	{
		bool b = join_transcripts();
		if(b == false) break;
	}
	return 0;
}

bool filter::join_transcripts()
{
	sort(trs.begin(), trs.end(), transcript_cmp);
	//print();

	int32_t mind = min_bundle_gap;
	int ki = -1, kj = -1;
	for(int i = 0; i < trs.size(); i++)
	{
		int j = locate_next_transcript(i);
		if(j == -1) continue;
		if(trs[i].exons.size() >= 2 && trs[j].exons.size() >= 2) continue;
		int32_t d = trs[j].get_bounds().first.p32 - trs[i].get_bounds().second.p32;
		if(d > mind) continue;
		mind = d;
		ki = i;
		kj = j;
	}
	if(ki == -1 || kj == -1) return false;
	if(mind > min_bundle_gap - 1) return false;

	if(verbose >= 3) printf("join transcript %d and %d\n", ki, kj);

	if(trs[ki].exons.size() >= 2)
	{
		assert(trs[kj].exons.size() == 1);
		as_pos32 p1 = trs[ki].get_bounds().second;
		as_pos32 p2 = trs[kj].get_bounds().second;
		trs[ki].add_exon(p1, p2);
		trs[kj].sort();
		trs[ki].shrink();
		trs.erase(trs.begin() + kj);
		return true;
	}
	else if(trs[kj].exons.size() >= 2)
	{
		assert(trs[ki].exons.size() == 1);
		as_pos32 p1 = trs[ki].get_bounds().first;
		as_pos32 p2 = trs[kj].get_bounds().first;
		trs[kj].add_exon(p1, p2);
		trs[kj].sort();
		trs[kj].shrink();
		trs.erase(trs.begin() + ki);
		return true;
	}
	else
	{
		assert(trs[ki].exons.size() == 1);
		assert(trs[kj].exons.size() == 1);
		as_pos32 p1 = trs[ki].get_bounds().first;
		as_pos32 p2 = trs[kj].get_bounds().first;
		trs[kj].add_exon(p1, p2);
		trs[kj].sort();
		trs[kj].shrink();
		double cov = 0;
		cov += trs[ki].coverage * trs[ki].length();
		cov += trs[kj].coverage * trs[kj].length();
		cov /= (trs[ki].length() + trs[kj].length());
		trs[kj].coverage = cov;
		trs.erase(trs.begin() + ki);
		return true;
	}

	return true;
}

int filter::locate_next_transcript(int t)
{
	if(t < 0 || t >= trs.size()) return -1;
	PI32 p = trs[t].get_bounds();
	int a = 0;
	int b = trs.size() - 1;
	if(trs[b].get_bounds().first.leftto(p.second)) return -1;
	while(true)
	{
		assert(a <= b);
		if(a == b) return a;
		int k = (a + b) / 2;
		if(trs[k].get_bounds().first.samepos(p.second)) return k;
		if(trs[k].get_bounds().first.leftsameto(p.second)) a = k + 1;
		if(trs[k].get_bounds().first.rightto(p.second)) b = k;
	}
	assert(false);
	return -1;
}

int filter::merge_single_exon_transcripts(vector<transcript> &trs0)
{
	// must be the same gt in order to merge
	genotype g = UNPHASED;
	for(const transcript & t: trs0)
	{
		if(t.gt == ALLELE1 || t.gt == ALLELE2)
		{
			assert(! gt_conflict(g, t.gt));
			g = t.gt;
		}
	}

	typedef pair<PI32, int> PPI;
	vector<PPI> vv;
	for(int i = 0; i < trs0.size(); i++)
	{
		vector<PI32> v = trs0[i].exons;
		for(int k = 0; k < v.size(); k++)
		{
			vv.push_back(PPI(v[k], i));
		}
	}

	sort(vv.begin(), vv.end());

	set<int> fb;
	for(int i = 0; i < vv.size(); i++)
	{
		as_pos32 p1 = vv[i].first.first;
		as_pos32 q1 = vv[i].first.second;
		int k1 = vv[i].second;
		transcript &t1 = trs0[k1];
		if(t1.exons.size() != 1) continue;
		if(t1.strand != '.') continue;

		bool b = false;
		for(int k = i - 1; k >= 0 && k >= i - 10; k--)
		{
			as_pos32 p2 = vv[k].first.first;
			as_pos32 q2 = vv[k].first.second;
			int k2 = vv[k].second;
			if(fb.find(k2) != fb.end()) continue;
			transcript &t2 = trs0[k2];
			if(t2.seqname != t1.seqname) continue;

			assert(p1.rightsameto(p2));
			if(q2.leftto(q1)) continue;

			//if(b == true) printf("AAA insert k1 = %d (%d, %d) to fb with k2 = %d (%d, %d)\n", k1, p1, q1, k2, p2, q2);

			b = true;
			break;
		}

		if(b == true) fb.insert(k1);
		if(b == true) continue;

		for(int k = i + 1; k < vv.size(); k++)
		{
			as_pos32 p2 = vv[k].first.first;
			as_pos32 q2 = vv[k].first.second;
			int k2 = vv[k].second;
			if(fb.find(k2) != fb.end()) continue;
			transcript &t2 = trs0[k2];
			if(t2.seqname != t1.seqname) continue;

			if(p2.rightto(p1)) break;
			assert(p2 == p1);
			if(q2.leftto(q1)) continue;
			b = true;

			//if(b == true) printf("BBB insert k1 = %d (%d, %d) to fb with k2 = %d (%d, %d)\n", k1, p1, q1, k2, p2, q2);

			break;
		}
		if(b == true) fb.insert(k1);
	}

	vector<transcript> v;
	for(int i = 0; i < trs0.size(); i++)
	{
		if(fb.find(i) != fb.end()) continue;
		v.push_back(trs0[i]);
	}
	trs0 = v;
	return 0;
}

int filter::merge_single_exon_transcripts()
{
	typedef vector<transcript> VT;
	typedef pair<string, VT> PSVT;
	typedef map<string, VT> MSVT;

	MSVT msvt;
	for(int i = 0; i < trs.size(); i++)
	{
		transcript &t = trs[i];
		string s = t.seqname;
		if(msvt.find(s) == msvt.end())
		{
			VT v;
			v.push_back(t);
			msvt.insert(PSVT(s, v));
		}
		else
		{
			msvt[s].push_back(t);
		}
	}

	trs.clear();
	for(MSVT::iterator it = msvt.begin(); it != msvt.end(); it++)
	{
		merge_single_exon_transcripts(it->second);
		trs.insert(trs.end(), it->second.begin(), it->second.end());
	}

	return 0;
}

int filter::print()
{
	for(int i = 0; i < trs.size(); i++)
	{
		transcript &t = trs[i];
		printf("transcript %d: exons = %lu, pos = %d-%d\n",
				i, t.exons.size(), t.get_bounds().first.p32, t.get_bounds().second.p32);
	}
	return 0;
}

bool transcript_cmp(const transcript &x, const transcript &y)
{
	if(x.exons[0].first < y.exons[0].first) return true;
	else return false;
}
