/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cassert>
#include <cstdio>
#include <cmath>
#include <climits>

#include "bundle_base.h"
#include "as_pos32.hpp"
#include "as_pos.hpp"

bundle_base::bundle_base()
{
	tid = -1;
	chrm = "";
	lpos = 1 << 30;
	rpos = 0;
	strand = '.';
	is_allelic = false;
}

bundle_base::~bundle_base()
{}

int bundle_base::add_hit(const hit &ht)
{
	// store new hit
	hits.push_back(ht);
	// calcuate the boundaries on reference
	if(ht.pos < lpos) lpos = ht.pos;
	if(ht.rpos > rpos) rpos = ht.rpos;

	// try to include more paired-end reads
	int32_t p = ht.rpos;
	if(ht.mpos > ht.rpos && ht.mpos <= ht.rpos + 100000) p = ht.mpos;
	if(p > rpos) rpos = p;

	// set tid
	if(tid == -1) tid = ht.tid;
	assert(tid == ht.tid);

	// set strand
	if(hits.size() <= 1) strand = ht.strand;
	assert(strand == ht.strand);

	// set bundle is_allelic
	if (! is_allelic)
	{
		if (ht.has_variant()) is_allelic = true;
	}	
	return 0;
}

int bundle_base::buildbase()
{
	map<as_pos, int> apos_count;  // count of AS pos 

	// count apos
	for(int i = 0; i < hits.size(); i++)
	{
		hit &ht = hits[i];
		for(int j = 0; j < ht.apos.size(); j++)
		{
			as_pos a = ht.apos[j];
			if(apos_count.find(a) != apos_count.end())
			{
				apos_count.find(a)->second += 1;
			}
			else
			{
				apos_count.insert(make_pair(a, 1));
			}
		}
	}

	
	for(int i = 0; i < hits.size(); i++)
	{
		hit &ht = hits[i];
		// DEBUG
		/*
		if(strand != ht.strand)
		{
			printf("strand = %c, ht.strand = %c, ht.xs = %c,\n", strand, ht.strand, ht.xs);
		}
		*/

		for(int k = 0; k < ht.itv_align.size(); k++)
		{
			as_pos32 s (high32(ht.itv_align[k]).p32, "$");
			as_pos32 t (low32(ht.itv_align[k]).p32, "$");
			mmap += make_pair(ROI(s, t), 1);
		}

		for(int k = 0; k < ht.itvi.size(); k++)
		{
			as_pos32 s = high32(ht.itvi[k]);
			as_pos32 t = low32(ht.itvi[k]);
			assert(s.ale == "$");
			assert(t.ale == "$");
			imap += make_pair(ROI(s, t), 1);
		}

		for(int k = 0; k < ht.itvd.size(); k++)
		{
			as_pos32 s = high32(ht.itvd[k]);
			as_pos32 t = low32(ht.itvd[k]);
			assert(s.ale == "$");
			assert(t.ale == "$");
			imap += make_pair(ROI(s, t), 1);
		}
	}
	return 0;
}

bool bundle_base::overlap(const hit &ht) const
{
	if(mmap.find(ROI(ht.pos, ht.pos + 1)) != mmap.end()) return true;
	if(mmap.find(ROI(ht.rpos - 1, ht.rpos)) != mmap.end()) return true;
	return false;
}

int bundle_base::clear()
{
	is_allelic = false;
	tid = -1;
	chrm = "";
	lpos = 1 << 30;
	rpos = 0;
	strand = '.';
	apos_count.clear();
	hits.clear();
	mmap.clear();
	imap.clear();
	return 0;
}

