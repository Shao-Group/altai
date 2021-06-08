/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __REGION_H__
#define __REGION_H__

#include <stdint.h>
#include <vector>
#include "interval_map.h"
#include "partial_exon.h"
#include "as_pos32.hpp"

using namespace std;

typedef pair<int, int> PI;

class region
{
public:
	region(as_pos32 _lpos, as_pos32 _rpos, int _ltype, int _rtype, const split_interval_map *_mmap, const split_interval_map *_imap);
	~region();

public:
	as_pos32 lpos;					// the leftmost boundary on reference
	as_pos32 rpos;					// the rightmost boundary on reference
	int ltype;						// type of the left boundary
	int rtype;						// type of the right boundary
	const split_interval_map *mmap;	// pointer to match interval map
	const split_interval_map *imap;	// pointer to indel interval map
	join_interval_map jmap;			// subregion intervals

	vector<partial_exon> pexons;	// generated partial exons

public:
	int print(int index) const;
	bool left_inclusive();
	bool right_inclusive();

private:
	int build_join_interval_map();
	int smooth_join_interval_map();
	bool empty_subregion(as_pos32 p1, as_pos32 p2);
	int build_partial_exons();
};

#endif
