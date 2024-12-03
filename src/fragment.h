/*
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __FRAGMENT_H__
#define __FRAGMENT_H__

#include <vector>
#include <stdint.h>

#include "hit.h"
#include "path.h"

using namespace std;

class fragment
{
public:
	fragment(hit *x1, hit *x2);

public:
	hit* h1;			// list of first mate
	hit* h2;			// list of second mate
	int cnt;			// count of the equal hits
	as_pos32 lpos;		// equals to hits[k1].pos
	as_pos32 rpos;		// equals to hits[k2].rpos
	int32_t k1l;		// k1-left-outside length
	int32_t k1r;		// k1-right-outside length
	int32_t k2l;		// k2-left-outside length
	int32_t k2r;		// k2-right-outside length
	bool b1;			// whether left mate can be shorten
	bool b2;			// whether right mate can be shorten
	vector<path> paths;	// possible connecting paths

	genotype gt;
	int type;		// 0: paird-end, 1: UMI-paird 2: both
	int next;		// index for next fragments in UMI-linked

public:
	// bool equal(const fragment &f) const;
	int append(const fragment &f);
	int print(int index);
	// int set_paired(bool b);
	int set_bridged(bool b);
	int clear();
};

// bool compare_fragment(const fragment &f1, const fragment &f2);

#endif
