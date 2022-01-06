/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __PARTIAL_EXON_H__
#define __PARTIAL_EXON_H__

#include <stdint.h>
#include <vector>
#include <string>
#include "as_pos32.hpp"

using namespace std;

class partial_exon
{
public:
	partial_exon(as_pos32 _lpos, as_pos32 _rpos, int _ltype, int _rtype);

public:
	as_pos32 lpos;					// the leftmost boundary on reference
	as_pos32 rpos;					// the rightmost boundary on reference
	int ltype;						// type of the left boundary
	int rtype;						// type of the right boundary
	bool is_allelic;				// whether it is allelic exon or not

	int rid;						// parental region id
	int pid;						// index in the parental pexons
	int type;						// label
	double ave;						// average abundance
	double max;						// maximum abundance
	double dev;						// standard-deviation of abundance
	bool operator < (partial_exon pe);
	bool operator < (const partial_exon pe) const;

public:
	string label() const;
	int print(int index) const;
};

#endif
