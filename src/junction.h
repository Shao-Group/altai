/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __BRIDGE_H__
#define __BRIDGE_H__

#include <stdint.h>
#include <string>
#include "as_pos.hpp"
#include "as_pos32.hpp"

using namespace std;

class junction
{
public:
	junction();
	junction(as_pos _p);
	junction(as_pos32, as_pos32, int);
	junction(as_pos _p, int _c);
	junction(const junction &p);

	

public:
	as_pos32 lpos;		// left position [left, right)
	as_pos32 rpos;		// right position
	int count;			// number of hits having this splice junction
	char strand;		// strandness of this junction

	int lexon;			// region index corresponds to lpos
	int rexon;			// region index corresponds to rpos
	
	bool operator<(const junction &x) const;

public:
	int print(const string &chrm, int index) const;
};

#endif
