/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cstdio>
#include "junction.h"
#include "config.h"
#include "as_pos.hpp"
#include "as_pos32.hpp"

junction::junction()
{}

junction::junction(as_pos _p)
{
	lpos = high32(_p);
	rpos = low32(_p);
	count = 0;
	strand = '.';
	lexon = -1;
	rexon = -1;
	lregion = -1;
	rregion = -1;
}

junction::junction(as_pos _p, int _c)
{
	lpos = high32(_p);
	rpos = low32(_p);
	count = _c;
	strand = '.';
	lexon = -1;
	rexon = -1;
	lregion = -1;
	rregion = -1;
}

junction::junction(as_pos32 _p1, as_pos32 _p2, int _c)
{
	lpos = _p1;
	rpos = _p2;
	count = _c;
	strand = '.';
	lexon = -1;
	rexon = -1;
	lregion = -1;
	rregion = -1;
}

junction::junction(const junction &sp)
{
	lpos = sp.lpos;
	rpos = sp.rpos;
	count = sp.count;
	lexon = sp.lexon;
	rexon = sp.rexon;
	strand = sp.strand;
	lregion = sp.lregion;
	rregion = sp.rregion;
}

bool junction::operator<(const junction &j) const
{
	if (this->lpos < j.lpos) return true;
	if (this->lpos == j.lpos) 
	{	
		if (this->rpos < j.rpos) return true;
	}
	return false;
}

int junction::print(const string &chrm, int index) const
{
	// printf("junction %d: region = %s:%d%s-%d%s, %d -> %d, length = %d, count = %d, strand = %c\n", 
			// index, chrm.c_str(), lpos.p32, lpos.ale.c_str(), rpos.p32, rpos.ale.c_str(), lexon, rexon, rpos - lpos, count, strand);
	printf("junction %d: region = %s:%d%s-%d%s, region = %d -> %d, pexon = %d -> %d, length = %d, count = %d, strand = %c\n", 
			index, chrm.c_str(), lpos.p32, lpos.ale.c_str(), rpos.p32, rpos.ale.c_str(), lregion, rregion, lexon, rexon, rpos - lpos, count, strand);
	return 0;
}
