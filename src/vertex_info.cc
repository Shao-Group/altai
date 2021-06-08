/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include "vertex_info.h"
#include "as_pos32.hpp"

vertex_info::vertex_info()
{
	stddev = 1.0;
	length = 0;
	sdist = -1;
	tdist = -1;
	type = -1;
	lpos = as_pos32(0);
	rpos = as_pos32(0);
	pos = 0;
	lstrand = '.';
	rstrand = '.';
	regional = false;
}

vertex_info::vertex_info(int l)
	: length(l)
{
	stddev = 1.0;
	sdist = -1;
	tdist = -1;
	type = -1;
	lpos = as_pos32(0);
	rpos = as_pos32(0);
	pos = 0;
	lstrand = '.';
	rstrand = '.';
	regional = false;
}

vertex_info::vertex_info(const vertex_info &vi)
{
	stddev = vi.stddev;
	length = vi.length;
	sdist = vi.sdist;
	tdist = vi.tdist;
	type = vi.type;
	lpos = vi.lpos;
	rpos = vi.rpos;
	pos = vi.pos;
	lstrand = vi.lstrand;
	rstrand = vi.rstrand;
	regional = vi.regional;
}
