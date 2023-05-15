/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include "vertex_info.h"
#include "as_pos32.hpp"
#include "config.h"

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
	gt = UNPHASED;
	as_type = 0;
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
	gt = UNPHASED;
	as_type = 0;
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
	gt = vi.gt;
	as_type = vi.as_type;
	if(DEBUG_MODE_ON) assert(as_type>= 0 && as_type <= 7);
	if(DEBUG_MODE_ON) assert(gt == ALLELE1 || gt == ALLELE2 || gt == UNPHASED || gt == NONSPECIFIC);
}

bool vertex_info::is_as_vertex() {return vertex_info::is_as_vertex(*this);}

bool vertex_info::is_adjacent_to_as_vertex() {return vertex_info::is_adjacent_to_as_vertex(*this);}

bool vertex_info::is_ordinary_vertex() {return vertex_info::is_ordinary_vertex(*this);}

bool vertex_info::is_as_vertex(vertex_info vi)
{
	if(DEBUG_MODE_ON) assert(vi.as_type>= 0 && vi.as_type <= 7);
	if (vi.as_type == AS_DIPLOIDVAR || vi.as_type == AS_MONOPLOIDVAR) return true;
	else return false;
}

bool vertex_info::is_adjacent_to_as_vertex(vertex_info vi)
{
	if(DEBUG_MODE_ON) assert(vi.as_type>= 0 && vi.as_type <= 7);
	if (vi.as_type == AJ_NSMONOVAR || vi.as_type == AJ_NONVAR) return true;
	else return false;
}

bool vertex_info::is_ordinary_vertex(vertex_info vi)
{
	if(DEBUG_MODE_ON) assert(vi.as_type>= 0 && vi.as_type <= 7);
	if (vi.as_type == NS_NONVAR || vi.as_type == NS_MONOVAR) return true;
	else return false;
}