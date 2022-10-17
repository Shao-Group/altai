/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cstdio>
#include "partial_exon.h"
#include "util.h"
#include "as_pos32.hpp"
#include "vcf_data.h"

partial_exon::partial_exon(as_pos32 _lpos, as_pos32 _rpos, int _ltype, int _rtype, genotype _gt)
	: lpos(_lpos), rpos(_rpos), ltype(_ltype), rtype(_rtype), gt(_gt)
{
	type = 0;
	rid = -1;
	pid = -1;
}

int partial_exon::print(int index) const
{
	printf("partial_exon %d: [%d%s-%d%s), rid = %d, pid = %d, type = (%d, %d), length = %d, ave-abd = %.1lf, max-abd = %.1lf, std-abd = %.1lf, gt=%s\n",
			index, rid, pid, lpos.p32, lpos.ale.c_str(), rpos.p32, rpos.ale.c_str(), ltype, rtype, rpos - lpos, ave, max, dev, gt_str(gt));
	return 0;
}

bool partial_exon::is_allelic() const
{
	if(lpos.ale != "$")
	{
		assert(rpos.ale != "$");
		return true;
	}
	assert(rpos.ale == "$");
	return false;
}

int partial_exon::assign_as_cov(double _ave, double _max, double _dev)
{
	ave = _ave;
	max = _max;
	dev = _dev;
}

bool partial_exon::operator < (partial_exon pe)
{
	if (this->lpos.leftto(pe.lpos)) return true;
	if (this->lpos.samepos(pe.lpos))
	{
		if(this->rpos.leftto(pe.rpos))  return true;
		else if(this->rpos.rightto(pe.rpos))  return false;
		else return this->lpos < pe.lpos;
	}
	return false;
}

bool partial_exon::operator < (const partial_exon pe) const
{
	if (this->lpos.leftto(pe.lpos)) return true;
	if (this->lpos.samepos(pe.lpos))
	{
		if(this->rpos.leftto(pe.rpos))  return true;
		else if(this->rpos.rightto(pe.rpos))  return false;
		else return this->lpos < pe.lpos;
	}
	return false;
}
