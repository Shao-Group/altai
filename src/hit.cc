/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cstring>
#include <cassert>
#include <cstdio>
#include <sstream>
#include <cmath>
#include <map>
#include "hit.h"
#include "util.h"
#include "config.h"
#include "as_pos.hpp"
#include "vcf_data.h"

/*
hit::hit(int32_t p)
{
	bam1_core_t::pos = p;
	strand = '.';
	xs = '.';
	ts = '.';
	hi = -1;
	nh = -1;
	nm = 0;
	qlen = 0;
	cigar = NULL;
}
*/

hit& hit::operator=(const hit &h)
{
	bam1_core_t::operator=(h);
	hid = h.hid;
	rpos = h.rpos;
	qlen = h.qlen;
	qname = h.qname;
	strand = h.strand;
	spos = h.spos;
	xs = h.xs;
	ts = h.ts;
	hi = h.hi;
	nh = h.nh;
	nm = h.nm;
	itvm = h.itvm;
	itvi = h.itvi;
	itvd = h.itvd;
	apos = h.apos;
	itvna = h.itvna;
	chrm = h.chrm;
	
	vlist = h.vlist;
	paired = h.paired;
	bridged = h.bridged;
	qhash = h.qhash;
	next = h.next;

	umi = h.umi;
	return *this;
}

hit::hit(const hit &h)
	:bam1_core_t(h)
{
	rpos = h.rpos;
	qlen = h.qlen;
	qname = h.qname;
	strand = h.strand;
	spos = h.spos;
	xs = h.xs;
	ts = h.ts;
	hi = h.hi;
	nh = h.nh;
	nm = h.nm;
	itvm = h.itvm;
	itvi = h.itvi;
	itvd = h.itvd;
	apos = h.apos;
	itvna = h.itvna;
	chrm = h.chrm;

	vlist = h.vlist;
	paired = h.paired;
	bridged = h.bridged;
	qhash = h.qhash;
	next = h.next;

	umi = h.umi;
}

hit::~hit()
{
}


// hit::hit(bam1_t *b, std::string chrm_name)
hit::hit(bam1_t *b, std::string chrm_name, int id) 
	:bam1_core_t(b->core), hid(id)
{
	chrm = chrm_name;
	// fetch query name
	qname = get_qname(b);
	qhash = string_hash(qname);
	paired = false;
	bridged = false;
	next = NULL;

	// compute rpos
	rpos = pos + (int32_t)bam_cigar2rlen(n_cigar, bam_get_cigar(b));
	qlen = (int32_t)bam_cigar2qlen(n_cigar, bam_get_cigar(b));

	// get cigar
	assert(n_cigar <= max_num_cigar);
	assert(n_cigar >= 1);
	uint32_t * cigar = bam_get_cigar(b);

	// build spos and apos
	spos.clear();
	apos.clear();
	int32_t p = pos;
	int32_t q = 0;
    for(int k = 0; k < n_cigar; k++)
	{
		if (bam_cigar_type(bam_cigar_op(cigar[k]))&2)							// p = consumes reference, excluded
			p += bam_cigar_oplen(cigar[k]);

		if (bam_cigar_type(bam_cigar_op(cigar[k]))&1)							// q = consumes query
			q += bam_cigar_oplen(cigar[k]);

		// build spos
		if(bam_cigar_op(cigar[k]) == BAM_CREF_SKIP)
		{
			if(bam_cigar_op(cigar[k-1]) != BAM_CMATCH) continue;
			if(bam_cigar_op(cigar[k+1]) != BAM_CMATCH) continue;
			// consider ALL splice positions
			// if(bam_cigar_oplen(cigar[k-1]) < min_flank_length) continue;
			// if(bam_cigar_oplen(cigar[k+1]) < min_flank_length) continue;
			int32_t s = p - bam_cigar_oplen(cigar[k]);							// s = previous consumed ref, included
			int64_t p64 = pack(s, p);
			spos.push_back(as_pos(p64, "$"));
		}		
		//build apos
		else if (bam_cigar_op(cigar[k]) == BAM_CMATCH)
		{
			int32_t s = p - bam_cigar_oplen(cigar[k]);
			int32_t itvm_st = s;
			int32_t itvm_ed = p;
						
			if (vcf_file == "") continue;
			if (vmap_chrm != chrm)
			{
				try
				{
					vcf_map_it = vcf_map.at(chrm).begin();
					vcf_map_end = vcf_map.at(chrm).end();
					vcf_map_len_it =  vcf_map_len.at(chrm).begin();
					vcf_map_len_end = vcf_map_len.at(chrm).end();
					vmap_chrm = chrm;
				}
				catch (std::out_of_range)
				{
					vmap_chrm = "";
					continue;
				}
			}

			// AS related
			auto it = vcf_map_it;
			auto it_len = vcf_map_len_it;
			for ( ; it != vcf_map_end && it_len!= vcf_map_len_end; vcf_data::increse_it(it, it_len))					// iterate through vcf
			{	
				if (it->first >= p)  break;
				if (it->first <  s)  
				{
					if (it->first < pos) { vcf_map_it = it; vcf_map_len_it = it_len;}
					continue;
				}
				int32_t alelpos = it->first; 									// 0-based ref pos, included
				int32_t qpos = alelpos - p + q;									// 0-based query pos, included

				std::string ale = "*";
				uint8_t *seq_ptr = bam_get_seq (b);
				for (auto iit = (it->second).begin(); iit != (it->second).end(); iit++) 								// Note: INS also take in here
				{
					// iit is the ale sequence
					for (int i = 0; i < (*iit).size(); ++i)
					{
						if (bam_seqi(seq_ptr, qpos + i) ==1 && (*iit)[i] != 'A')  break;
						if (bam_seqi(seq_ptr, qpos + i) ==2 && (*iit)[i] != 'C')  break;
						if (bam_seqi(seq_ptr, qpos + i) ==4 && (*iit)[i] != 'G')  break;
						if (bam_seqi(seq_ptr, qpos + i) ==8 && (*iit)[i] != 'T')  break;

						if (i == (*iit).size() - 1)  ale = string(*iit);
					}
				}				
				int32_t alerpos = alelpos + it_len->second;						// 0-based query pos, excluded
				if (alelpos < s)  continue;
				if (alerpos > p)
					if (alerpos - p != it_len->second - 1 || k+1 >= n_cigar || bam_cigar_op(cigar[k+1]) != BAM_CDEL)
						continue;												// the later condition checks whether it is DEL 				
				apos.push_back(as_pos(pack(alelpos, alerpos), ale));
			}
		}
	}

	// open for scallop+coral
	itvi.clear();
	itvd.clear();
	p = pos;
    for(int k = 0; k < n_cigar; k++)
	{
		if (bam_cigar_type(bam_cigar_op(cigar[k]))&2)
		{
			p += bam_cigar_oplen(cigar[k]);
		}

		if(bam_cigar_op(cigar[k]) == BAM_CMATCH)
		{
			int32_t s = p - bam_cigar_oplen(cigar[k]);
			itvm.push_back(as_pos(pack(s, p), "$"));
		}

		if(bam_cigar_op(cigar[k]) == BAM_CINS)
		{
			itvi.push_back(as_pos(pack(p - 1, p + 1), "$"));
		}

		if(bam_cigar_op(cigar[k]) == BAM_CDEL)
		{
			int32_t s = p - bam_cigar_oplen(cigar[k]);
			itvd.push_back(as_pos(pack(s, p), "$"));
		}
	}

	if (apos.size() != 0) make_itvna();
	else itvna = itvm;
}


int hit::make_itvna()
{
	itvna.clear();
	if (apos.size() == 0) 
	{
		itvna = itvm;
		return 0;
	}
	auto it1 = itvm.begin();
	auto it2 = apos.begin();
	int32_t l = high32(*it1);
	int32_t r = low32(*it1);
	int32_t a = high32(*it2);
	int32_t b = low32(*it2);
	while (it2 != apos.end() && it1 != itvm.end())	
	{
		if (r <= a) 
		{
			if(l < r) itvna.push_back(as_pos(pack(l, r), "$"));
			++it1;
			if (it1 != itvm.end()) {l = high32(*it1); r = low32(*it1);}
		}
		else if (b <= l)
		{
			++it2;
			if (it2 != apos.end()) {a = high32(*it2); b = low32(*it2);}
		} 
		// now apos intersect itvm (a < r && b > l)
		else if (a <= l && b <= r)
		{
			l = b;
			++it2;
			if (it2 != apos.end()) {a = high32(*it2); b = low32(*it2);}
		}
		else if (a > l && b <= r)
		{
			itvna.push_back(as_pos(pack(l, a), "$"));
			l = b;
			++it2;
			if (it2 != apos.end()) {a = high32(*it2); b = low32(*it2);}
		}
		else if (a > l && b > r)
		{
			itvna.push_back(as_pos(pack(l, a), "$"));
			++it1;
			if (it1 != itvm.end()) {l = high32(*it1); r = low32(*it1);}
		}
		else if (a <= l && b > r)
		{
			if (DEBUG_MODE_ON) printf("a %d, b %d, l %d, r %d", a, b, l, r);
			++it1;
			if (it1 != itvm.end()) {l = high32(*it1); r = low32(*it1);}
		}
		else assert(0);
	}

	while (it1 != itvm.end())
	{
		itvna.push_back(as_pos(pack(l, r), "$"));
		++it1;
		if (it1 != itvm.end()) {l = high32(*it1); r = low32(*it1);}
	}
	itvna.push_back(as_pos(pack(l, r), "$"));

	return 0;
}

int hit::get_aligned_intervals(vector<as_pos> &v) const
{
	v.clear();
	int32_t p1 = pos;
	for(int k = 0; k < spos.size(); k++)
	{
		int32_t p2 = high32(spos[k]);
		v.push_back(as_pos(pack(p1, p2), spos[k].ale));
		p1 = low32(as_pos(pack(p1, p2), spos[k].ale));
	}
	v.push_back(as_pos(pack(p1, rpos), "$"));
	return 0;
}

string hit::get_qname(bam1_t *b)
{
	char buf[1024];
	char *q = bam_get_qname(b);
	int l = strlen(q);
	memcpy(buf, q, l);
	buf[l] = '\0';
	return string(buf);
}

int hit::set_tags(bam1_t *b)
{
	ts = '.';
	uint8_t *p0 = bam_aux_get(b, "ts");
	if(p0 && (*p0) == 'A') ts = bam_aux2A(p0);
	if(p0 && (*p0) == 'a') ts = bam_aux2A(p0);

	xs = '.';
	uint8_t *p1 = bam_aux_get(b, "XS");
	if(p1 && (*p1) == 'A') xs = bam_aux2A(p1);
	if(p1 && (*p1) == 'a') xs = bam_aux2A(p1);

	if(xs == '.' && ts != '.')
	{
		// convert ts to xs
		if((flag & 0x10) >= 1 && ts == '+') xs = '-';
		if((flag & 0x10) >= 1 && ts == '-') xs = '+';
		if((flag & 0x10) <= 0 && ts == '+') xs = '+';
		if((flag & 0x10) <= 0 && ts == '-') xs = '-';
	}

	hi = -1;
	uint8_t *p2 = bam_aux_get(b, "HI");
	if(p2 && (*p2) == 'C') hi = bam_aux2i(p2);
	if(p2 && (*p2) == 'c') hi = bam_aux2i(p2);

	nh = -1;
	uint8_t *p3 = bam_aux_get(b, "NH");
	if(p3 && (*p3) == 'C') nh = bam_aux2i(p3);
	if(p3 && (*p3) == 'c') nh = bam_aux2i(p3);

	nm = 0;
	uint8_t *p4 = bam_aux_get(b, "nM");
	if(p4 && (*p4) == 'C') nm = bam_aux2i(p4);
	if(p4 && (*p4) == 'c') nm = bam_aux2i(p4);

	uint8_t *p5 = bam_aux_get(b, "NM");
	if(p5 && (*p5) == 'C') nm = bam_aux2i(p5);
	if(p5 && (*p5) == 'c') nm = bam_aux2i(p5);


	// set umi
	umi = "";
	uint8_t *p6 = bam_aux_get(b, "UB");
	if(p6 && (*p6) == 'H') umi = bam_aux2Z(p6);
	if(p6 && (*p6) == 'Z') umi = bam_aux2Z(p6);

	return 0;
}

// int hit::set_concordance()
// {
// 	bool concordant = false;
// 	if((flag & 0x10) <= 0 && (flag & 0x20) >= 1 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) concordant = true;		// F1R2
// 	if((flag & 0x10) >= 1 && (flag & 0x20) <= 0 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) concordant = true;		// R1F2
// 	if((flag & 0x10) <= 0 && (flag & 0x20) >= 1 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) concordant = true;		// F2R1
// 	if((flag & 0x10) >= 1 && (flag & 0x20) <= 0 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) concordant = true;		// R2F1
// 	return 0;
// }

int hit::set_strand()
{
	strand = '.';
	
	if(library_type == FR_FIRST && ((flag & 0x1) >= 1))
	{
		if((flag & 0x10) <= 0 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) strand = '-';
		if((flag & 0x10) >= 1 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) strand = '+';
		if((flag & 0x10) <= 0 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) strand = '+';
		if((flag & 0x10) >= 1 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) strand = '-';
	}

	if(library_type == FR_SECOND && ((flag & 0x1) >= 1))
	{
		if((flag & 0x10) <= 0 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) strand = '+';
		if((flag & 0x10) >= 1 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) strand = '-';
		if((flag & 0x10) <= 0 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) strand = '-';
		if((flag & 0x10) >= 1 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) strand = '+';
	}

	if(library_type == FR_FIRST && ((flag & 0x1) <= 0))
	{
		if((flag & 0x10) <= 0) strand = '-';
		if((flag & 0x10) >= 1) strand = '+';
	}

	if(library_type == FR_SECOND && ((flag & 0x1) <= 0))
	{
		if((flag & 0x10) <= 0) strand = '+';
		if((flag & 0x10) >= 1) strand = '-';
	}

	return 0;
}

bool hit::operator<(const hit &h) const
{
	if(qname < h.qname) return true;
	if(qname > h.qname) return false;
	if(hi != -1 && h.hi != -1 && hi < h.hi) return true;
	if(hi != -1 && h.hi != -1 && hi > h.hi) return false;
	return (pos < h.pos);
}

int hit::print() const
{
	// print basic information
	// printf("Hit %s: chrm %s [%d-%d), mpos = %d, flag = %d, quality = %d, strand = %c, xs = %c, ts = %c, isize = %d, qlen = %d, hi = %d\n", 
			// qname.c_str(), chrm.c_str(), pos, rpos, mpos, flag, qual, strand, xs, ts, isize, qlen, hi);
	printf("Hit %s: hid = %d, chrm %s [%d-%d), mpos = %d, flag = %d, quality = %d, strand = %c, xs = %c, ts = %c, isize = %d, qlen = %d, hi = %d, nh = %d, umi = %s, bridged = %c\n", 
			qname.c_str(), hid, chrm.c_str(), pos, rpos, mpos, flag, qual, strand, xs, ts, isize, qlen, hi, nh, umi.c_str(), bridged ? 'T' : 'F');

	printf(" start position [%d - )\n", pos);
	for(int i = 0; i < spos.size(); i++)
	{
		as_pos p = spos[i];
		as_pos32 p1 = high32(p);
		as_pos32 p2 = low32(p);
		printf(" splice position [%d - %d) ale %s \n", p1.p32, p2.p32, p1.ale.c_str());
	}
	for (auto i: apos)
	{
		as_pos32 p1 = high32(i);
		as_pos32 p2 = low32(i);
		printf(" apos position [%d - %d) ale %s \n", p1.p32, p2.p32, p1.ale.c_str());
	}
	for (auto i: itvm)
	{
		as_pos32 p1 = high32(i);
		as_pos32 p2 = low32(i);
		printf(" itvm position [%d - %d) ale %s \n", p1.p32, p2.p32, p1.ale.c_str());
	}
	for (auto i: itvna)
	{
		as_pos32 p1 = high32(i);
		as_pos32 p2 = low32(i);
		printf(" itvna position [%d - %d) ale %s \n", p1.p32, p2.p32, p1.ale.c_str());
	}
	printf(" end position [ - %d)\n", rpos);
	return 0;
}

vector<int> encode_vlist(const vector<int> &v)
{
	vector<int> vv;
	if(v.size() <= 0) return vv;

	int p = v[0];
	int k = 1;
	for(int i = 1; i < v.size(); i++)
	{
		if(v[i] == v[i - 1] + 1)
		{
			k++;
		}
		else
		{
			assert(k >= 1);
			vv.push_back(p);
			vv.push_back(k);
			p = v[i];
			k = 1;
		}
	}
	vv.push_back(p);
	vv.push_back(k);

	/*
	printf("encode: (");
	printv(v);
	printf(") -> (");
	printv(vv);
	printf(")\n");
	*/
	return vv;
}

vector<int> decode_vlist(const vector<int> &v)
{
	vector<int> vv;
	assert(v.size() % 2 == 0);
	if(v.size() <= 0) return vv;

	for(int i = 0; i < v.size() / 2; i++)
	{
		int p = v[i * 2 + 0];
		int k = v[i * 2 + 1];
		for(int j = p; j < p + k; j++)
		{
			vv.push_back(j);
		}
	}
	return vv;
}
