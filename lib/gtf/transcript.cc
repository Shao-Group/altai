/*
Part of aletsch
(c) 2020 by  Mingfu Shao, The Pennsylvania State University.
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <map>
#include "util.h"
#include "transcript.h"
#include "src/config.h"

transcript::transcript()
{
}

transcript::transcript(const item &e)
{
	throw "transcript::transcript should not be constructed with an item object";
	assign(e);
	exons.clear();
	as_exons.clear();
}

transcript::~transcript()
{
}

int transcript::assign(const item &e)
{
	//assert(e.feature == "transcript");
	seqname = e.seqname;
	source = e.source;
	feature = e.feature;
	gene_id = e.gene_id;
	transcript_id = e.transcript_id;
	transcript_type = e.transcript_type;
	gene_type = e.gene_type;
	start = e.start;
	end = e.end;
	strand = e.strand;
	frame = e.frame;
	score = e.score;
	coverage = e.coverage;
	RPKM = e.RPKM;
	FPKM = e.FPKM;
	TPM = e.TPM;
	return 0;
}

bool transcript::operator< (const transcript &t) const
{
	assert(gt_implicit_same(gt, t.gt));
	int b = seqname.compare(t.seqname);
	if(b < 0) return true;
	if(b > 0) return false;
	if(exons.size() == 0) return true;
	if(t.exons.size() == 0) return false;
	if(exons[0].first < t.exons[0].first) return true;
	else return false;
}

int transcript::clear()
{
	exons.clear();
	as_exons.clear();
	seqname = "";
	source = "";
	feature = "";
	gene_id = "";
	transcript_id = "";
	transcript_type = "";
	gene_type = "";
	start = as_pos32(0);
	end = as_pos32(0);
	strand = '.';
	frame = -1;
	coverage = 0;
	RPKM = 0;
	TPM = 0;
	gt = UNPHASED;
	return 0;
}

int transcript::add_exon(as_pos32 s, as_pos32 t)
{
	assert(s.ale == "$");
	assert(t.ale == "$");
	exons.push_back(PI32(s, t));
	return 0;
}

int transcript::add_exon(const item &e)
{
	assert(e.transcript_id == transcript_id);
	add_exon(e.start, e.end);
	return 0;
}

int transcript::add_as_exons(as_pos32 s, as_pos32 t)
{
	assert(s.ale != "$");
	assert(t.ale != "$");
	as_exons.push_back(PI32(s, t));
	return 0;
}

int transcript::sort()
{
	std::sort(exons.begin(), exons.end());
	std::sort(as_exons.begin(), as_exons.end());
	return 0;
}

int transcript::shrink()
{
	if(exons.size() == 0) return 0;
	vector<PI32> v;
	PI32 p = exons[0];
	for(int i = 1; i < exons.size(); i++)
	{
		PI32 q = exons[i];
		if(p.second.samepos(q.first))
		{
			if(DEBUG_MODE_ON)
			{
				assert(q.first.ale == "$"); 
				assert(q.second.ale == "$"); 
				assert(p.first.ale == "$");			
				assert(p.second.ale == "$");
			}
			p.second= q.second;
		}
		else
		{
			// assert(p.second < q.first);
			v.push_back(p);
			p = q;
		}
	}
	v.push_back(p);
	exons = v;
	return 0;
}

int transcript::assign_RPKM(double factor)
{
	RPKM = coverage * factor;
	return 0;
}

int transcript::length() const
{
	int s = 0;
	for(int i = 0; i < exons.size(); i++)
	{
		assert(exons[i].second > exons[i].first);
		s += exons[i].second - exons[i].first;
	}
	return s;
}

int transcript::make_non_specific()
{
	gt = NONSPECIFIC;
	as_exons.clear();
	for (const PI32 & e: exons)
	{
		assert(e.first.ale == "$");
		assert(e.second.ale == "$");
	}
	return 0;
}

int transcript::assign_gt(genotype g)
{	
	assert(!gt_conflict(gt, g));
	gt = g;
	return 0;
}

/*
*	transform gt and as_exon
*	@return	!gt_explicit_same(this->gt, g)
*/
bool transcript::transform_gt(genotype g)
{	
	if(gt_explicit_same(gt, g)) return true;
	
	gt = g;

	if(g == NONSPECIFIC || g == UNPHASED) 
	{
		make_non_specific();
		return false;
	}

	vector<PI32> as_new;
	for(PI32 a: as_exons)
	{
		assert(a.first.ale == a.second.ale);
		int p = a.first.p32;
		int q = a.second.p32;
		string s = "$";
		for(const auto & string_gt_pair: asp.vcf_pos_map[seqname][p])
		{
			if(string_gt_pair.second == g) s = string_gt_pair.first;
		}
		as_pos32 p2(p, s);
		as_pos32 q2(q, s);
		as_new.push_back({p2, q2});
	}
	as_exons = as_new;

	return false;
}

PI32 transcript::get_bounds() const
{
	if(exons.size() == 0) return PI32(-1, -1);
	as_pos32 p = exons[0].first;
	as_pos32 q = exons[exons.size() - 1].second;
	return PI32(p, q);
}

vector<PI32> transcript::get_intron_chain() const
{
	vector<PI32> v;
	if(exons.size() <= 1) return v;

	as_pos32 p = exons[0].second;
	for(int k = 1; k < exons.size(); k++)
	{
		as_pos32 q = exons[k].first;
		v.push_back(PI32(p, q));
		p = exons[k].second;
	}
	return v;
}

size_t transcript::get_intron_chain_hashing() const
{
	if(exons.size() == 0) return 0;

	if(exons.size() == 1)
	{
		assert(get_intron_chain().size() == 0);
		vector<as_pos32> vv;
		vv.push_back(as_pos32{-2, "SingleExon"});
		vv.push_back(exons[0].first);
		vv.push_back(exons[0].second);
		return vector_hash(vv) + 1;
	}

	vector<as_pos32> vv;
	vector<PI32> v = get_intron_chain();
	for(int i = 0; i < v.size(); i++) 
	{
		vv.push_back(v[i].first);
		vv.push_back(v[i].second);
	}
	return vector_hash(vv) + 1;
}

int transcript::intron_chain_compare(const transcript &t) const
{
	assert(!gt_conflict(gt, t.gt));
	if (DEBUG_MODE_ON)  for(const auto & e: exons)  assert(e.first.ale == "$" && e.second.ale == "$" );
	
	if(exons.size() < t.exons.size()) return +1;
	if(exons.size() > t.exons.size()) return -1;
	if(exons.size() <= 1) return 0;

	int n = exons.size() - 1;
	if(exons[0].second < t.exons[0].second) return +1;
	if(exons[0].second > t.exons[0].second) return -1;
	for(int k = 1; k < n - 1; k++)
	{
		if(exons[k].first < t.exons[k].first) return +1;
		if(exons[k].first > t.exons[k].first) return -1;
		if(exons[k].second < t.exons[k].second) return +1;
		if(exons[k].second > t.exons[k].second) return -1;
	}
	if(exons[n].first < t.exons[n].first) return +1;
	if(exons[n].first > t.exons[n].first) return -1;
	return 0;
}

int transcript::compare1(const transcript &t, double single_exon_overlap) const
{
	assert(!gt_conflict(gt, t.gt));
	if(exons.size() < t.exons.size()) return +1;
	if(exons.size() > t.exons.size()) return -1;

	if(seqname < t.seqname) return +1;
	if(seqname > t.seqname) return -1;
	if(strand < t.strand) return +1;
	if(strand > t.strand) return -1;

	if (DEBUG_MODE_ON)  for(const auto & e: exons)  assert(e.first.ale == "$" && e.second.ale == "$" );
		
	if(exons.size() == 1)
	{
		// int32_t p1 = exons[0].first < t.exons[0].first ? exons[0].first.p32 : t.exons[0].first.p32;
		// int32_t p2 = exons[0].first < t.exons[0].first ? t.exons[0].first.p32 : exons[0].first.p32;
		// int32_t q1 = exons[0].second > t.exons[0].second ? exons[0].second.p32 : t.exons[0].second.p32;
		// int32_t q2 = exons[0].second > t.exons[0].second ? t.exons[0].second.p32 : exons[0].second.p32;
		as_pos32 p1 = exons[0].first < t.exons[0].first ? exons[0].first.p32 : t.exons[0].first.p32;
		as_pos32 p2 = exons[0].first < t.exons[0].first ? t.exons[0].first.p32 : exons[0].first.p32;
		as_pos32 q1 = exons[0].second > t.exons[0].second ? exons[0].second.p32 : t.exons[0].second.p32;
		as_pos32 q2 = exons[0].second > t.exons[0].second ? t.exons[0].second.p32 : exons[0].second.p32;

		int32_t overlap = q2.p32 - p2.p32;
		if(overlap >= single_exon_overlap * length()) return 0;
		if(overlap >= single_exon_overlap * t.length()) return 0;

		//double overlap = (q2 - p2) * 1.0 / (q1 - p1);
		//if(overlap >= 0.8) return 0;

		if(exons[0].first < t.exons[0].first) return +1;
		if(exons[0].first > t.exons[0].first) return -1;
		if(exons[0].second < t.exons[0].second) return +1;
		if(exons[0].second > t.exons[0].second) return -1;
	}

	return intron_chain_compare(t);
}

int transcript::extend_bounds(const transcript &t)
{
	assert(gt_implicit_same(gt, t.gt));
	if(exons.size() == 0) return 0;
	if(t.exons.front().first < exons.front().first) exons.front().first = t.exons.front().first;
	if(t.exons.back().second > exons.back().second) exons.back().second = t.exons.back().second;
	return 0;
}

int transcript::write(ostream &fout, double cov2, int count) const
{
	fout.precision(4);
	fout<<fixed;

	if(exons.size() == 0) return 0;
	
	PI32 p = get_bounds();

	fout	<<	seqname.c_str()		<<"\t";										// chromosome name
	fout	<<	source.c_str()		<<"\t";										// source
	fout	<<	"transcript"		<<"\t";										// feature
	fout	<<	p.first.p32 + 1		<<"\t";										// left position; 0-based to 1-based by adding 1
	fout	<<	p.second.p32 		<<"\t";										// right position
	fout	<<	1000				<<"\t";										// score, now as expression
	fout	<<	strand				<<"\t";										// strand
	fout	<<	"."					<<"\t";															// frame
	fout	<<	"gene_id \""		<<	gene_id.c_str()				<<"\"; ";
	fout	<<	"transcript_id \""	<<	transcript_id.c_str()		<<"\"; ";
	if(gene_type != "") 
		fout<<	"gene_type \""		<<	gene_type.c_str()			<<"\"; ";
	if(transcript_type != "") 
		fout<<"transcript_type \""	<<	transcript_type.c_str()		<<"\"; ";
	
	//fout<<"RPKM \""<<RPKM<<"\"; ";
	fout	<<	"cov \""			<<	coverage					<<"\"; ";

	if(cov2 >= -0.5) 
		fout<<	"cov2 \""			<<	cov2						<<"\"; ";
	
	if(count >= -0.5) 
		fout<<	"count \""			<<	count						<<"\"; ";
	
	fout 	<< 	endl;

	for(int k = 0; k < exons.size(); k++)
	{
		fout<<seqname.c_str()<<"\t";		// chromosome name
		fout<<source.c_str()<<"\t";			// source
		fout<<"exon\t";						// feature
		fout<<exons[k].first.p32 + 1<<"\t";		// left position; 0-based to 1-based by adding 1
		fout<<exons[k].second.p32<<"\t";		// right position
		fout<<1000<<"\t";					// score, now as expression
		fout<<strand<<"\t";					// strand
		fout<<".\t";						// frame
		fout<<"gene_id \""<<gene_id.c_str()<<"\"; ";
		fout<<"transcript_id \""<<transcript_id.c_str()<<"\"; ";
		fout<<"exon \""<<k + 1<<"\"; "<<endl;
	}
	return 0;
}

int transcript::write_gvf(ostream &fout, double cov2, int count) const
{
	fout.precision(4);
	fout<<fixed;

	if(exons.size() == 0) return 0;
	vector<PI32> exons_and_as_exons(exons);
	exons_and_as_exons.insert(exons_and_as_exons.end(), as_exons.begin(), as_exons.end());
	std::sort(exons_and_as_exons.begin(), exons_and_as_exons.end());

	PI32 p = get_bounds();

	fout<<seqname.c_str()<<"\t";				// chromosome name
	fout<<source.c_str()<<"\t";					// source
	fout<<"transcript\t";						// feature
	fout<<p.first.p32 + 1<<"\t";					// left position
	fout<<p.second.p32 <<"\t";						// right position
	fout<<1000<<"\t";							// score, now as expression
	fout<<strand<<"\t";							// strand
	fout<<".\t";								// frame
	fout<<"gene_id \""<<gene_id.c_str()<<"\"; ";
	fout<<"transcript_id \""<<transcript_id.c_str()<<"\"; ";
	fout<<"allele \""<< gt_str(gt) <<"\"; ";
	if(gene_type != "") fout<<"gene_type \""<<gene_type.c_str()<<"\"; ";
	if(transcript_type != "") fout<<"transcript_type \""<<transcript_type.c_str()<<"\"; ";
	//fout<<"RPKM \""<<RPKM<<"\"; ";
	fout<<"cov \""<<coverage<<"\"; ";
	if(cov2 >= -0.5) fout<<"cov2 \""<<cov2<<"\"; ";
	if(count >= -0.5) fout<<"count \""<<count<<"\"; ";
	fout << endl;

	for(int k = 0, exon_num = 0; k < exons_and_as_exons.size(); k++)
	{
		fout<<seqname.c_str()<<"\t";		// chromosome name
		fout<<source.c_str()<<"\t";			// source
		
		string a = exons_and_as_exons[k].first.ale;
		if (a == "$") 
		{
			fout<<"exon\t";						// feature
			exon_num ++;
		}
		else 
		{
			fout<<"variant\t";	
		}
		fout<<exons_and_as_exons[k].first.p32 + 1<<"\t";		// left position
		fout<<exons_and_as_exons[k].second.p32<<"\t";		// right position
		fout<<1000<<"\t";					// score, now as expression
		fout<<strand<<"\t";					// strand
		fout<<".\t";						// frame
		fout<<"gene_id \""<<gene_id.c_str()<<"\"; ";
		fout<<"transcript_id \""<<transcript_id.c_str()<<"\"; ";
		fout<<"exon \""<<exon_num <<"\"; ";
		if (a != "$") 
		{
			fout<<"seq \""<<a.c_str()<<"\"; ";
			fout<<"allele \""<< gt_str(gt) <<"\"; ";
		}
		fout<<endl;
	}
	return 0;
}

int transcript::write_fasta(ostream &fout, int line_len, faidx_t *fai) const
{
	assert(0 && "Not implemented");		// TODO: pull sequences from genome fasta, replace allele


	if(exons.size() == 0) return 0;
		
	string sequence;
	for(int k = 0; k < exons.size(); ++k)
	{
		assert(exons[k].first.ale == "$");
		assert(exons[k].second.ale == "$");
		int seqlen;
		char* s = faidx_fetch_seq(fai, seqname.c_str(), exons[k].first.p32, exons[k].second.p32-1, &seqlen);	//both [seq_begin, seq_end] included; both exons and fai are 0-based
		sequence += s;
		free(s);	
	}

	// int k = 0;
	// int i = 0;
	// while(i < as_exons.size() && k < exons.size())
	// {
	// 	int p1 = as_exons[i].first.p32;
	// 	int p2 = as_exons[i].second.p32;
	// 	int q1 = exons[k].first.p32;
	// 	int q2 = exons[k].second.p32;
		
	// 	if() 
	// 	{
			
	// 	}
	// 	else if ()
	// 	{

	// 	}
		

	// 	assert(k >= 1);
	// 	assert(exons[k].first.ale != "$")
	// 	assert(exons[k].first.ale == exons[k].second.ale);
	// 	assert(exons[k].first.rightsameto(exons[k-1].first));
	// 	int l = sequence.length();
	// 	string ale = exons[k].first.ale;
	// 	int exon_len = exons[k-1].second.p32 - exons[k-1].first.p32;
	// 	int p = l - exon_len + (exons[k].first.p32 - exons[k-1].first.p32);
	// 	int as_len = exons[k].second.p32 - exons[k].first.p32;
	// 	sequence.replace(p, as_len, exons[k].first.ale);

	// 	i ++;
	// }
	// if(i < as_exons.size())
	// {

	// }


	// string rc;
	// reverse_complement_DNA(rc, sequence);
	// if (strand == "-")
	// {
	// }
	// else if (strand == '.'){	//TODO: use both
	// }
	// else 
	// {
	// 	assert (strand == '+'); 	//trand must be one of ./+/- 
	// }


	string name = transcript_id;	

	fout << ">" << name << endl;
	int i = 0;
	for(; i < sequence.length() / line_len; ++i)  fout << sequence.substr(line_len * i, line_len) << endl;
	fout << sequence.substr(line_len * i, sequence.length() - line_len * i) << endl;
	return 0;	
}


int transcript::reverse_complement_DNA(string &rc, const string s)
{
	char c;
	for(int i = s.length() - 1; i >= 0; --i)
	{
		switch (s[i])
		{
		case 'A':
			c = 'T';
			break;
		case 'T':
			c = 'A';
			break;
		case 'C':
			c = 'G';
			break;
		case 'G':
			c = 'C';		
			break;
		default:
			c = 'N'; //TODO: specify other degenerate nt
			break;
		}
		rc += c;
	}
	return 0;
}