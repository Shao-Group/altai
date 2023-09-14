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
#include "config.h"

transcript::transcript()
{
}

transcript::transcript(const item &e)
{
	assign(e);
	exons.clear();
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
	return 0;
}

int transcript::add_exon(as_pos32 s, as_pos32 t)
{
	exons.push_back(PI32(s, t));
	return 0;
}

int transcript::add_exon(const item &e)
{
	assert(e.transcript_id == transcript_id);
	add_exon(e.start, e.end);
	return 0;
}

int transcript::sort()
{
	std::sort(exons.begin(), exons.end());
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
			if (q.second.ale != "$" || p.second.ale != "$" )
			{
				extern faidx_t* fai;
				int seqlen1, seqlen2;
				char * s  = faidx_fetch_seq(fai, seqname.c_str(), p.first.p32, p.second.p32-1, &seqlen1);	//both [seq_begin, seq_end] included
				char * s2 = faidx_fetch_seq(fai, seqname.c_str(), q.first.p32, q.second.p32-1, &seqlen2);
				if (seqlen1 < 0 || seqlen2 < 0) 
					;//TODO: when seq N/A
				string seq;
				seq += s;
				seq += s2;
				free(s);
				free(s2);
				cout << "shrink triggered" << seq << endl;
				p.second.p32 = q.second.p32;
				p.second.ale = seq;
				p.first.ale = seq;
			}
			else p.second.p32 = q.second.p32;
		}
		else
		{
			//assert(p.second < q.first);
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

PI32 transcript::get_bounds() const
{
	if(exons.size() == 0) return PI32(-1, -1);
	as_pos32 p = exons[0].first;
	as_pos32 q = exons[exons.size() - 1].second;
	return PI32(p, q);
}

// PI32 transcript::get_first_intron() const
// {
// 	if(exons.size() <= 1) return PI32(-1, -1);
// 	as_pos32 p = exons[0].second;
// 	as_pos32 q = exons[1].first;
// 	return PI32(p, q);
// }

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

// bool transcript::intron_chain_match(const transcript &t) const
// {
// 	if(exons.size() != t.exons.size()) return false;
// 	if(exons.size() <= 1) return false;
// 	int n = exons.size() - 1;
// 	if(exons[0].second != (t.exons[0].second)) return false;
// 	if(exons[n].first != (t.exons[n].first)) return false;
// 	for(int k = 1; k < n - 1; k++)
// 	{
// 		if(exons[k].first != (t.exons[k].first)) return false;
// 		if(exons[k].second != (t.exons[k].second)) return false;
// 	}
// 	return true;
// }

int transcript::intron_chain_compare(const transcript &t) const
{
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

// bool transcript::equal1(const transcript &t, double single_exon_overlap) const
// {
// 	if(exons.size() != t.exons.size()) return false;

// 	if(seqname != t.seqname) return false;
// 	if(strand == '+' && t.strand == '-') return false;
// 	if(strand == '-' && t.strand == '+') return false;

// 	if(exons.size() == 1)
// 	{
// 		// single exon transcripts have no allele info
// 		int32_t p1 = exons[0].first < t.exons[0].first ? exons[0].first.p32 : t.exons[0].first.p32;
// 		int32_t p2 = exons[0].first < t.exons[0].first ? t.exons[0].first.p32 : exons[0].first.p32;
// 		int32_t q1 = exons[0].second > t.exons[0].second ? exons[0].second.p32 : t.exons[0].second.p32;
// 		int32_t q2 = exons[0].second > t.exons[0].second ? t.exons[0].second.p32 : exons[0].second.p32;

// 		int32_t overlap = q2 - p2;
// 		if(overlap >= single_exon_overlap * length()) return true;
// 		if(overlap >= single_exon_overlap * t.length()) return true;
// 		return false;

// 		/*
// 		double overlap = (q2 - p2) * 1.0 / (q1 - p1);
// 		if(overlap < 0.8) return false;
// 		else return true;
// 		*/
// 	}

// 	return intron_chain_match(t);
// }

// TODO: need to properly handle AS-single-exon transcripts, which are represented as 3 "exons"
int transcript::compare1(const transcript &t, double single_exon_overlap) const
{
	if(exons.size() < t.exons.size()) return +1;
	if(exons.size() > t.exons.size()) return -1;

	if(seqname < t.seqname) return +1;
	if(seqname > t.seqname) return -1;
	if(strand < t.strand) return +1;
	if(strand > t.strand) return -1;

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
	if(exons.size() == 0) return 0;
	if(t.exons.front().first < exons.front().first) exons.front().first = t.exons.front().first;
	if(t.exons.back().second > exons.back().second) exons.back().second = t.exons.back().second;
	return 0;
}

// string transcript::label() const
// {
// 	char buf[10240];
// 	PI32 p = get_bounds();
// 	sprintf(buf, "%s:%d-%d", seqname.c_str(), p.first.p32, p.second.p32);
// 	return string(buf);
// }

int transcript::write(ostream &fout, double cov2, int count) const
{
	fout.precision(4);
	fout<<fixed;

	if(exons.size() == 0) return 0;
	
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
	if(gene_type != "") fout<<"gene_type \""<<gene_type.c_str()<<"\"; ";
	if(transcript_type != "") fout<<"transcript_type \""<<transcript_type.c_str()<<"\"; ";
	//fout<<"RPKM \""<<RPKM<<"\"; ";
	fout<<"cov \""<<coverage<<"\"; ";
	if(cov2 >= -0.5) fout<<"cov2 \""<<cov2<<"\"; ";
	if(count >= -0.5) fout<<"count \""<<count<<"\"; ";
	fout << endl;

	for(int k = 0; k < exons.size(); k++)
	{
		fout<<seqname.c_str()<<"\t";		// chromosome name
		fout<<source.c_str()<<"\t";			// source
		fout<<"exon\t";						// feature
		fout<<exons[k].first.p32 + 1<<"\t";		// left position
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
	if(gene_type != "") fout<<"gene_type \""<<gene_type.c_str()<<"\"; ";
	if(transcript_type != "") fout<<"transcript_type \""<<transcript_type.c_str()<<"\"; ";
	//fout<<"RPKM \""<<RPKM<<"\"; ";
	fout<<"cov \""<<coverage<<"\"; ";
	if(cov2 >= -0.5) fout<<"cov2 \""<<cov2<<"\"; ";
	if(count >= -0.5) fout<<"count \""<<count<<"\"; ";
	fout << endl;

	for(int k = 0; k < exons.size(); k++)
	{
		fout<<seqname.c_str()<<"\t";		// chromosome name
		fout<<source.c_str()<<"\t";			// source
		
		string a = exons[k].first.ale;
		if (a == "$") fout<<"exon\t";						// feature
		else fout<<"variant\t";	
		fout<<exons[k].first.p32 + 1<<"\t";		// left position
		fout<<exons[k].second.p32<<"\t";		// right position
		fout<<1000<<"\t";					// score, now as expression
		fout<<strand<<"\t";					// strand
		fout<<".\t";						// frame
		fout<<"gene_id \""<<gene_id.c_str()<<"\"; ";
		fout<<"transcript_id \""<<transcript_id.c_str()<<"\"; ";
		fout<<"exon \""<<k + 1<<"\"; ";
		if (a != "$") 
		{
			fout<<"seq \""<<a.c_str()<<"\"; ";
			fout<<"allele \""<<allele.c_str()<<"\"; ";
		}
		fout<<endl;
	}
	return 0;
}

int transcript::write_fasta(ostream &fout, int line_len, faidx_t *fai) const
{
	if(exons.size() == 0) return 0;

	string name = transcript_id;	
	
	int seqlen;
	if (strand != '.')
	{
		string sequence;
		if (strand == '+')
		{
			for(int k = 0; k < exons.size(); ++k)
			{
				if (exons[k].first.ale == "$")
				{
					char * s = faidx_fetch_seq(fai, seqname.c_str(), exons[k].first.p32, exons[k].second.p32-1, &seqlen);	//both [seq_begin, seq_end] included
					sequence += s;
					free(s);	
				}
				else 
				{
					if (exons[k].first.ale == "*") sequence += "N";
					else sequence += exons[k].first.ale;
				}
			}
		}
		else if (strand == '-')
		{
			for(int k = exons.size() - 1; k >= 0 ; --k)
			{
				string rc;
				if (exons[k].first.ale == "$")
				{
					char * s = faidx_fetch_seq(fai, seqname.c_str(), exons[k].first.p32, exons[k].second.p32-1, &seqlen);	//both [seq_begin, seq_end] included
					reverse_complement_DNA(rc, s);
					free(s);
				}
				else 
				{
					string s = exons[k].first.ale;
					if (s == "*") s = "N";
					reverse_complement_DNA(rc, s);
				}
				sequence += rc;
			}
		}
		else {cout << "\n"<< strand; assert(0 && "Error strand must be one of ./+/- "); }

		fout << ">" << name << endl;
		int i = 0;
		for(; i < sequence.length() / line_len; ++i)  fout << sequence.substr(line_len * i, line_len) << endl;
		fout << sequence.substr(line_len * i, sequence.length() - line_len * i) << endl;
	}
	else  // When strand un-determined, unlikely to happen
	{
		string sequence;
		string sequence_rc;
		for(int k = 0; k < exons.size(); ++k)
		{
			string seq;
			if (exons[k].first.ale == "$")
			{
				char * s = faidx_fetch_seq(fai, seqname.c_str(), exons[k].first.p32, exons[k].second.p32-1, &seqlen);	//both [seq_begin, seq_end] included
				seq = s;
				free(s);
			}
			else 
			{
				if(exons[k].first.ale == "*") seq = "N";
				else seq = exons[k].first.ale;
			}
			sequence += seq;
			string rc;
			reverse_complement_DNA(rc, seq);
			sequence_rc.insert(0, rc);
		}
		fout << ">" << name << "_forward" << endl;
		int i = 0;
		for(; i < sequence.length() / line_len; ++i)  fout << sequence.substr(line_len * i, line_len) << endl;
		fout << sequence.substr(line_len * i, sequence.length() - line_len * i) << endl;

		fout << ">" << name  << "_reverse" << endl;
		i = 0;
		for(; i < sequence_rc.length() / line_len; ++i)  fout << sequence_rc.substr(line_len * i, line_len) << endl;
		fout << sequence_rc.substr(line_len * i, sequence_rc.length() - line_len * i) << endl;
	}

	return 0;	
}


int transcript::write_fasta_AS_only(ostream &fout, int line_len, faidx_t *fai) const
{
	if(exons.size() == 0) return 0;
	// is transcript AS //TODO: time complexity optimize
	bool _b = false;
	for(int i = 0; i < exons.size(); ++i) if(exons[i].first.ale != "$") {_b = true; break;}
	if (! _b) return 0;

	string name = transcript_id;	
	
	int seqlen;
	if (strand != '.')
	{
		string sequence;
		if (strand == '+')
		{
			for(int k = 0; k < exons.size(); ++k)
			{
				if (exons[k].first.ale == "$")
				{
					char * s = faidx_fetch_seq(fai, seqname.c_str(), exons[k].first.p32, exons[k].second.p32-1, &seqlen);	//both [seq_begin, seq_end] included
					sequence += s;
					free(s);	
				}
				else 
				{
					if (exons[k].first.ale == "*") sequence += "N";
					else sequence += exons[k].first.ale;
				}
			}
		}
		else if (strand == '-')
		{
			for(int k = exons.size() - 1; k >= 0 ; --k)
			{
				string rc;
				if (exons[k].first.ale == "$")
				{
					char * s = faidx_fetch_seq(fai, seqname.c_str(), exons[k].first.p32, exons[k].second.p32-1, &seqlen);	//both [seq_begin, seq_end] included
					reverse_complement_DNA(rc, s);
					free(s);
				}
				else 
				{
					string s = exons[k].first.ale;
					if (s == "*") s = "N";
					reverse_complement_DNA(rc, s);
				}
				sequence += rc;
			}
		}
		else assert(0 && "Error strand must be one of ./+/- "); 

		fout << ">" << name << endl;
		int i = 0;
		for(; i < sequence.length() / line_len; ++i)  fout << sequence.substr(line_len * i, line_len) << endl;
		fout << sequence.substr(line_len * i, sequence.length() - line_len * i) << endl;
	}
	else  // When strand un-determined, unlikely to happen
	{
		string sequence;
		string sequence_rc;
		for(int k = 0; k < exons.size(); ++k)
		{
			string seq;
			if (exons[k].first.ale == "$")
			{
				char * s = faidx_fetch_seq(fai, seqname.c_str(), exons[k].first.p32, exons[k].second.p32-1, &seqlen);	//both [seq_begin, seq_end] included
				seq = s;
				free(s);
			}
			else 
			{
				if(exons[k].first.ale == "*") seq = "N";
				else seq = exons[k].first.ale;
			}
			sequence += seq;
			string rc;
			reverse_complement_DNA(rc, seq);
			sequence_rc.insert(0, rc);
		}
		fout << ">" << name << "_forward" << endl;
		int i = 0;
		for(; i < sequence.length() / line_len; ++i)  fout << sequence.substr(line_len * i, line_len) << endl;
		fout << sequence.substr(line_len * i, sequence.length() - line_len * i) << endl;

		fout << ">" << name  << "_reverse" << endl;
		i = 0;
		for(; i < sequence_rc.length() / line_len; ++i)  fout << sequence_rc.substr(line_len * i, line_len) << endl;
		fout << sequence_rc.substr(line_len * i, sequence_rc.length() - line_len * i) << endl;
	}

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
