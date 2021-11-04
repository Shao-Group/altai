/*
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

PI32 transcript::get_first_intron() const
{
	if(exons.size() <= 1) return PI32(-1, -1);
	as_pos32 p = exons[0].second;
	as_pos32 q = exons[1].first;
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

bool transcript::intron_chain_match(const transcript &t) const
{
	if(exons.size() != t.exons.size()) return false;
	if(exons.size() <= 1) return false;
	int n = exons.size() - 1;
	if(!exons[0].second.samepos(t.exons[0].second)) return false;
	if(!exons[n].first.samepos(t.exons[n].first)) return false;
	for(int k = 1; k < n - 1; k++)
	{
		if(!exons[k].first.samepos(t.exons[k].first)) return false;
		if(!exons[k].second.samepos(t.exons[k].second)) return false;
	}
	return true;
}

string transcript::label() const
{
	char buf[10240];
	PI32 p = get_bounds();
	sprintf(buf, "%s:%d-%d", seqname.c_str(), p.first.p32, p.second.p32);
	return string(buf);
}

int transcript::write(ostream &fout) const
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
	fout<<"RPKM \""<<RPKM<<"\"; ";
	fout<<"cov \""<<coverage<<"\";"<<endl;

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

int transcript::write_gvf(ostream &fout) const
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
	fout<<"RPKM \""<<RPKM<<"\"; ";
	fout<<"cov \""<<coverage<<"\";"<<endl;

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