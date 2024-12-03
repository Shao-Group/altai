/*
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <iostream>
#include <fstream>
#include <map>
#include <vector> 
#include <string>
#include <algorithm> 
#include <cstdlib>
#include <cctype>
#include <boost/algorithm/string.hpp>
#include "config.h"
#include "vcf_data.h"


using namespace std;

const char* gt_str(genotype gt)	
{
	// if(!(gt >= 0 && gt <= 3)) cout << "\n" << "wired gt value " << gt << endl;
	assert(gt >= 0 && gt <= 3);
	size_t i = gt;
	vector<const char*>ss {"UNPHASED", "ALLELE1", "ALLELE2", "NONSPECIFIC"};
	return ss[i];
}

bool gt_conflict(genotype g1, genotype g2) 
{
	set<genotype> gts {UNPHASED, ALLELE1, ALLELE2, NONSPECIFIC};
	assert(gts.find(g1) != gts.end());
	assert(gts.find(g2) != gts.end());
	if ((g1 == ALLELE1 && g2 == ALLELE2) || (g1 == ALLELE2 && g2 == ALLELE1))
		return true;
	else return false;
}

bool gt_explicit_same(genotype g1, genotype g2) 
{
	if (g1 == g2) 
		if (g1 == ALLELE1 || g1 == ALLELE2 || g1 == NONSPECIFIC) 
			return true;
	return false;
}

// true if explicit_same or (UNPHASE, NONSPECIFIC) same. (UNPHASE, ALE1) is false
bool gt_implicit_same(genotype g1, genotype g2) 
{
	if (gt_explicit_same(g1, g2)) return true;
	if (g1 == UNPHASED && g2 == UNPHASED) return true;
	if (g1 == UNPHASED && g2 == NONSPECIFIC) return true;
	if (g2 == UNPHASED && g1 == NONSPECIFIC) return true;
	return false;
}

bool gt_as(genotype g)
{
	if (g == ALLELE1 || g == ALLELE2) 
	{
		return true;
	}	
	else 
	{
		assert( g == NONSPECIFIC || g == UNPHASED);
		return false;
	}
	assert(0); // should never happen
}

vcf_data::vcf_data() {}

vcf_data::vcf_data(std::string file_name) 
{
	read_as_counts(file_name);
}

/*
** @return ALLELE1, ALLELE2, NONSPECIFIC, if not found, UNPHASED 
*/
genotype vcf_data::get_genotype(string chrm, int pos, string ale)
{
	auto vcf1 = vcf_map.find(chrm);
	if(vcf1 == vcf_map.end()) return UNPHASED;
	
	auto vcf2 = vcf1->second.find(pos);
	if(vcf2 == vcf1->second.end()) return UNPHASED;
	
	auto vcf3 = vcf2->second.find(ale);
	if(vcf3 == vcf2->second.end()) return UNPHASED;
	
	genotype gt = vcf3->second;
	return gt;
}

/*
** map < str chrm, map <int pos, map <str nt, genotype> vcf_pos_map; 
** assumption: var.len == 1, sorted; lf[8-9] must have GT first, per vcf specification
** parameters: use_phased_var_only == true
** TODO: check phase set (PS) to ensure phased variants are in the same PS
** TODO: assert vcf header chr overlap w. bam header chr, otherwise throw check chr names
*/
int vcf_data::read_as_counts(const string & name) 
{
	ifstream vcf_count_file(name);
	if (vcf_count_file.is_open()) 
	{
		string l;
		while (!getline(vcf_count_file, l).eof()) 
		{
			if (l.find("#") == 0 || l.empty()) continue;

			// lf := line.split('\t')
			vector<string> lf;  
			boost::trim(l);
			size_t delim_pos = 0;
			while ((delim_pos = l.find("\t")) != string::npos)  
			{
				lf.push_back(l.substr(0, delim_pos));
				l.erase(0, delim_pos + 1);
			}
			lf.push_back(l);
			
			string chrm = lf[0];
			int pos = stoi(lf[1]) -1;	// 0-based

			// get alleles
			vector<string> alleles;
			string a = toupperstring(lf[3]);
			alleles.push_back(toupperstring(lf[3]));	// reference allele
			string lfa = lf[4];					// alt alleles with comma-seperation
			while((delim_pos = lfa.find(",")) != string::npos) 
			{
				alleles.push_back(toupperstring(lfa.substr(0, delim_pos)));
				lfa.erase(0, delim_pos + 1);
			}
			alleles.push_back(toupperstring(lfa));

			// get gt
			// the first filed must always be "GT"
			assert (lf[8].substr(0,2) == "GT");
			map<string, genotype> ng;
			string gt_fields = lf[9];  // = "1|0:555:97,59:31,24:548:0/1:.:PATMAT"

			// if (use_phased_var_only && gt_fields[1] != '|') continue;

			int i1 = gt_fields[0] - '0'; // allele in gt1
			int i2 = gt_fields[2] - '0'; // allele in gt2
			assert (alleles.size() >= 1);
			assert (alleles[i1].length() >= 1);
			if (i1 == i2)
			{
				ng.insert({alleles[i1], NONSPECIFIC});
			}
			else
			{
				if (gt_fields[1] == '|')
				{					
					ng.insert({alleles[i1], ALLELE1});
					ng.insert({alleles[i2], ALLELE2});						
				}
				else
				{
					for (auto&& a: alleles) ng.insert({a, UNPHASED});
				}
			}					

			vcf_pos_map[chrm][pos] = ng;
			vcf_ale_len[chrm][pos] = alleles[0].size();
		}
		vcf_count_file.close();
	}
	else
	{
		cerr << "Unable to open vcf file " << name << endl;
		throw runtime_error("Unable to open vcf file.");
	} 

	return 0;
}


int vcf_data::increse_it(map <int, map <string, genotype> >::iterator &it1, map <int, int >::iterator &it2)
{
	++it1;
	++it2;
	return 0;
}


int vcf_data::print()
{
	cout << "alleles should be in laxico order. 0-based cord" << endl;
	
	// variants and gt
	for (auto it = vcf_pos_map.begin(); it != vcf_pos_map.end(); it ++)
	{
		for (auto it2 = (it->second).begin(); it2 != (it->second).end(); it2 ++)
		{
			cout << it->first << ":" << it2->first << "  ";
			auto ii = it2->second;
			for(auto&& aa: ii) cout << aa.first << "-" << gt_str(aa.second) << "  ";
			cout << endl;
		}

	}

	// ref variant length
	/*
	for (auto it = vcf_ale_len.begin(); it != vcf_ale_len.end(); it ++)
	{
		for (auto it2 = (it->second).begin(); it2 != (it->second).end(); it2 ++)
		{
			cout << it->first << ":" << it2->first << "  " << it2->second << endl;
		}
	}
	*/
	return 0;
}

string vcf_data::graphviz_gt_color_shape(genotype gt, int vertex_type)
{
	string color = "";
	string shape = "";
	if(vertex_type == EMPTY_VERTEX) shape = "\", shape=\"polygon";

	if 		(vertex_type == PSEUDO_AS_VERTEX && gt == ALLELE1) 	color = "purple1";
	else if (vertex_type == PSEUDO_AS_VERTEX && gt == ALLELE2) 	color = "lime";
	else if (gt == ALLELE1) 									color = "red";
	else if (gt == ALLELE2) 									color = "green4";
	else if (gt == NONSPECIFIC) 								color = "gray";
	else if (gt == UNPHASED) 									color = "black";
	else throw runtime_error("graphviz_gt_color: gt not existing");

	return color + shape;
}