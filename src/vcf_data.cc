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
#include <boost/algorithm/string.hpp>
#include "config.h"
#include "vcf_data.h"


using namespace std;

const char* gt_str(genotype gt)	
{
	assert(gt >= 0 && gt <= 3);
	size_t i = gt;
	vector<char*>ss {"UNPHASED", "ALLELE1", "ALLELE2", "NONSPECIFIC"};
	return ss[gt];
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
	if ((g1 == ALLELE1 && g1 == g2) || (g1 == ALLELE2 && g1 == g2))
		return true;
	else return false;
}

vcf_data::vcf_data() {}

vcf_data::vcf_data(std::string file_name) 
{
	read_as_counts(file_name);
}

int vcf_data::read_as_counts(const string & name) 
{
	// map < str chrm, map <int pos, map <str nt, genotype> vcf_pos_map; 
	// assumption: var.len == 1, sorted; lf[8-9] must have GT first, per vcf specification
	// parameters: use_phased_var_only == true
	// TODO: check phase set (PS) to ensure phased variants are in the same PS
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
			alleles.push_back(lf[3]);	// reference allele
			string lfa = lf[4];	// alt alleles with comma-seperation
			while((delim_pos = lfa.find(",")) != string::npos) 
			{
				alleles.push_back(lfa.substr(0, delim_pos));
				lfa.erase(0, delim_pos + 1);
			}
			alleles.push_back(lfa);

			// get gt
			map<string, genotype> ng;
			string gt_fields = lf[9];  // = "1|0:555:97,59:31,24:548:0/1:.:PATMAT"

			if (use_phased_var_only && gt_fields[1] != '|') continue;

			int i1 = gt_fields[0] - '0'; // allele in gt1
			int i2 = gt_fields[2] - '0'; // allele in gt2
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
			// if (vcf_pos_map.find(chrm) == vcf_pos_map.end())
			// {
			// 	map <int, std::vector <std::string> > _m_al;
			// 	_m_al.insert(pair<int, std::vector <std::string> > (pos, alleles));
			// 	vcf_pos_map.insert(pair<std::string, map <int, std::vector <std::string> > > (chrm, _m_al));
			// 	map <int, int> _m_al_len;
			// 	_m_al_len.insert(pair<int, int> (pos, lf[3].size()));
			// 	vcf_ale_len.insert(pair<std::string, map <int, int> > (chrm, _m_al_len));
			// }
			// else 
			// {
			// 	vcf_pos_map[chrm].insert(pair<int, std::vector <std::string> > (pos, alleles));
			// 	vcf_ale_len[chrm].insert(pair<int, int> (pos, lf[3].size()) );
			// }
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


int vcf_data::print(int k)
{
	cout << "alleles should be in laxico order. 0-based cord" << endl;
	for (auto it = vcf_pos_map.begin(); it != vcf_pos_map.end(); it ++)
	{
		cout << "chrm " << it->first << endl;
		for (auto it2 = (it->second).begin(); it2 != (it->second).end(); it2 ++)
		{
			cout << it2->first << "  ";
			auto ii = it2->second;
			for(auto&& aa: ii) cout << aa.first << "-" << gt_str(aa.second) << "  ";
			cout << endl;
		}

	}
	for (auto it = vcf_ale_len.begin(); it != vcf_ale_len.end(); it ++)
	{
		cout << "chrm " << it->first << endl;
		for (auto it2 = (it->second).begin(); it2 != (it->second).end(); it2 ++)
		{
			cout << it2->first << "  " << it2->second << endl;
		}

	}
	return 0;
}
