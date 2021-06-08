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
#include <boost/algorithm/string.hpp>
#include "config.h"
#include "vcf_data.h"


using namespace std;

vcf_data::vcf_data() {}

vcf_data::vcf_data(std::string file_name) 
{
	read_as_counts(file_name);
}

vcf_data::~vcf_data() {}


int vcf_data::read_as_counts(const string & name) 
{
	// map < std::string, map <int, std::vector <std::string> > > vcf_pos_map;
	// assume vcf is cleaned and sorted
	ifstream vcf_count_file(name);
	if (vcf_count_file.is_open()) 
	{
		string l;
		while (!getline(vcf_count_file, l).eof()) 
		{
			if (l.find("#") == 0 || l.empty()) continue;
			boost::trim(l);
			size_t delim_pos = 0;
			vector<string> lf;
			// split l by '\t' and push_back in lf
			while ((delim_pos = l.find("\t")) != string::npos) 
			{
				lf.push_back(l.substr(0, delim_pos));
				l.erase(0, delim_pos + 1);
			}
			lf.push_back(l);
			
			string chrm = lf[0];

			int pos = stoi(lf[1]) -1;	// 0-based
			vector<string> alleles;
			alleles.push_back(lf[3]);	// reference allele
			std::string lfa = lf[4];	// alt alleles with comma-seperation
			while((delim_pos = lfa.find(",")) != string::npos) 
			{
				alleles.push_back(lfa.substr(0, delim_pos));
				lfa.erase(0, delim_pos + 1);
			}
			alleles.push_back(lfa);
			std::sort (alleles.begin(), alleles.end());		// sort alleles in lexico order
			if (vcf_pos_map.find(chrm) == vcf_pos_map.end())
			{
				map <int, std::vector <std::string> > _m_al;
				_m_al.insert(pair<int, std::vector <std::string> > (pos, alleles));
				vcf_pos_map.insert(pair<std::string, map <int, std::vector <std::string> > > (chrm, _m_al));
				map <int, int> _m_al_len;
				_m_al_len.insert(pair<int, int> (pos, lf[3].size()));
				vcf_ale_len.insert(pair<std::string, map <int, int> > (chrm, _m_al_len));
			}
			else 
			{
				vcf_pos_map[chrm].insert(pair<int, std::vector <std::string> > (pos, alleles));
				vcf_ale_len[chrm].insert(pair<int, int> (pos, lf[3].size()) );
			}
			
		}
		vcf_count_file.close();
	}
	else cout << "Unable to open file" << name << endl;  // TODO: catch IO error

	return 0;
}


int vcf_data::increse_it(map <int, vector <std::string> >::iterator &it1, map <int, int >::iterator &it2)
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
			for(int i = 0; i < (it2->second).size(); i++)
			{
				cout << (it2->second)[i] << "  ";
			}
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
