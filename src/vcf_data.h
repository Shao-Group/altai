/*
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __VCF_DATA_H__
#define __VCF_DATA_H__

#include <stdint.h>
#include <map>
#include <vector> 
#include <string>
#include "htslib/vcf.h"



class vcf_data
{
public:  						
	vcf_data();								
	vcf_data(std::string);
	virtual ~vcf_data();

public:
	map < string, map <int, vector <string> > > vcf_pos_map; // map <string chrm, map<int pos, vector<string var> > >   map of variant posisions and vector_of_variant_sequences
	map < string, map <int, int > > vcf_ale_len;			 // map <string chrm, map<int pos, int length > >		    map of variant positions and lengths_on_reference
	map < string, map <int, vector <int> > > vcf_pos_phase;	 // map <string chrm, map<int pos, vector<int phase> > >    map of variant positions and 0|1|-1 (on which allele or not present)

public:  																		
	static int increse_it(map <int, vector <std::string> >::iterator &it1, map <int, int >::iterator &it2);
	int read_as_counts(const std::string &);																			// read .asf file, make vcf_map and vcf_pos_map
	int find_var_between(string chrm, int lpos, int rpos, vector < map< string, map<int, int> >::iterator > v);			// TODO: find variant positions between lpos and rpos log (n), returned value is in v
	int print(int);
};


#endif