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
	std::map < std::string, std::map <int, std::vector <std::string> > > vcf_pos_map;									
	std::map < std::string, std::map <int, int > > vcf_ale_len;

public:  																		
	static int increse_it(map <int, vector <std::string> >::iterator &it1, map <int, int >::iterator &it2);
	int read_as_counts(const std::string &);									// read .asf file, make vcf_map and vcf_pos_map
	int print(int);
};


#endif