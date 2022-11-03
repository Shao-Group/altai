/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __VERTEX_INFO__
#define __VERTEX_INFO__

#include <stdint.h>
#include "as_pos32.hpp"
#include "vcf_data.h"

#define NS_NONVAR 0
#define	NS_MONOVAR 1
#define AS_DIPLOIDVAR 2
#define AS_MONOPLOIDVAR 3
#define	UHPHASED_MONOVAR 4
#define AJ_NSMONOVAR 5
#define AJ_NONVAR 6
#define START_OR_SINK 7

class vertex_info
{
public:
	vertex_info();
	vertex_info(int l);
	vertex_info(const vertex_info &vi);

public:
	int32_t pos;		// position - not used
	as_pos32 lpos;		// left position
	as_pos32 rpos;		// right position
	double stddev;		// standard deviation of read coverage
	int length;			// length of this partial exon
	int sdist;			// shortest distance to s
	int tdist;			// shortest distance to t
	int type;			// for various usage
	char lstrand;		// left side strand
	char rstrand;		// right side strand	
	bool regional;		// if a vertex is regional
	genotype gt;

	/*  
		3 bits type :	
		0-NS non-var	
		1-NS monoploid var
		2-AS diploid var 
		3-AS monoploid  var 
		4-UHPHASED monoploid var (assert no when assuming all phased)
		5-AJ NS monoploid var adjacent to AS var
		6-AJ Non-var adjacent to AS var		
		7-start/sink
	*/
	unsigned int as_type: 3;

public:
	bool is_as_vertex();
	bool is_adjacent_to_as_vertex();
	bool is_ordinary_vertex();
	
	static bool is_as_vertex(vertex_info vi);
	static bool is_adjacent_to_as_vertex(vertex_info vi);
	static bool is_ordinary_vertex(vertex_info vi);
};

#endif
