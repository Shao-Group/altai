/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __INTERVAL_MAP_H__
#define __INTERVAL_MAP_H__

// boost::interval map
#include "boost/icl/interval_map.hpp"
#include "boost/icl/split_interval_map.hpp"
#include "as_pos32.hpp"

#include <vector>

using namespace boost;
using namespace std;

typedef icl::right_open_interval<as_pos32> ROI;

// join interval map
typedef icl::interval_map<as_pos32, int, icl::partial_absorber, less, icl::inplace_plus, icl::inter_section, ROI> join_interval_map;
typedef join_interval_map::const_iterator JIMI;
typedef pair<JIMI, JIMI> PJIMI;
typedef icl::interval_map<int32_t, int, icl::partial_absorber, less, icl::inplace_plus, icl::inter_section, ROI> join_interval_map_int;
typedef join_interval_map_int::const_iterator JIMI_int;

// split interval map
typedef icl::split_interval_map<as_pos32, int, icl::partial_absorber, less, icl::inplace_plus, icl::inter_section, ROI> split_interval_map; 
typedef split_interval_map::const_iterator SIMI;
typedef pair<SIMI, SIMI> PSIMI;


// return the overlap at position p
int compute_overlap(const split_interval_map &imap, as_pos32 p);

// find the leftmost iterator whose upper posistion <= x
SIMI locate_right_iterator(const split_interval_map &imap, as_pos32 x);

// find the rightmost interval whose lower position >= x
SIMI locate_left_iterator(const split_interval_map &imap, as_pos32 x);

// locate boundary iterators
PSIMI locate_boundary_iterators(const split_interval_map &imap, as_pos32 x, as_pos32 y);

// return the sum of the lengths of intervals from p to q (include q)
int compute_coverage(const split_interval_map &imap, SIMI &p, SIMI &q);

// return the maximum overlap of the intervals from p to q (include q)
int compute_max_overlap(const split_interval_map &imap, SIMI &p, SIMI &q);

// return the sum of the overlap of the intervals from p to q (include q)
int compute_sum_overlap(const split_interval_map &imap, SIMI &p, SIMI &q);

// evaluate a region
int evaluate_rectangle(const split_interval_map &imap, as_pos32 ll, as_pos32 rr, double &ave, double &dev);
int evaluate_triangle(const split_interval_map &imap, as_pos32 ll, as_pos32 rr, double &ave, double &dev);

// testing
int test_split_interval_map();

#endif
