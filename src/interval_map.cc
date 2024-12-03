/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include "interval_map.h"
#include "as_pos32.hpp"


int32_t compute_overlap(const split_interval_map &imap, as_pos32 p)
{
	SIMI it = imap.find(ROI(p, p+1));
	if(it == imap.end()) return 0;
	return it->second;
}

// locate iterator to the right of x - 1. If x is a left boundary, returned ROI's lower = x  
// if (x-1, x) is in the middle of an interval, should return this interval (one left upper_bound)
// lower(it->first) <= x < upper(it->first)
SIMI locate_right_iterator(const split_interval_map &imap, as_pos32 x)
{
	SIMI it = imap.upper_bound(ROI(as_pos32(x - 1), as_pos32(x))); 
	return it;
}

// locate iterator to the left of x - 1, aka contains x in right-CLOSED-interval
// lower(it->first) <= x <= upper(it->first)
SIMI locate_left_iterator(const split_interval_map &imap, as_pos32 x)
{
	SIMI it = imap.lower_bound(ROI(as_pos32(x - 1), as_pos32(x)));
	if(it == imap.end() && it == imap.begin()) return it;
	if(it == imap.end()) it--;

	while(upper(it->first) > x)
	{
		if(it == imap.begin()) return imap.end();
		it--;
	}
	return it;
}

PSIMI locate_boundary_iterators(const split_interval_map &imap, as_pos32 x, as_pos32 y)
{
	SIMI lit, rit;
	lit = locate_right_iterator(imap, x);
	if(lit == imap.end() || upper(lit->first) > y) lit = imap.end();

	rit = locate_left_iterator(imap, y);
	if(rit == imap.end() || lower(rit->first) < x) rit = imap.end();

	if(lit == imap.end()) assert(rit == imap.end());
	if(rit == imap.end() && lit != imap.end()) 
	{
		// printf("x = %d, y = %d, lit = [%d, %d)\n", x, y, lower(lit->first), upper(lit->first));
		assert(lit == imap.end());
	}

	return PSIMI(lit, rit); 
}

int32_t compute_max_overlap(const split_interval_map &imap, SIMI &p, SIMI &q)
{
	if(p == imap.end()) return 0;

	int32_t s = 0;
	for(SIMI it = p; it != q; it++)
	{
		int32_t x = it->second;
		if(x > s) s = x;
	}

	if(q != imap.end())
	{
		int32_t x = q->second;
		if(x > s) s = x;
	}

	return s;
}

int32_t compute_sum_overlap(const split_interval_map &imap, SIMI &p, SIMI &q)
{
	if(p == imap.end()) return 0;

	int32_t s = 0;
	for(SIMI it = p; it != q; it++)
	{
		int l = lower(it->first);
		int u = upper(it->first);
		assert(u >= l);

		//printf(" AA add [%d, %d) : %d\n", lower(it->first), upper(it->first), it->second);

		s += (u - l) * it->second;
	}
	if(q != imap.end())
	{
		//printf(" BB add [%d, %d) : %d\n", lower(q->first), upper(q->first), q->second);
		s += (upper(q->first) - lower(q->first)) * q->second;
	}
	return s;
}

int evaluate_rectangle(const split_interval_map &imap, as_pos32 ll, as_pos32 rr, double &ave, double &dev, double &max)
{
	ave = 0;
	dev = 1.0;

	PSIMI pei = locate_boundary_iterators(imap, ll, rr); //typedef pair<SIMI, SIMI> PSIMI;
	SIMI lit = pei.first, rit = pei.second; // typedef split_interval_map::const_iterator SIMI;

	if(lit == imap.end()) return 0;
	if(rit == imap.end()) return 0;

	ave = 1.0 * compute_sum_overlap(imap, lit, rit) / (rr - ll);
	//printf("compute average %d-%d = %.2lf\n", ll, rr, ave);

	double var = 0;
	for(SIMI it = lit; ; it++)
	{
		assert(upper(it->first) > lower(it->first));
		var += (it->second - ave) * (it->second - ave) * (upper(it->first) - lower(it->first));

		if(it == rit) break;
	}

	dev = sqrt(var / (rr - ll));
	//if(dev < 1.0) dev = 1.0;
	
	max = 1.0 * compute_max_overlap(imap, lit, rit);
	return 0;
}

/**
 * find (x, y): returns first ROI with any overlapping, or end()
 * lowerbound(ROI(x, y)): returns first ROI strictly not less (ROI's upper > x) maybe with overlapping or end()
 * upperbound(ROI(x, y)): returns first ROI strictly larger (ROI's lower >= y) without overlapping or end()
 */
int test_split_interval_map()
{
	split_interval_map imap;

	imap += make_pair(ROI(7, 10), 3);
	imap += make_pair(ROI(6, 7), 3);
	imap += make_pair(ROI(1, 3), 3);
	imap += make_pair(ROI(1, 2), 1);
	imap += make_pair(ROI(2, 5), 2);

	SIMI it;
	
	for(it = imap.begin(); it != imap.end(); it++)
	{
		printf("interval: [%s,%s) -> %d\n", lower(it->first).aspos32string().c_str(), upper(it->first).aspos32string().c_str(), int32_t(it->second));
	}

	for(int i = 1; i <= 7; i++)
	{
		it = imap.find(ROI(i, i+1));
		if(it == imap.end())
		{
			printf("find %d-%d: does not exist\n", i,  i+1);
		}
		else
		{
			printf("find %d: [%s,%s) -> %d\n", i, lower(it->first).aspos32string().c_str(), upper(it->first).aspos32string().c_str(), it->second);
		}
	}

	for(int i = 1; i <= 7; i++)
	{
		it = imap.find(ROI(i, (i+2)));
		if(it == imap.end())
		{
			printf("find %d: does not exist\n", i);
		}
		else
		{
			printf("find %d: [%s,%s) -> %d\n", i, lower(it->first).aspos32string().c_str(), upper(it->first).aspos32string().c_str(), it->second);
		}
	}

	for(int i = 1; i <= 8; i++)
	{
		it = imap.lower_bound(ROI(i, i + 1));

		if(it == imap.end())
		{
			printf("lower bound %d: does not exist\n", i);
		}
		else
		{
			printf("lower bound %d: [%s,%s) -> %d\n", i, lower(it->first).aspos32string().c_str(), upper(it->first).aspos32string().c_str(), it->second);
		}
	}
	
	for(int i = 1; i <= 8; i++)
	{
		it = imap.lower_bound(ROI(i, i+2));

		if(it == imap.end())
		{
			printf("lower bound %d: does not exist\n", i);
		}
		else
		{
			printf("lower bound %d: [%s,%s) -> %d\n", i, lower(it->first).aspos32string().c_str(), upper(it->first).aspos32string().c_str(), it->second);
		}
	}

	for(int i = 1; i <= 8; i++)
	{
		it = imap.upper_bound(ROI(i, i + 1));

		if(it == imap.end())
		{
			printf("upper bound %d: does not exist\n", i);
		}
		else
		{
			printf("upper bound %d: [%s,%s) -> %d\n", i, lower(it->first).aspos32string().c_str(), upper(it->first).aspos32string().c_str(), it->second);
		}
	}

	for(int i = 1; i <= 8; i++)
	{
		it = imap.upper_bound(ROI(i, i+2));

		if(it == imap.end())
		{
			printf("upper bound %d: does not exist\n", i);
		}
		else
		{
			printf("upper bound %d: [%s,%s) -> %d\n", i, lower(it->first).aspos32string().c_str(), upper(it->first).aspos32string().c_str(), it->second);
		}
	}

	return 0;
}

