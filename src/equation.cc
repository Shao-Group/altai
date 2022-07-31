/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "equation.h"
#include "util.h"
#include <cstdio>

equation::equation()
{
	e = 0;
	f = 0;
	a = 0;
	d = 0;
	w = 0;
}

equation::equation(double _e)
	:e(_e)
{
	f = 0;
	a = 0;
	d = 0;
	w = 0;
}

equation::equation(const vector<int> &_s, const vector<int> &_t)
	: s(_s), t(_t)
{
	e = 0;
	f = 0;
	d = 0;
	a = 0;
	w = 0;
}

equation::equation(const vector<int> &_s, const vector<int> &_t, double _e)
	: s(_s), t(_t), e(_e)
{
	f = 0;
	d = 0;
	a = 0;
	w = 0;
}

int equation::clear()
{
	s.clear();
	t.clear();
	e = 0;
	f = 0;
	a = 0;
	d = 0;
	w = 0;
	return 0;
}

int equation::print(int index) const
{
	printf("equation %3d: (%2lu, %2lu) edges, error = %.2lf, f = %d, w = %d, adjacent = %2d, distant = %2d. ", 
			index, s.size(), t.size(), e, f, w, a, d);

	printf("S = ( ");
	printv(s);
	printf("), T = ( ");
	printv(t);
	printf(")\n");
	
	return 0;
}
