/*
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __INTERVAL_2D_MAP_HPP__
#define __INTERVAL_2D_MAP_HPP__

#include <map>
#include <vector>
#include <string>
// #include <tuple>
#include "as_pos32.hpp"

using namespace std;
 

class branched_map
{
public:
    branched_map();
    ~branched_map();

public:
    map<pair<int32_t, int32_t>, map<string, int> > linear_interval;

public:
    branched_map add(pair<as_pos32, as_pos32>, int);
    branched_map add(as_pos32 l, as_pos32 r, int);
    branched_map add(int32_t, int32_t, string, int);
    branched_map operator+= (pair< pair<as_pos32, as_pos32>, int >);


    pair <pair<as_pos32, as_pos32>, int>* find(as_pos32, as_pos32);   // return count
    pair <pair<as_pos32, as_pos32>, int>* find(int32_t, int32_t, string);
    map<pair<int32_t, int32_t>, map<string, int> >::iterator end()  {return linear_interval.end();}
    void clear();
    bool is_no_overlap();
    void print();

private:
    bool validate_pexon(pair<as_pos32, as_pos32> pp);

};

#endif
