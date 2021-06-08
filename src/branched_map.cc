/*
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <iostream>
#include <map>
#include <vector>
#include <string>
// #include <tuple>
#include <cassert>
#include "as_pos32.hpp"
#include "branched_map.hpp"
#include "util.h"

using namespace std;

branched_map::branched_map()
{}

branched_map::~branched_map()
{}

branched_map branched_map::add(pair<as_pos32, as_pos32> pp, int i)
{
    // assert(validate_pexon(pp));
    if (!validate_pexon(pp)) throw BundleError();
    string ale = pp.first.ale;
    int32_t p1 = pp.first.p32;
    int32_t p2 = pp.second.p32;
    return this->add(p1, p2, ale, i);
}

branched_map branched_map::add(as_pos32 l, as_pos32 r, int i)
{
    return this->add(make_pair(l, r), i);
}

branched_map branched_map::add(int32_t p1, int32_t p2, string ale, int i)
{
    assert(p1<p2);
    pair<int32_t, int32_t> pp = make_pair(p1, p2);
    auto it = linear_interval.find(pp);
    if (it == linear_interval.end())
    {
        map <string, int> sim;
        sim.insert(make_pair(ale, i));
        linear_interval.insert(make_pair(make_pair(p1, p2), sim));
    }
    else
    {
        auto sim_it = (it->second).find(ale);
        if (sim_it == (it->second).end())
        {
            (it->second).insert(make_pair(ale, i));
        }
        else
        {
            sim_it->second = sim_it->second + i;
        }
    }
    return *this;    
}

branched_map branched_map::operator+= (pair< pair<as_pos32, as_pos32>, int > ppi)
{
    return this->add(ppi.first, ppi.second);
}

pair <pair<as_pos32, as_pos32>, int>* branched_map::find(as_pos32 l, as_pos32 r)
{
    assert(validate_pexon(make_pair(l, r)));
    string ale = l.ale;
    int32_t p1 = l.p32;
    int32_t p2 = r.p32;
    return this->find(p1, p2, ale);
}

pair <pair<as_pos32, as_pos32>, int>* branched_map::find(int32_t l, int32_t r, string ale)
{
    pair<int32_t, int32_t> lr = make_pair(l, r);
    auto up_it = linear_interval.upper_bound(lr);           // key greater than lr
    auto low_it = linear_interval.lower_bound(lr);          // key no larger than lr.

    // cout << "1. "<< l << "-" << r << ale << " " << low_it->first.first << low_it->first.second << endl;

    if (low_it == linear_interval.end() || low_it->first.first > l )
    {
        assert(low_it == up_it);
        --low_it;
        // cout << "2. " << l << "-" << r << ale << " " << low_it->first.first << low_it->first.second << endl;
    } 

    if (low_it == linear_interval.end()) return nullptr;
    if (!low_it->first.first <= l) return nullptr;
    // assert(low_it->first.second >= r); // low_it can be smaller or equal
    if (up_it != linear_interval.end())
    {
        assert(up_it->first.second >= r);
        assert(up_it->first.first >= r || up_it->first.first == l);
    }
    if(!(low_it->first.first == l || (low_it->first.second == r))) return nullptr;            // because what to find is a pexon //TODO: change this to assert after implementing DEL in spos&itvm
    if (ale == "$")
    {
        if (low_it->second.find(ale) == low_it->second.end())  ale = "~";
    }
    if(!(low_it->second.find(ale) != low_it->second.end())) return nullptr;                   // ale present // TODO: change this to assert after taking care of DEL and INS in itvm

    int c = (low_it->second)[ale];

    as_pos32 aslpos = as_pos32(low_it->first.first);
    as_pos32 asrpos = as_pos32(low_it->first.second);
    if (ale == "$") {aslpos.ale = "$"; asrpos.ale = "$";}
    else if (ale == "~") {aslpos.ale = "~"; asrpos.ale = "$";}
    else {aslpos.ale = ale; asrpos.ale = ale;}

    pair <pair<as_pos32, as_pos32>, int> ptr_con;
    ptr_con = make_pair(make_pair(aslpos, asrpos), c);
    pair <pair<as_pos32, as_pos32>, int> * ptr  = &ptr_con;
    return ptr;
}


void branched_map::clear()
{
    linear_interval.clear();
    assert(linear_interval.size() == 0);
}

bool branched_map::validate_pexon(pair<as_pos32, as_pos32> pp) 
{
    if (pp.first.ale == pp.second.ale)  return true;
    if (pp.first.ale == "~" || pp.second.ale == "$")  return true;
    cerr << "Error in making pmap! Not valid pexon: ";
    cerr << pp.first.p32 << pp.first.ale << "-";
    cerr << pp.second.p32 << pp.second.ale << endl;
    return false;
}

bool branched_map::is_no_overlap()
{
    int32_t prev_end = -1;
    for (auto i: linear_interval)
    {
        auto lr = i.first;
        int32_t l = lr.first;
        int32_t r = lr.second;
        if (prev_end > l) return false;
        if (l > r) return false; 
        prev_end = r;
    }
    return true;
}

void branched_map::print()
{
    cout << "pmap:\n Position\tAle=Idx\n";
    for (auto i: linear_interval)
    {
        cout << i.first.first << "\t" << i.first.second << "\t";
        for (auto j: i.second)
        {
            cout << j.first << "=" << j.second<< "\t";
        }
        cout << endl;
    }
}