/*
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <string>
#include <iostream>
#include "util.h"
#include "as_pos.hpp"
#include "as_pos32.hpp"

as_pos::as_pos(){}

as_pos::as_pos(int64_t p, std::string a) 
{
    p64 = p;
    ale = a;  // "$" means any string
}

as_pos::as_pos(int32_t p1, int32_t p2, std::string a) 
{
    p64 = pack(p1, p2);
    ale = a;  // "$" means any string
}

as_pos::as_pos(const as_pos &a)
{
    p64 = a.p64;
    ale = a.ale;
}


as_pos& as_pos::operator=(const as_pos &a)
{
    ale = a.ale;
    p64 = a.p64;
    return (*this);
}

bool as_pos::operator< (as_pos _a) 
{
    if (p64 < _a.p64) return true;
    if ((p64 == _a.p64) && (ale.compare(_a.ale) < 0)) return true; 
    return false;
}

bool as_pos::operator< (const as_pos& _a) const
{
    if (p64 < _a.p64) return true;
    if ((p64 == _a.p64) && (ale.compare(_a.ale) < 0)) return true; 
    return false;
}

bool as_pos::operator> (as_pos _a) 
{
    if (p64 > _a.p64) return true;
    if ((p64 == _a.p64) && (ale.compare(_a.ale) > 0)) return true; 
    return false;
}

bool as_pos::outside(as_pos a)                           { return high32(*this).leftsameto(high32(a)) && low32(*this).rightsameto(low32(a)); }
bool as_pos::outside_strict(as_pos a)                    { return high32(*this).leftto(high32(a)) && low32(*this).rightto(low32(a)); }
bool as_pos::inside(as_pos a)                            { return high32(*this).rightsameto(high32(a)) && low32(*this).leftsameto(low32(a)); }
bool as_pos::inside_strict(as_pos a)                     { return high32(*this).rightto(high32(a)) && low32(*this).leftto(low32(a)); }

bool as_pos::sameasitv(as_pos a)                         { return a.p64 == (this->p64);}
