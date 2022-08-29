/*
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <as_pos.hpp>
#include <as_pos32.hpp>
#include <stdexcept>
#include "util.h"


as_pos32::as_pos32()
{
    p32 = -1;
    ale = "$";
}


as_pos32::as_pos32(const as_pos &as, bool isl, bool isr)
{
    if (isl == isr) throw std::invalid_argument("cannot be both/neither lpos && rpos") ;
    if (isl) 
    {
        p32 = (int32_t)((as.p64) >> 32);
    }
    if (isr) 
    {
        p32 = (int32_t)(((as.p64) << 32) >> 32);
    }
    ale = as.ale;
}

as_pos32::as_pos32(const int32_t i)
{
    p32 = i;
    ale = "$";
}

as_pos32::as_pos32(int32_t i, string s)
{
    p32 = i;
    ale = s;
}


as_pos32 as_pos32::operator=(const as_pos32 &a)
{
    p32 = a.p32;
    ale = a.ale;
    return *this;
}


bool as_pos32::operator==(as_pos32 a) const
{
    if (p32 == a.p32) 
    {
        // if (ale == "$" || a.ale == "$") return true;
        if (ale == a.ale) return true;
    }
    return false;
}


bool as_pos32::operator<(const as_pos32& a) const
{
    if (p32 < a.p32) return true;
    if ((p32 == a.p32) && (ale.compare(a.ale) < 0)) return true;
    return false;
}
       
bool as_pos32::operator>(as_pos32 a) const
{
    if (p32 > a.p32) return true;
    if ((p32 == a.p32) && (ale.compare(a.ale) > 0)) return true;
    return false;
}

string as_pos32::aspos32string() const
{
    string s = to_string(p32) + ale;
    return s;
}