/*
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __AS_POS32_HPP__
#define __AS_POS32_HPP__

#include <string>
#include "as_pos.hpp"

using namespace std;

class as_pos32
{
    public:
        as_pos32();
        as_pos32(const as_pos&, bool, bool);
        as_pos32(const int32_t);
        as_pos32(int32_t, std::string);
        as_pos32 operator=(const as_pos32 &a);        

    public:
        int32_t p32;                                        // position on ref 
        std::string ale;                                    // alias of sequence in string vector, to the right of the pos p32
    
    public:
        bool leftto(const as_pos32 a) const                 { return (this->p32 < a.p32); }
        bool rightto(const as_pos32 a) const                { return (this->p32 > a.p32); }
        bool samepos(const as_pos32 a) const                { return (this->p32 == a.p32); }
        bool leftsameto(const as_pos32 a) const             { return (this->leftto(a) || this->samepos(a)); }
        bool rightsameto(const as_pos32 a) const            { return (this->rightto(a) || this->samepos(a)); }
        int operator+(int i) const                          { return int(p32+i); }
        int operator+(as_pos32 as) const                    { return int(p32+as.p32); }
        int operator-(int i) const                          { return int(p32-i); }
        double operator-(double d) const                    { return double(p32) - d; }
        int operator-(as_pos32 a) const                     { return int(p32 - a.p32); }
        bool operator==(int _i) const                       { return p32 == _i; }
        bool operator!=(int _i) const                       { return p32!= _i; }
        bool operator==(as_pos32) const;                    // "$" represents any seq
        bool operator!=(as_pos32 a) const                   { return !((*this)==a); }
        bool operator<(int i) const                         { return this->p32 < i; }
        bool operator>(int i) const                         { return p32 > i; }
        bool operator<=(int i) const                        { return p32 <= i; }
        bool operator>=(int i) const                        { return p32 >= i; }
        bool operator<(const as_pos32& _a) const;
        bool operator>(as_pos32) const;
        bool operator<=(as_pos32 a) const                   { return (*this)< a || (*this) ==a; }
        bool operator>=(as_pos32 a) const                   { return (*this)> a || (*this) ==a; }     
        operator int32_t() const                            { return p32;}

        string aspos32string() const; 
        static bool inside_strict(as_pos32 l1, as_pos32 r1, as_pos32 l2, as_pos32 r2); // true if l2/r2 is strictly inside l1/r1
};


inline as_pos32 high32 (as_pos _a) { return as_pos32(_a, true, false);}
inline as_pos32 low32 (as_pos _a) { return as_pos32(_a, false, true);}

inline int32_t high32(int64_t x) { return (int32_t)((x) >> 32); }
inline int32_t low32(int64_t x) { return (int32_t)(((x) << 32) >> 32); }


#endif 
