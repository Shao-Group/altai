/*
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __AS_POS_HPP__
#define __AS_POS_HPP__

#include <string>

class as_pos 
{
    public:
        as_pos();
        as_pos(int64_t, std::string);
        as_pos(int32_t, int32_t, std::string);                   
        as_pos(const as_pos &a);     
        as_pos& operator=(const as_pos &a);              

    public:
        int64_t p64;
        std::string ale;

    public: 
        bool outside(as_pos a);
        bool outside_strict(as_pos a);
        bool inside(as_pos a);
        bool inside_strict(as_pos a);
        bool sameasitv(as_pos a);
        bool operator==(int _i)                          { return p64 == _i; }
        bool operator!=(int _i)                          { return p64 != _i; }
        bool operator==(as_pos _a)                       { return (p64 == _a.p64) && (ale == _a.ale); }
        bool operator==(const as_pos _a) const           { return (p64 == _a.p64) && (ale == _a.ale); }
        bool operator!=(as_pos _a)                       { return !((*this) == _a); }
        bool operator!=(const as_pos _a) const           { return !((*this) == _a); }
        bool operator<(int _i)                           { return p64 < _i; }
        bool operator<(int _i) const                     { return p64 < _i; }
        bool operator>(int _i)                           { return p64 > _i; }
        bool operator<=(int _i)                          { return (*this) < _i || (*this) == _i; }
        bool operator>=(int _i)                          { return (*this) > _i || (*this) == _i; }
        bool operator<(as_pos);
        bool operator<(const as_pos& _a) const;
        bool operator>(as_pos);
        bool operator<=(as_pos _a)                       { return (*this) < _a || (*this) == _a; }
        bool operator>=(as_pos _a)                       { return (*this) > _a || (*this) == _a; }
};

#endif 
