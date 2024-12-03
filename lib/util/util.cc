/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <string>
#include <cstring>
#include "util.h"

string toupperstring(const string s)
{
	string s2 = "";
	for (int i = 0; i < s.length(); i++)
	{
		s2 = s2 + char(toupper(s[i]));
	}
	return s2;
}

// mimic python string.split() 
vector<string> split(const string s, const string sep)
{
	vector<string> v;
	string s2 = s;
	int delim_pos = 0;
	while((delim_pos = s2.find(sep)) != string::npos) 
	{
		v.push_back(s2.substr(0, delim_pos));
		s2.erase(0, delim_pos + sep.length()); //also remove sep
	}
	v.push_back(s);
	return v;
}

size_t string_hash(const std::string& str)
{
	size_t hash = 1315423911;
	for(std::size_t i = 0; i < str.length(); i++)
	{
		hash ^= ((hash << 5) + str[i] + (hash >> 2));
	}

	return (hash & 0x7FFFFFFF);
}

// size_t vector_hash(const vector<int32_t> & vec)
// {
// 	size_t seed = vec.size();
// 	for(int i = 0; i < vec.size(); i++)
// 	{
// 		seed ^= (size_t)(vec[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
// 	}
// 	return (seed & 0x7FFFFFFF);
// }

size_t vector_hash(const vector<as_pos32> & vec)
{
	size_t seed = vec.size();
	string s = "";
	for(int i = 0; i < vec.size(); i++)
	{
		seed ^= (size_t)(vec[i].p32 + string_hash(vec[i].ale)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		s += vec[i].ale;
	}
	return (seed & 0x7FFFFFFF);
}


// reverse-complement DNA sequence w. IUPAC symbols
int reverse_complement_DNA(string &rc, const string s)
{
	char c;
	for(int i = s.length() - 1; i >= 0; --i)
	{
		switch (toupper(s[i]))
		{
		case 'A':
			c = 'T';
			break;
		case 'T':
			c = 'A';
			break;
		case 'C':
			c = 'G';
			break;
		case 'G':
			c = 'C';		
			break;
		case 'U':
			c = 'A';
			break;
		case 'N':
			c = 'N';
			break;
		case 'W':
			c = 'W';
			break;
		case 'S':
			c = 'S';
			break;
		case 'M':
			c = 'K';
			break;
		case 'K':
			c = 'M';
			break;
		case 'R':
			c = 'Y';
			break;
		case 'B':
			c = 'V';
			break;
		case 'D':
			c = 'H';
			break;
		case 'H':
			c = 'D';
			break;
		case 'V':
			c = 'B';
			break;
		case '-':
			c = '-';
			break;
		default:
			c = 'N';
			break;
		}
		rc += c;
	}
	return 0;
}