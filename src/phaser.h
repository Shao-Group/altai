/*
Part of Altai
(c) 2022 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __PHASER_H__
#define __PHASER_H__

#include <vector>
#include "splice_graph.h"
#include "hyper_set.h"
#include "as_pos32.hpp"
#include "bundle.h"

class phaser
{
public:
	phaser();


private:
    splice_graph gr;
    hyper_set hs;
    bundle bd;
    string chrm;
    vector<int> var_sites;  // pos of sites with hetero variants
    vector<int> allele1;    // variants on allele 1 <ref 0, alt 1, 2, ...>
    vector<int> allele2;    // variants on allele 2 <ref 0, alt 1, 2, ...>
                            // allele1.size() == allele2.size() == var_sites.size()
                            
private:
    int phase(const splice_graph gr, const hyper_set hs, const bundle bd);
    int clear();
};
#endif
