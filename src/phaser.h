/*
Part of Altai
(c) 2022 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __PHASER_H__
#define __PHASER_H__

#include <vector>
#include <map>
#include "splice_graph.h"
#include "hyper_set.h"
#include "as_pos32.hpp"
#include "bundle.h"
#include "scallop.h"
#define MEPD map<edge_descriptor, pair<double, double> >

/*
*   phaser takes scallop object as an input and does:
*   1. split it into two allelic scallop instances (objects splitted: sc, gr, hs);
*   2. assemble allelic scallop instances (sc1, sc2)
*       2.1 assemble anyway regardless of having only one allele or two. The unexpressed allele will be filtered by filter class.
*   3. transcripts and non-full-length transcripts are stored as public, and will be collected + processed in assembler class
*/
class phaser
{
public:
	phaser(scallop& sc, bool is_allelic);

private:
    const scallop& sc;
    splice_graph& gr;
    MED ewrt1;              // edge weight in allele1
    MED ewrt2;              // edge weight in allele2
    vector<double> vwrt1;   // vertex weight in allele1
    vector<double> vwrt2;   // vertex weight in allele2

    double vwrtbg1;        // sum bg weights of allele 1
    double vwrtbg2;        // sum bg weights of allele 2
    double ewrtbg1;        // normalized bg ratio of allele 1
    double ewrtbg2;        // normalized bg ratio of allele 2
    double ewrtratiobg1;   // normalized bg ratio of allele 1
    double ewrtratiobg2;   // normalized bg ratio of allele 2

    MEE x2y_1;             // use x2y to map original edge to new edge, in allele 1
	MEE y2x_1;             // use y2x to map new edge to original edge, in allele 1
    MEE x2y_2;             // use x2y to map original edge to new edge, in allele 2
	MEE y2x_2;             // use y2x to map new edge to original edge, in allele 2

    splice_graph* pgr1;    // pointer to sg of allele1
    splice_graph* pgr2;    // pointer to sg of allele2
    hyper_set* phs1;       // pointer to hs of allele1
    hyper_set* phs2;       // pointer to hs of allele2

    bool is_allelic;       // should be true

public: 
    vector<transcript> trsts1;              // transcripts of allele1
    vector<transcript> trsts2;              // transcripts of allele2
    vector<transcript> non_full_trsts1;		// predicted non full length transcripts of allele1
    vector<transcript> non_full_trsts2;		// predicted non full length transcripts of allele2

private:
    string strategy;                    
    double epsilon = 0.01;  

private:
    int init();
    int assign_gt();
    int split_gr();
    int refine_allelic_graphs();
    int split_hs();
    int assemble_allelic_scallop();     
    int allelic_transform(scallop& sc, splice_graph* pgr, MEE& x2y);
    int assign_transcripts_gt();

private:
    pair<double, double> get_as_ratio(int i);
    vector<int> sort_nodes_by_currecnt_mae(const set<int>& s);
    bool split_local(int i);
    bool split_global(int i);
    bool split_by_ratio(int v, const PEEI& in, const PEEI& out, double r1);                  // split edges of vertex v, by ratio 
    int split_by_phasing(int v, const PEEI& in, const PEEI& out, double r1);                // split edges of vertex v, by phasing path
    int split_by_min_parsimony(int v, const PEEI& in, const PEEI& out, double r1);          // split edges of vertex v, by parsimony
    pair<double, double> normalize_epsilon(double x, double y);                             // adjustment to allele ratio. new_r1 = (r1+eps) / (r1+r2+2*eps). returns<-1, -1> if both input are 0

};

#endif
