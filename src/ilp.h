#ifndef __ILP_H__
#define __ILP_H__

#include "gurobi_c++.h"
#include "hyper_set.h"
#include "splice_graph.h"
#include "path.h"

#include <vector>
#include <set>
#include <map>

using namespace std;

class ilp
{
public:
    ilp(const splice_graph &_gr, const hyper_set &_HS, GRBEnv *e);
    virtual ~ilp();
    
public:
    GRBEnv * env;
    GRBModel * model;
    
    splice_graph gr;
    hyper_set HS;
    vector < vector <int> > H;
    vector < vector <int> > Hv;
    vector < int > H_weight;
    MEI e2i;// edge map, from edge to index
    VE i2e;// edge map, from index to ed
    map< int, set<int> > map_e_to_phasing;
    
    int M; //phasing path
    int N;//result path
    int N1;
    int N2;
    double MAX_EDGE_WEIGHT;
    int FLOW_NUM;
    double PHASING_PATH_THRESHOLD;
    
    vector< vector<GRBVar> > X;
    vector< vector<GRBVar> > T;
    vector< GRBVar > EDGE_SUPPORT;
    vector< GRBVar > VERTEX_SUPPORT;

    vector< vector<GRBVar> > Y;
    vector< vector<GRBVar> > S;
    vector< GRBVar > Z;
    vector< GRBVar > PHASING_SUPPORT;
    vector< vector< vector<GRBVar> > > UBJ;
    vector< vector<GRBVar> > D; 
    vector< vector<GRBVar> > NPE;

    vector<path> paths;//results
    set<int> AntiChain;
    set<int> POSSIBLE_INTRON;
    int INTRON_PATH;
    int MINOR_PATH;
    double PHASING_COVER_RATIO;
    set<int> INTRON_PHASING_PATH;
    vector < vector <int> > H_conflict;
    vector < vector <int> > H_X_conflict;
    map<int, int> H_nest;
    vector<double> SQRTDIFF;
    vector<double> SQRTDIFF_E;
    set<int> NON_OVLP_EDGE;

    vector< GRBVar > Ratio;
    set<int> minor_phasing_path;
    int PHASING_COUNT;
    double PHASING_WEIGHT_RATIO;
    double PHASING_COUNT_RATIO;

    GRBVar obj_UBJ;
    GRBVar obj_edge_weight;
    GRBVar obj_phasing_count;
    GRBVar obj_phasing_cover;
    
private:
    directed_graph bgr;        // bipartite graph
    int max_antichain();
    int greedy_antichain();
    int max_weight_antichain();

    int remove_nest_phasing_path();
    int remove_possible_intron();
    int add_edge_not_in_H();
    int set_one_spot_phasing_path(int pred, int pid);
    bool whether_compatible(vector<int> &h1, vector<int> &h2);
    int draw_bipartite_graph(const string &file);
    bool sort_pair_decreasing(pair<int, int> p1, pair<int, int> p2);
    int sqrtdiff_trivial_vertex();
    int edge_start_mid_end();
    int topological_balance();

public:
    
    int solve(int suggestN);

    int add_edge_support_variables();
    int add_edge_abundance_variables();
    int add_predicted_edge_coverage_variables();
    int add_predicted_vertex_coverage_variables();
    
    int add_phasing_support_variables();
    int add_path_abundance_variables();
    int add_predicted_path_coverage_variables();
    int add_phasing_cover_variables();

    int add_valid_path_constraint();
    int add_valid_path_constraint_st();
    int add_x_to_y_constraint();
    int add_phasing_cover_constraint();
    //int add_anti_chain_constraint();
    int add_one_spot_constraint();
    int add_edge_conflict_constraint();
    int add_conflict_nest_constraint();
    
    int add_t0x0_constraint();
    int add_x_to_t_constraint();
    
    int add_s0y0_constraint();
    int add_y_to_s_constraint();
    int add_s_to_t_constraint();

    int add_unbalanced_junction_variables();
    int add_unbalanced_junction_constraint();

    int add_duplicate_path_variables();
    int add_duplicate_path_constraint();

    int add_non_phasing_edge_variables();

    int add_adjacent_edge_constraint();

    int allele_phase_constraint();

    int set_objective();
    int collect_results(int suggestN);
    int print_phasing_path(const vector<int> &phasing_edge);
};

#endif
