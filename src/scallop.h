/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __SCALLOP3_H__
#define __SCALLOP3_H__

#include "splice_graph.h"
#include "hyper_set.h"
#include "equation.h"
#include "router.h"
#include "path.h"
#include "transcript.h"

typedef map< edge_descriptor, vector<int> > MEV;
typedef pair< edge_descriptor, vector<int> > PEV;
typedef pair< vector<int>, vector<int> > PVV;
typedef pair<PEE, int> PPEEI;
typedef map<PEE, int> MPEEI;
typedef pair<int, int> PI;
typedef map<int, int> MI;

// for noisy splice graph
class scallop
{
public:
	scallop(splice_graph &gr, const hyper_set &hs, bool r = false, bool keep_as = true);                      // sc w. both alleles 
	scallop(splice_graph* gr, const hyper_set &hs, const scallop &sc, bool r = false, bool keep_as = false);  // sc for each allele
	virtual ~scallop();

public:
	int assemble(bool is_allelic);
	int assemble_continue(bool is_allelic);
	int transform(splice_graph* pgr, const VE& old_i2e, const MEE& x2y);  // allelic transform

public:
	splice_graph& gr;					// splice graph
	MEI e2i;							// edge map, from edge to index
	VE i2e;								// edge map, from index to edge
	bool random_ordering;				// whether using random ordering
	MEV mev;							// super edges
	vector<int> v2v;					// vertex map
	hyper_set hs;						// hyper edges
	int round;							// iteration
	bool keep_as_nodes;					// true: decompose ns nodes only, false: decompose all nodes

	set<int> asnonzeroset;			    // vertices with degree >= 1 && !vi.is_ordinary_vertex()
	set<int> nsnonzeroset;				// vertices with degree >= 1 &&  vi.is_ordinary_vertex()
	vector<path> paths;					// predicted paths
	vector<transcript> trsts;			// predicted transcripts
	vector<transcript> non_full_trsts;		// predicted non full length transcripts

private:
	// init
	int classify();
	int init_vertex_map();
	int init_super_edges();
	int init_inner_weights();
	int init_vertex_astype();
	int init_nonzeroset(bool keep_as_nodes);

	// resolve iteratively
	bool resolve_trivial_vertex(int type, double jump_ratio);
	bool resolve_trivial_vertex_fast(double jump_ratio);
	bool resolve_single_trivial_vertex_fast(int i, double jump_ratio);
	bool resolve_smallest_edges(double max_ratio);
	bool resolve_negligible_edges(bool extend, double max_ratio);
	bool resolve_splittable_vertex(int type, int degree, double max_ratio);
	bool resolve_unsplittable_vertex(int type, int degree, double max_ratio);
	bool resolve_hyper_edge(int fsize);

	// smooth vertex
	int balance_vertex(int x);
	double compute_balance_ratio(int x);

	// decomposing subroutines
	int compute_smallest_edge(int x, double &ratio);
	int decompose_trivial_vertex(int v);
	int decompose_vertex_extend(int v, MPID &pe2w);
	int decompose_vertex_replace(int v, MPID &pe2w);
	int classify_trivial_vertex(int v, bool fast);
	int exchange_sink(int old_sink, int new_sink);
	int split_vertex(int x, const vector<int> &xe, const vector<int> &ye);
	int split_edge(int exi, double w);
	int merge_adjacent_edges(int x, int y);
	int merge_adjacent_edges(int x, int y, double ww);
	int merge_adjacent_equal_edges(int x, int y);
	int remove_edge(int e);
	int split_merge_path(const VE &p, double w);
	int split_merge_path(const vector<int> &p, double w);
	int collect_existing_st_paths();
	int collect_path(int e);
	int compute_length(const path &p);
	int greedy_decompose();

	// debug
	bool assert_debug();
	bool assert_mev_gr_edge_descriptor_bijection();
	bool assert_mev_super_set_gr_edge_descriptor();

	// stats, print, and draw
	int print();
	int stats();
	int draw_splice_graph(const string &file);
	vector<int> topological_sort();
	int clear();
};

#endif
