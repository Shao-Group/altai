/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __BUNDLE_H__
#define __BUNDLE_H__

#include "interval_map.h"
#include "bundle_base.h"
#include "bundle_bridge.h"
#include "junction.h"
#include "region.h"
#include "partial_exon.h"
#include "splice_graph.h"
#include "hyper_set.h"
#include "path.h"
#include "gene.h"
#include "transcript.h"

using namespace std;

class bundle
{
public:
	bundle(bundle_base &bb);
	virtual ~bundle();

public:
	bundle_base &bb;															// input bundle base	
	bundle_bridge br;															// contains fragments
	split_interval_map fmap;													// matched interval map, not AS. (alleles collapsed)

	/* 
	** re-use from bundle_bridge
	** vector<junction> junctions;													// splice junctions, not needed, use h.get_aligned_interval
	** vector<junction> allelic_junctions; 										    // allelic junctions, not used
	** map<as_pos, vector<int> > allelic_itv;										// allelic aspos intervals and hits containing them
	** vector<region> regions;														// regions, use br.regions instead
	*/
	
	vector<partial_exon> pexons;												// partial exons
	map<pair<int32_t, int32_t>, vector<int> > pos_pids;							// pos pair to partial exon ids, allelic pexons are put in vector
	vector<bool> regional;														// if a pe is regional
	map<pair<int, int>, int > jset;											    // < <pid-to-pid, hit-counts>
	splice_graph gr;															// splice graph
	hyper_set hs;																// hyper set

public:
	virtual int build(int mode, bool revise);
	int print(int index);

public:
	int prepare();
	// check and init
	int check_left_ascending();
	int check_right_ascending();
	int compute_strand();

	// splice graph
	int build_intervals();
	int build_partial_exons();
	int build_pos_pids_map();
	int build_pseudo_variant_exon();
	int pexon_jset(map<pair<int, int>, int >& pexon_jset);
	int build_splice_graph(int mode);

	// revise splice graph 
	int revise_splice_graph();
	bool remove_false_boundaries();
	bool tackle_false_boundaries();
	
	// hyper set
	int build_hyper_set();

private:
	int build_splice_graph_vertices(int mode);
	int build_splice_graph_edges(int mode);
	int add_pseudo_as_in_edge(int mode, int pse_id, int counter_v_id);
	int add_pseudo_as_out_edge(int mode, int pse_id, int counter_v_id);
	int build_splice_graph_vertices_as_type(int mode);
	vector<int> align_hit(hit &h);
	vector<int> align_fragment(fragment &f);
};

#endif
