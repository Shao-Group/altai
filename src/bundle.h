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
	bundle_base &bb;				// input bundle base	
	bundle_bridge br;				// contains fragments
	split_interval_map fmap;		// matched interval map
	vector<junction> junctions;													// splice junctions
	vector<junction> allelic_junctions; 										// allelic junctions
	map<as_pos, vector<int> > allelic_itv;										// allelic aspos intervals and hits containing them
	vector<region> regions;														// regions
	vector<partial_exon> pexons;												// partial exons
	vector<bool> regional;														// if a pe is regional
	map<pair<as_pos32, as_pos32>, int> pmap;									// partial exon map, save pexon and its index
	map<pair<as_pos32, as_pos32>, int> pmap_na;									// non allelic partial exon map, save pexon and its index
	map<pair<as_pos32, as_pos32>, int> pmap_a;									// allelic partial exon map, save pexon and its index
	splice_graph gr;															// splice graph
	hyper_set hs;																// hyper set

public:
	virtual int build(int mode, bool revise);
	// int count_junctions() const;
	int print(int index);

public:
	int prepare();
	// check and init
	int check_left_ascending();
	int check_right_ascending();
	int compute_strand();

	// splice graph
	int build_intervals();
	int build_junctions();
	int build_allelic_junctions();
	int build_regions();
	int build_partial_exons();
	vector<partial_exon> partition_allelic_partial_exons(const partial_exon&);
	int link_partial_exons();
	// int build_splice_graph();
	int build_splice_graph(int mode);
	int build_partial_exon_map();
	int locate_left_partial_exon(as_pos32 x);
	int locate_right_partial_exon(as_pos32 x);
	vector<int> align_hit(hit &h);
	vector<int> align_fragment(fragment &f);

	// revise splice graph
	VE compute_maximal_edges();
	int revise_splice_graph();
	int refine_splice_graph();
	bool keep_surviving_edges();
	bool extend_boundaries();
	bool extend_start_boundaries();
	bool extend_end_boundaries();
	bool remove_small_junctions();
	bool remove_small_exons();
	bool remove_inner_boundaries();
	bool remove_intron_contamination();
	bool remove_false_boundaries();
	bool tackle_false_boundaries();
	// int find_contamination_chain();

	// super edges
	// int build_hyper_edges2();			// paired end
	// bool bridge_read(int x, int y, vector<int> &s);

	// hyper set
	int build_hyper_set();
};

#endif
