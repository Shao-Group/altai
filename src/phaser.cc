/*
Part of Altai
(c) 2022 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include "phaser.h"
#include "splice_graph.h"
#include "hyper_set.h"
#include "vertex_info.h"
#include "as_pos32.hpp"
#include <limits.h>

phaser::phaser(scallop& _sc, splice_graph* _gr1, hyper_set* _hs1, splice_graph* _gr2, hyper_set* _hs2)
	: sc(_sc), gr(_sc.gr), pgr1(_gr1), pgr2(_gr2), phs1(_hs1), phs2(_hs2)
{
	assert(sc.asnonzeroset.size() != 0); // throw runtime_error("does not have AS nodes");
	
	init();
	assign_gt();
	split_gr();
	refine_allelic_graphs();
	split_hs_by_rebuild();
	// refine_hyper_sets();
}

// init ewrt1/2, countbg1/2, normalize ratiobg1/2
int phaser::init()
{
	pgr1->clear();
	pgr2->clear();
	phs1->clear();
	phs2->clear();

	strategy = "split_by_ratio";

	for(auto i: gr.ewrt)
	{
		edge_descriptor e = i.first;
		ewrt1.insert({e, -1});
		ewrt2.insert({e, -1});
	}

	vwrtbg1 = 0;
	vwrtbg2 = 0;
	ewrtbg1 = 0;       
    ewrtbg2 = 0;
	for(int i = 0; i < gr.vinf.size(); i++)
	{
		if (gr.vinf[i].gt == ALLELE1)
		{
			PEEI in = gr.in_edges(i);
			PEEI out = gr.out_edges(i);
			for (auto e = in.first; e!= in.second; e++)	
			{
				ewrt1[*e] = gr.ewrt[*e];
				ewrt2[*e] = 0;
				ewrtbg1 += gr.ewrt[*e];
			}
			for (auto e = out.first; e!= out.second; e++)	
			{
				ewrt1[*e] = gr.ewrt[*e];
				ewrt2[*e] = 0;
				ewrtbg1 += gr.ewrt[*e];
			}
			vwrtbg1 += gr.get_vertex_weight(i);
		}
		else if (gr.vinf[i].gt == ALLELE2)
		{
			PEEI in = gr.in_edges(i);
			PEEI out = gr.out_edges(i);
			for (auto e = in.first; e!= in.second; e++)	
			{
				ewrt1[*e] = 0;
				ewrt2[*e] = gr.ewrt[*e];
				ewrtbg2 += gr.ewrt[*e];
			}
			for (auto e = out.first; e!= out.second; e++)	
			{
				ewrt1[*e] = 0;
				ewrt2[*e] = gr.ewrt[*e];
				ewrtbg2 += gr.ewrt[*e];
			}
			vwrtbg2 += gr.get_vertex_weight(i);
		}
	}
	ewrtratiobg1 = normalize_epsilon(ewrtbg1, ewrtbg1);
	ewrtratiobg2 = 1 - ewrtratiobg1;
	
	//TODO: what if only one allele is expressed?
	// if (countbg1 == 0 || countbg2 == 0){;} {	// return one empty graph}	

	return 0;
}

// assign edges to different gt
int phaser::assign_gt()
{
	set<int> asnodes;  // := as nodes only
	set<int> nsnodes;  // := ns & ns-related nodes only
	for(const int i: sc.asnonzeroset)
	{
		if(gr.vinf[i].is_as_vertex()) 
		{
			asnodes.insert(i);
		}
		else 
		{
			nsnodes.insert(i);
		}
	}
	for(const int i: sc.nsnonzeroset) 
	{
		nsnodes.insert(i);
	}
	assert(asnodes.size() >= 1);
	
	// rank ns nodes by as_neightbor_degree and % decomposed neighbors, split_local if small degree 
	if (nsnodes.size() <= max_num_exons)
	{
		int as_neighbor_degree = 1;
		bool is_traversed = false;
		set<int> dcnodes; 								// := splitted nodes
		set<int> nsnodes_with_last_degree{asnodes};		// := start from as nodes neighbors, deg 1

		int flag = 0; // debug print

		while(nsnodes.size() >= 1 && !is_traversed && as_neighbor_degree <= 3)
		{
			set<int> nsnodes_with_this_degree;
			nsnodes_with_this_degree.clear();
			for(int i: nsnodes_with_last_degree)
			{	
				const PEEI in = gr.in_edges(i);
				const PEEI out = gr.out_edges(i);
				for (auto e = in.first; e!= in.second; e++)	nsnodes_with_this_degree.insert((*e)->source());
				for (auto e = out.first; e!= out.second; e++)	nsnodes_with_this_degree.insert((*e)->target());
			}

			if (flag == 0 && DEBUG_MODE_ON)
			{
				cout << "flag" << flag << endl;
				cout << "nsnodes_with_this_degree" ;
				for (int k : nsnodes_with_this_degree) cout << k << " ";
				cout << endl;
				cout << "nsnodes_with_this_degree ";
				for (int k : nsnodes_with_this_degree) cout << k << " ";
				cout << endl;
				flag ++;
			}

			// remove ns node w. current deg and highest decomposed rate of neighbors
			nsnodes_with_last_degree = nsnodes_with_this_degree;
			bool removed_one = true;
			while(nsnodes_with_this_degree.size() >= 2 && removed_one) 
			{
				assert(as_neighbor_degree <= 3);
				removed_one = false;

				int num_as_neighrbor = -1;
				int num_dc_neighbor = -1;
				double ratio_dc_neighbor = -1;
				int node_to_split = -1;
				for(int i: nsnodes_with_this_degree)
				{	
					const PEEI in = gr.in_edges(i);
					const PEEI out = gr.out_edges(i);
					set<int> neighbors;
					for (auto e = in.first; e!= in.second; e++)	neighbors.insert((*e)->source());
					for (auto e = out.first; e!= out.second; e++)	neighbors.insert((*e)->target());
					int a = 0;
					int b = 0;
					for (int j: neighbors)
					{
						if(asnodes.find(j) != asnodes.end()) 
						{
							a++;
						}
						else if (dcnodes.find(j) != dcnodes.end()) 
						{
							b++;
						}
						b = b + a;
					}
					double c = double(b) / neighbors.size();

					if  (c > ratio_dc_neighbor ||
						(c >= ratio_dc_neighbor && a > num_as_neighrbor) || 
					    (c >= ratio_dc_neighbor && a >= num_as_neighrbor && b > num_dc_neighbor))
					{
						num_as_neighrbor = a;
						num_dc_neighbor = b;
						ratio_dc_neighbor = c;
						node_to_split = i;
					}
				}

				if (node_to_split!= -1)	
				{
					split_local(node_to_split);
					nsnodes_with_this_degree.erase(node_to_split);
					nsnodes.erase(node_to_split);
					dcnodes.insert(node_to_split);
					removed_one = true;
					is_traversed = false;
				}
			}

			if (nsnodes_with_this_degree.size() >= 1)	
			{
				assert(nsnodes_with_this_degree.size() == 1);
				int i = *nsnodes_with_this_degree.begin();
				split_local(i);
				nsnodes.erase(i);
				dcnodes.insert(i);
				is_traversed = false;
			}

			as_neighbor_degree += 1;
		}
	}

	// split all based on background, if as_neighbor_degree large or too many nodes
	for (int i: nsnodes)
	{
		split_global(i);
	}
	
	return 0;
}

int phaser::split_local(int i)
{
	// compute ratio
	vertex_info v = gr.get_vertex_info(i);
	const PEEI in = gr.in_edges(i);
	const PEEI out = gr.out_edges(i);

	double local1, local2;
	for (auto e = in.first; e!= in.second; e++)	
	{
		if(ewrt1[*e] > 0)	local1 += ewrt1[*e];
		if(ewrt2[*e] > 0)	local2 += ewrt2[*e]; 
	}
	for (auto e = out.first; e!= out.second; e++)
	{
		if(ewrt1[*e] > 0)	local1 += ewrt1[*e];
		if(ewrt2[*e] > 0)	local2 += ewrt2[*e]; 
	}

	double ratio1 = normalize_epsilon(local1, local2);
	double ratio2 = 1 - ratio1;
	
	if (strategy == "split_by_ratio")
	{
		split_by_ratio(i, in, out, ratio1);
	}
	else if (strategy == "split_by_phasing")
	{
		split_by_phasing(i, in, out, ewrtratiobg1);
	}
	else if (strategy == "min_parsimony")
	{
		split_by_min_parsimony(i, in, out, ratio1);
	}
	else
	{
		cerr << "split strategy = " << strategy << endl;
		throw runtime_error("split strategy not defined?");
	}

	return 0;
}

int phaser::split_global(int i)
{
	vertex_info v = gr.get_vertex_info(i);
	const PEEI in = gr.in_edges(i);
	const PEEI out = gr.out_edges(i);
	if (strategy == "split_by_ratio")
	{
		split_by_ratio(i, in, out, ewrtratiobg1);
	}
	else if (strategy == "split_by_phasing")
	{
		split_by_phasing(i, in, out, ewrtratiobg1);
	}
	else if (strategy == "min_parsimony")
	{
		split_by_min_parsimony(i, in, out, ewrtratiobg1);
	}
	else
	{
		cerr << "split strategy = " << strategy << endl;
		throw runtime_error("split strategy not defined?");
	}
	return 0;
}


//	split edges of vertex v, by ratio.
//  directly modify ewrt1, ewrt2, if ewrt1[e] or ewrt2[e] <= -1 (not assigned)
int phaser::split_by_ratio(int v, const PEEI& in, const PEEI& out, double ratio_allele1)
{
	assert(ratio_allele1 > 0); // ratio normalized, won't equal
	assert(ratio_allele1 < 1);
	vwrt1[v] = gr.get_vertex_weight(v) * ratio_allele1;
	vwrt2[v] = gr.get_vertex_weight(v) * (1 - ratio_allele1);
	for (auto e = in.first; e!= in.second; e++)	
	{
		double w = gr.ewrt[*e];
		if(ewrt1[*e] < 0)
		{
			ewrt1[*e] = w * ratio_allele1;
		} 
		if(ewrt2[*e] < 0)
		{
			ewrt2[*e] = w * (1 - ratio_allele1);
		} 
	}
	for (auto e = out.first; e!= out.second; e++)
	{
		double w = gr.ewrt[*e];
		if(ewrt1[*e] < 0)
		{
			ewrt1[*e] = w * ratio_allele1;
		} 
		if(ewrt2[*e] < 0)
		{
			ewrt2[*e] = w * (1 - ratio_allele1);
		}
	}
	
	return 0;
}

// split edges of vertex v, by phasing path
int phaser::split_by_phasing(int v, const PEEI& in, const PEEI& out, double r1)
{
	throw runtime_error("split_by_phasing not implemented yet");
	return 0;
}

int phaser::split_by_min_parsimony(int v, const PEEI& itr_in_edges, const PEEI& itr_out_edges, double ratio_allele1)
{
	throw runtime_error("split_by_parsimony not defined yet");
	return -1;
}

// split sg into two pairs of sg1/hs1 and sg2/hs2
int phaser::split_gr()
{
	splice_graph tmp_sg1(gr);

	tmp_sg1.vwrt = vwrt1;
	tmp_sg1.ewrt = ewrt1;
	MEE x2y_1;
	MEE y2x_1;
	pgr1->copy(tmp_sg1, x2y_1, y2x_1);

	
	splice_graph tmp_sg2(gr);

	tmp_sg2.vwrt = vwrt2;
	tmp_sg2.ewrt = ewrt2;	
	MEE x2y_2;
	MEE y2x_2;
	pgr2->copy(tmp_sg2, x2y_2, y2x_2);
			
	return 0;
}


// remove edges < min_guaranteed_edge_weight, 
// remove nodes in/out-degree == 0
int phaser::refine_allelic_graphs()
{
	vector<splice_graph*> gr_pointers{pgr1, pgr2};
	for (splice_graph* pgr: gr_pointers)
	{
		for (auto ew: pgr->get_edge_weights())
		{	
			edge_descriptor e = ew.first;
			double w = ew.second;
			if (w < min_guaranteed_edge_weight) pgr->remove_edge(e);
		}

		while(true)
		{
			bool b = false;
			for(int i = 1; i < pgr->num_vertices() - 1; i++)
			{
				if(pgr->degree(i) == 0) continue;
				if(pgr->in_degree(i) >= 1 && pgr->out_degree(i) >= 1) continue;
				pgr->clear_vertex(i);
				b = true;
			}
			if(b == false) break;
		}
	}
	
	return 0;
}

int phaser::split_hs_by_rebuild()
{
	//FIXME: some AS hs should have been cut off 
	vector<splice_graph*> gr_pointers{pgr1, pgr2};
	vector<hyper_set*> hs_pointers{phs1, phs2};
	for (int i = 0; i < gr_pointers.size(); i++)
	{
		splice_graph* pgr = gr_pointers[i];
		hyper_set* phs = hs_pointers[i];
		
		//ewrt to < <int, int>, double> stewrt
		map<pair<int, int>, double> stewrt;
		for (auto&& ew: pgr->get_edge_weights())
		{	
			edge_descriptor e = ew.first;
			double w = ew.second;
			stewrt.insert({{e->source(), e->target()}, w});
		}

		// phs->add_node_list iff every edge in hs has weight > 0 in allele
		// new hs count is c or min edge weight in allele
		for (auto&& nc: sc.hs.nodes)
		{
			auto&& nodelist = nc.first;
			int c = nc.second;
			double bottleneck = c;
			bool is_removed = false;
			for(int j = 0; j < nodelist.size() - 1; j++)
			{
				pair<int, int> st{nodelist[j], nodelist[j+1]};
				auto k = stewrt.find(st);
				if(k != stewrt.end()) 
				{
					if(k->second < bottleneck) bottleneck = k->second;
				}
				else
				{
					is_removed = true;
					break;
				}
			}
			if (!is_removed && int(bottleneck) > 0)
			{
				phs->add_node_list(nodelist, int(bottleneck));
			}
		}
	}
	
	return 0;
}

/* 
*  normalize value of x and y with epsilon. x, y are ratio/counts
*  @return z = (x + eps) / (x + y + 2 * eps) 
*/
double phaser::normalize_epsilon(double x, double y)
{
	x = x / (x + y);
	y = 1 - x;
	double z = (x + epsilon) / (x + y + 2 * epsilon);
	return z;
}