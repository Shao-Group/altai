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

phaser::phaser(scallop& _sc, splice_graph* _gr1, hyper_set* _hs1, splice_graph* _gr2, hyper_set* _hs2, scallop* _sc1, scallop* _sc2)
	: sc(_sc), gr(_sc.gr), pgr1(_gr1), pgr2(_gr2), phs1(_hs1), phs2(_hs2), sc1(_sc1), sc2(_sc2)
{
	assert(sc.asnonzeroset.size() != 0); // throw runtime_error("does not have AS nodes");
	
	init();
	assign_gt();
	split_gr();
	refine_allelic_graphs();
	split_hs();
	populate_allelic_scallop(); 
}

// init ewrt1/2, countbg1/2, normalize ratiobg1/2
int phaser::init()
{
	pgr1->clear();
	pgr2->clear();
	phs1->clear();
	phs2->clear();

	strategy = "split_by_ratio";

	vwrt1.resize(gr.vwrt.size(), -1);
	vwrt2.resize(gr.vwrt.size(), -1); 

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
		// cout << "inside for loop: " << i << gt_str(gr.vinf[i].gt) << endl;
		if (gr.vinf[i].gt == ALLELE1)
		{
			// cout << "inside if1 statement" << endl;
			PEEI in = gr.in_edges(i);
			PEEI out = gr.out_edges(i);
			for (auto e = in.first; e!= in.second; e++)	
			{
				assert(gr.ewrt.find(*e) != gr.ewrt.end());
				ewrt1[*e] = gr.ewrt[*e];
				ewrt2[*e] = 0;
				ewrtbg1 += gr.ewrt[*e];
				// cout << "ale 1 in edge weight" << gr.ewrt[*e] << endl;
			}
			for (auto e = out.first; e!= out.second; e++)	
			{
				assert(gr.ewrt.find(*e) != gr.ewrt.end());
				ewrt1[*e] = gr.ewrt[*e];
				ewrt2[*e] = 0;
				ewrtbg1 += gr.ewrt[*e];
			}
			vwrtbg1 += gr.get_vertex_weight(i);
		}
		else if (gr.vinf[i].gt == ALLELE2)
		{
			// cout << "inside if2 statement" << endl;
			PEEI in = gr.in_edges(i);
			PEEI out = gr.out_edges(i);
			for (auto e = in.first; e!= in.second; e++)	
			{
				assert(gr.ewrt.find(*e) != gr.ewrt.end());
				ewrt1[*e] = 0;
				ewrt2[*e] = gr.ewrt[*e];
				ewrtbg2 += gr.ewrt[*e];
				// cout << "ale 2 in edge weight" << gr.ewrt[*e] << endl;
			}
			for (auto e = out.first; e!= out.second; e++)	
			{
				assert(gr.ewrt.find(*e) != gr.ewrt.end());
				ewrt1[*e] = 0;
				ewrt2[*e] = gr.ewrt[*e];
				ewrtbg2 += gr.ewrt[*e];
			}
			vwrtbg2 += gr.get_vertex_weight(i);
		}
	}
	pair<double, double> r1r2 = normalize_epsilon(ewrtbg1, ewrtbg2);
	ewrtratiobg1 = r1r2.first;
	ewrtratiobg2 = r1r2.second;
	assert(ewrtratiobg1 + ewrtratiobg2 < 1.001);
	assert(ewrtratiobg1 + ewrtratiobg2 > 0.999);
	assert(ewrtratiobg1 >= 0);
	assert(ewrtratiobg2 >= 0);

	if(DEBUG_MODE_ON && print_phaser_detail)
	{
		cout << "phaser ratio bg" << ewrtbg1 << "--" << ewrtbg2 << "--";
		cout << ewrtratiobg1 << "--" << ewrtratiobg2 << endl;
	}	

	//TODO: what if only one allele is expressed?
	// if (countbg1 == 0 || countbg2 == 0){;} {	// return one empty graph}	

	return 0;
}

/*
	assign edges to different gt
 	TODO: did not consider hs
*/ 
int phaser::assign_gt()
{
	// get nsnodes, sort by AS ratio
	set<int> asnodes;					// := as nodes only
	set<int> nsnodes;  					// := ns nodes only
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
	assert(asnodes.size() >= 1);
	assert(nsnodes.size() >= 1);
	
	// split local
	if (nsnodes.size() + asnodes.size() < max_num_exons)
	{
		while(nsnodes.size() >= 1)
		{
			// split nsnodes wrt descending AS ratio
			vector<int> vi = sort_nodes_by_currecnt_mae(nsnodes);
			// if(DEBUG_MODE_ON) cout << "returned sort_nodes_by_currecnt_mae" << endl;
			for(int i : vi)
			{
				// if(DEBUG_MODE_ON) cout << i << " " << endl;
				if(split_local(i))
				{
					assert(nsnodes.find(i) != nsnodes.end());
					nsnodes.erase(i);
				} 
				else break;
			}
		}
	}
	
	// split global
	while(nsnodes.size() >= 1)
	{
		// split nsnodes wrt descending AS ratio
		vector<int> vi = sort_nodes_by_currecnt_mae(nsnodes);
		for(int i : vi)
		{
			if(split_global(i))
			{
				assert(nsnodes.find(i) != nsnodes.end());
				nsnodes.erase(i);
			} 
			else break;
		}
	}
	
	return 0;
}

// returns nodes with valid mae only
vector<int> phaser::sort_nodes_by_currecnt_mae(const set<int>& s)
{
	vector< pair<double, int> > nodes_mae;
	for(int i : s)
	{
		pair<double, double> r1r2 = get_as_ratio(i);
		double mae = max(r1r2.first, r1r2.second);
		if (mae <= 0) continue;
		nodes_mae.push_back({mae, i});
	}

	sort(nodes_mae.begin(), nodes_mae.end());
	vector<int> nodes_descending_mae;

	for(auto i = nodes_mae.end(); i != nodes_mae.begin(); i--)
	{
		auto j = prev(i, 1);
		if(i!= nodes_mae.end()) 
		{
			cout << j->second << " " << j->first << " " << i->second  << " " << i->first << endl;
			assert(j->first <= i->first);
		}
		assert(j->first > 0);
		nodes_descending_mae.push_back(j->second);
	}
	
	return nodes_descending_mae;
}


/** 
 *	@param	i	node index
 *	@return		<ratio1, ratio2>, if abnormal <-1, -1>
 */ 
pair<double, double> phaser::get_as_ratio(int i)
{
	vertex_info v = gr.get_vertex_info(i);
	const PEEI in = gr.in_edges(i);
	const PEEI out = gr.out_edges(i);
	double local1 = 0;
	double local2 = 0;
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

	assert(local1 >= 0);
	assert(local2 >= 0);

	if(local1 + local2 <= 0) 
		return {-1, -1}; // w/o splitting, not being 1-degree neighbor of decomposed nodes
	else 
		return normalize_epsilon(local1, local2);
	
	assert(0);
}

bool phaser::split_local(int i)
{
	// get normalized as local ratio
	pair<double, double> r1r2 = get_as_ratio(i);
	double ratio1 = r1r2.first, ratio2 = r1r2.second;
	if(ratio1 < 0) assert(ratio2 < 0);
	if(ratio2 < 0) assert(ratio1 < 0);
	if(ratio1 + ratio2 <= 0) return false;
	
	const PEEI in = gr.in_edges(i);
	const PEEI out = gr.out_edges(i);
	
	if (strategy == "split_by_ratio")
	{
		return split_by_ratio(i, in, out, ratio1);
	}
	else
	{
		assert(0); // other strategy not implemented
		return false;
	}	
}

bool phaser::split_global(int i)
{
	// vertex_info v = gr.get_vertex_info(i);
	const PEEI in = gr.in_edges(i);
	const PEEI out = gr.out_edges(i);
	if (strategy == "split_by_ratio")
	{
		return split_by_ratio(i, in, out, ewrtratiobg1);
	}
	else
	{
		assert(0); // other strategy not implemented
		return false;
	}
}


//	split edges of vertex v, by ratio.
//  directly modify ewrt1, ewrt2, if ewrt1[e] or ewrt2[e] <= -1 (not assigned)
bool phaser::split_by_ratio(int v, const PEEI& in, const PEEI& out, double ratio_allele1)
{
	assert(ratio_allele1 > 0); // ratio normalized, won't equal. NaN handled bef calling
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
	
	return true;
}

// split edges of vertex v, by phasing path
int phaser::split_by_phasing(int v, const PEEI& in, const PEEI& out, double r1)
{
	assert(0);
	throw runtime_error("split_by_phasing not implemented yet");
	return 0;
}

int phaser::split_by_min_parsimony(int v, const PEEI& itr_in_edges, const PEEI& itr_out_edges, double ratio_allele1)
{
	assert(0);
	throw runtime_error("split_by_parsimony not defined yet");
	return -1;
}

// split sg into two pairs of sg1/hs1 and sg2/hs2
int phaser::split_gr()
{	
	x2y_1.clear();// use x2y to map original edge to new edge
	y2x_1.clear();
    x2y_2.clear();
	y2x_2.clear();

	//copy MEV(this is edge_discro), v2v
	gr.vwrt = vwrt1;
	gr.ewrt = ewrt1;
	pgr1->copy(gr, x2y_1, y2x_1);

	gr.vwrt = vwrt2;
	gr.ewrt = ewrt2;	
	pgr2->copy(gr, x2y_2, y2x_2);

	if(DEBUG_MODE_ON) 
	{
		cout << "pgr1" << pgr1->ewrt.size();
		for (auto i:pgr1->ewrt) cout << i.second << " " << endl;
	}
	if(DEBUG_MODE_ON) 
	{
		cout << "pgr2" << pgr2->ewrt.size();
		for (auto i:pgr2->ewrt) cout << i.second << " " << endl;
	}

	return 0;
}


// remove edges < min_guaranteed_edge_weight, 
// remove nodes in/out-degree == 0
int phaser::refine_allelic_graphs()
{
	vector<splice_graph*> gr_pointers{pgr1, pgr2};
	for (splice_graph* pgr: gr_pointers)
	{
		PEEI pei;
		edge_iterator it1, it2;

		// remove edges < min_guaranteed_edge_weight
		for (pei = pgr->edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{	
			edge_descriptor e = *it1;
			double w = pgr->get_edge_weight(e);
			if (w < min_guaranteed_edge_weight)
			{
				pgr->remove_edge(e);	
			} 
		}

		// recursively remove nodes in/out-degree == 0
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

	if(DEBUG_MODE_ON) 
	{
		cout << "pgr1-refine" << pgr1->ewrt.size();
		for (auto i:pgr1->ewrt) cout << i.second << " " << endl;
	}
	if(DEBUG_MODE_ON) 
	{
		cout << "pgr2-refine" << pgr2->ewrt.size();
		for (auto i:pgr2->ewrt) cout << i.second << " " << endl;
	}

	return 0;
}

/*
**	split hs0 to two allelic hs1/hs2
**	via keeping hyper_edge whose all edges' weight >= 1 in each allelic graph
**	TODO: break hyper_edge into pieces if the drop of weight is at AS pos
*/
int phaser::split_hs()
{
	for (int i = 0; i < 2; i++)
	{
		// only two potential alleles 
		assert (i == 0 || i == 1); 
		if (i == 0)
		{
			splice_graph* pgr = pgr1;
			hyper_set* phs = phs1;
			MED& ewrt_cur = ewrt1;
		}
		else if (i == 1)
		{
			splice_graph* pgr = pgr2;
			hyper_set* phs = phs2;
			MED& ewrt_cur = ewrt2;
		}
		
		// copy hs0 to hs1/hs2; remove undesired edges
		MVII edges_w_count;
		for (int j = 0; j < sc.hs.edges.size(); j++)
		{
			vector<int>& edge_idx_list = sc.hs.edges[j];
			int c = sc.hs.ecnts[j];
			double bottleneck = c;
			bool is_removed = false;
			for(int edge_idx : edge_idx_list)
			{
				edge_descriptor e = sc.i2e[edge_idx];
				try 
				{
					double w = ewrt_cur[e];
					assert(w >= 0);
					if(w < bottleneck) bottleneck = w;
				}
				catch (EdgeWeightException ex) 
				{
					is_removed = true;
					break;
				}
			}
			// add hyper_edge if all edges have AS weight > 1 (hs will be transformed)
			if (!is_removed && int(bottleneck) >= 1)
			{
				int allelic_c = int(bottleneck);
				assert(edges_w_count.find(edge_idx_list) == edges_w_count.end());
				edges_w_count.insert({edge_idx_list, allelic_c})
			}
		}
		phs->clear();
		phs->add_edge_list(edges_w_count);
	}
	return 0;
}

/*
** populate & build sc1, sc2; transform hs1, hs2;
** at the end, sc1, sc2 are ready to assemble
*/
int phaser::populate_allelic_scallop()
{
	for (int i = 0; i < 2; i++)
	{
		// only two potential alleles 
		assert (i == 0 || i == 1); 
		if (i == 0)
		{
			splice_graph* pgr = pgr1;
			hyper_set* phs = phs1;
			scallop* psc = sc1;
			// MED& ewrt_cur = ewrt1;
			MEE& x2y = x2y_1;
		}
		else if (i == 1)
		{
			splice_graph* pgr = pgr2;
			hyper_set* phs = phs2;
			scallop* psc = sc2;
			// MED& ewrt_cur = ewrt2;
			MEE& x2y = x2y_2;
		}
		//FIXME: this sc constructor will use *move* constructor of gr
		//FIXME: potential memory leak?
		*psc = scallop(splice_graph&& (*pgr), const hyper_set & (*phs), bool r = true, bool keep_as = false, const scallop &sc); 
		allelic_transform(psc, pgr, x2y);
	}
	return 0;
}

/*
** transforms edge_descriptor and other pointers from sc0/hs0 to new pointers, using x2y
** objects transformed: sc, hs
*/
int phaser::allelic_transform(scallop* psc, splice_graph* pgr, MEE& x2y)
{
	psc->transform(pgr, sc.i2e, x2y); // hs.transform called in sc
}


/* 
**  normalize value of x and y with epsilon. x, y are ratio/counts
**  @return z = (x + eps) / (x + y + 2 * eps) 
**  @return	<-1, -1> if both input are 0
*/
pair<double, double> phaser::normalize_epsilon(double x, double y)
{
	assert(x >= 0);
	assert(y >= 0);
	
	if(x+y<= 0) return {-1,-1}; // neighbors not splitted yet

	double z = (x + epsilon) / (x + y + 2 * epsilon);
	assert(z > 0 && z < 1);
	return {z, 1.0 - z};
}