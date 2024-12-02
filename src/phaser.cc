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

phaser::phaser(scallop& _sc, bool _is_allelic)
	: sc(_sc), gr(_sc.gr), is_allelic(_is_allelic)
{
	sc.gr.edge_integrity_examine();
	sc.gr.edge_integrity_enforce();

	splice_graph gr1, gr2;
	hyper_set hs1, hs2;
	pgr1 = &gr1;  
	pgr2 = &gr2;
	phs1 = &hs1;
	phs2 = &hs2;

	if (!is_allelic)
	{
		if(verbose >= 1)  printf("splice graph is not allelic, assembly of which is fully completed\n");
	}
	else if(sc.asnonzeroset.size() == 0) 
	{
		if(verbose >= 1)  printf("splice graph no longer has allelic vertices, assembly of which is resumed\n");
		assemble_scallop0(_sc);			// non-const sc0
	}
	else
	{
		init();	
		
		if (ewrtbg1 >= -0.01 && ewrtbg1 <= 0.01 && ewrtbg2 >= -0.01 && ewrtbg2 <= 0.01 && ewrtbg1 + ewrtbg2 < 0.01 && ewrtbg1 + ewrtbg2 > -0.01)
		{
			if(verbose >= 1)  printf("splice graph no longer has allelic edges, assembly of which is resumed\n");
			assemble_scallop0(_sc); 	 // non-const sc0
		}
		else
		{
			if(verbose >= 1)  
			{
				printf("partition graph %s to two allelic splice graphs, AS-vertices = %lu, overall allele frequency (%.2lf, %.2lf)\n", 
						gr.gid.c_str(), sc.asnonzeroset.size(), ewrtratiobg1, ewrtratiobg2);
			}
			assert(ewrtratiobg1 + ewrtratiobg2 < 1.001);
			assert(ewrtratiobg1 + ewrtratiobg2 > 0.999);
			assert(ewrtratiobg1 >= 0);
			assert(ewrtratiobg2 >= 0);
			
			assign_gt();
			split_gr();
			refine_allelic_graphs();
			//FIXME: revise graph
			split_hs();
			assemble_allelic_scallop(); 
			assign_transcripts_gt();
		}
	} 
}


/*
*	init ewrt1/2, countbg1/2, normalize ratiobg1/2
*
*   NOTE: assemble anyway regardless of having only one allele or two. The unexpressed allele will be filtered by filter class.
*/
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
	for(int i = 0; i < gr.num_vertices(); i++)
	{ 
		// cout << "inside for loop: " << i << gt_str(gr.vinf[i].gt) << endl;
		if (gr.get_vertex_info(i).gt == ALLELE1)
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
		else if (gr.get_vertex_info(i).gt == ALLELE2)
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
		else assert(!gr.get_vertex_info(i).is_as_vertex() || gr.get_vertex_info(i).gt == UNPHASED);
	}
	pair<double, double> r1r2 = normalize_epsilon(ewrtbg1, ewrtbg2);
	ewrtratiobg1 = r1r2.first;
	ewrtratiobg2 = r1r2.second;

	if(DEBUG_MODE_ON && print_phaser_detail)
	{
		cout << "phaser ratio bg" << ewrtbg1 << "--" << ewrtbg2 << "--";
		cout << ewrtratiobg1 << "--" << ewrtratiobg2 << endl;
	}	

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
		if(gr.vinf[i].is_as_vertex()) asnodes.insert(i);
		else nsnodes.insert(i);
	}
	for(const int i: sc.nsnonzeroset)
	{
		if(gr.vinf[i].is_as_vertex()) asnodes.insert(i);
		else nsnodes.insert(i);
	}
	assert(asnodes.size() >= 1);
	assert(nsnodes.size() >= 1);

	// split local; split nsnodes wrt descending AS ratio
	if (nsnodes.size() + asnodes.size() < max_num_exons)
	{
		while(nsnodes.size() >= 1)
		{
			vector<int> vi = sort_nodes_by_currecnt_mae(nsnodes);
			if(vi.size() == 0) break;
			bool b = false;

			for(int i : vi)
			{
				if(split_local(i))
				{
					assert(nsnodes.find(i) != nsnodes.end());
					nsnodes.erase(i);
				} 
				else 
				{
					b = true;
					break;
				}
			}
			if (b) break;
		}
	}
	
	// split global by vertex
	while(nsnodes.size() >= 1)
	{
		int i = *(nsnodes.begin());
		split_global(i);
		assert(nsnodes.find(i) != nsnodes.end());
		nsnodes.erase(i);
	}
	assert(nsnodes.size() == 0);

	// split global by edge
	// in rare cases some edges are left b/c their nodes are isolated or edges remained while incidental nodes removed
	for(const auto & ed: gr.ewrt)
	{
		split_global(ed.first);
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
		if(DEBUG_MODE_ON) if(i!= nodes_mae.end()) 
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
	assert(ewrtratiobg1 >= 0);
	assert(ewrtratiobg2 >= 0);
	assert(ewrtratiobg1 + ewrtratiobg2 > 0);

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

bool phaser::split_global(edge_descriptor e)
{
	if(strategy != "split_by_ratio") assert (0);

	if(DEBUG_MODE_ON)
	{
		assert(gr.ewrt.find(e) != gr.ewrt.end());
		assert(ewrt1.find(e) != ewrt1.end());
		assert(ewrt2.find(e) != ewrt2.end());
	}

	double ratio_allele1 = ewrtratiobg1;
	double w = gr.ewrt[e];
	assert(w >= 0);
	assert((ewrt1[e] < 0 && ewrt2[e] < 0) || (ewrt1[e] >= 0 && ewrt2[e] >= 0));

	if(ewrt1[e] < 0)
	{
		ewrt1[e] = w * ratio_allele1;
	} 
	if(ewrt2[e] < 0)
	{
		ewrt2[e] = w * (1 - ratio_allele1);
	} 
	return true;
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
		if(DEBUG_MODE_ON)
		{
			assert(gr.ewrt.find(*e) != gr.ewrt.end());
			assert(ewrt1.find(*e) != ewrt1.end());
			assert(ewrt2.find(*e) != ewrt2.end());
		}

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
		if(DEBUG_MODE_ON)
		{
			assert(gr.ewrt.find(*e) != gr.ewrt.end());
			assert(ewrt1.find(*e) != ewrt1.end());
			assert(ewrt2.find(*e) != ewrt2.end());
		}

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
	sc.gr.edge_integrity_examine();
	MED gr0_ewrt_copy;
	if(DEBUG_MODE_ON) 
	{	
		gr0_ewrt_copy = gr.ewrt;
		for (auto && ei0: gr0_ewrt_copy) assert(ei0.second >= 0);
		for (auto && ei1: ewrt1) assert(ei1.second >= 0);
		for (auto && ei2: ewrt2) assert(ei2.second >= 0);
	}

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

	if(DEBUG_MODE_ON && print_phaser_detail) 
	{
		cout << "DEBUG phaser::split_gr()" << endl;
		cout << "ewrt size:" << ewrt1.size() << endl;
		cout << "edge\tgr0.ewrt\tewrt1\tewrt2" << endl;
		assert (ewrt1.size() == gr0_ewrt_copy.size());
		assert (ewrt1.size() == ewrt2.size());

		for (int j = 0; j < ewrt1.size(); j ++)
		{
			auto i = next(ewrt1.begin(), j);
			auto k = next(ewrt2.begin(), j);
			auto l = next(gr0_ewrt_copy.begin(), j);
			assert (i->first == k->first);  // all edge_descriptors are the same before transform
			assert (l->first == i->first);
			cout << "edge " << i->first->source() << "->" << i->first->target();
			cout << "\t" << l->first << ": " << l->second;
			cout << "\t" << i->second << "\t"  << k->second << " " << endl;
		}	

		cout << "pgr1(order of ewrt may be different)\tsize: " << pgr1->ewrt.size() << "addr-" << pgr1 << endl;
		for (auto i:pgr1->ewrt) cout << "\t" << i.first << ": " << i.second << " " << endl;
		pgr1->edge_integrity_examine();

		cout << "pgr2(order of ewrt may be different)\tsize: " << pgr2->ewrt.size() << "addr-" << pgr2 << endl;
		for (auto i:pgr2->ewrt) cout << "\t" << i.first << ": " << i.second << " " << endl;
		pgr2->edge_integrity_examine();
	}
	return 0;
}


// remove edges < min_guaranteed_edge_weight, 
// remove edges incident to nodes in/out-degree == 0
int phaser::refine_allelic_graphs()
{
	vector<splice_graph*> gr_pointers{pgr1, pgr2};
	for (splice_graph* pgr: gr_pointers)
	{
		PEEI pei;
		edge_iterator it1, it2;

		// To avoid boundary error during removal, add edges into a set a prior.
		set<edge_descriptor> edges_1;
		for (pei = pgr->edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++) edges_1.insert(*it1);
		for (edge_descriptor e: edges_1)
		{
			if(e == null_edge) pgr->remove_edge(e);
			if(pgr->get_edge_weight(e) < min_guaranteed_edge_weight) pgr->remove_edge(e);
		}

		// recursively remove edges incident to nodes in/out-degree == 0
		// nodes will always remain in graph (maybe as isolated)
		while(true)
		{
			bool b = false;
			for(int i = 1; i < pgr->num_vertices() - 1; i++)
			{
				if(pgr->degree(i) == 0) continue;
				if(pgr->in_degree(i) >= 1 && pgr->out_degree(i) >= 1) continue;
				// pgr->clear_vertex(i);
				vertex_info vi = pgr->get_vertex_info(i);
				vi.type = REVIVAL_AS_UNDERSEQ;
				pgr->set_vertex_info(i, vi);
				pgr->add_edge(i, pgr->num_vertices() -1);
				b = true;
			}
			if(b == false) break;
		}
		pgr->edge_integrity_enforce();
	}
	
	if(DEBUG_MODE_ON && print_phaser_detail) 
	{
		cout << "phaser::refine_allelic_graphs done" << endl;
		pgr1->edge_integrity_examine();
		pgr2->edge_integrity_examine();
		
		cout << "pgr1-refine\tsize:" << pgr1->ewrt.size() << "\taddr-" << pgr1 << endl;
		set<edge_descriptor> gr1edges;
		for (auto i:pgr1->ewrt) 
		{
			cout << "\t" << i.first << ": " << i.second << " " << endl;
			gr1edges.insert(i.first);
		}
		if (gr1edges.size() == 0 && ewrtbg1 > 0.05) cerr << pgr1->gid << "(ale1) is empty after refining but has non-empty AS weight" << endl;

		cout << "pgr2-refine\tsize" << pgr2->ewrt.size() << "\taddr-" << pgr2 << endl;
		set<edge_descriptor> gr2edges;
		for (auto i:pgr2->ewrt) 
		{
			cout << "\t" << i.first << ": " << i.second << " " << endl;
			assert(gr1edges.find(i.first) == gr1edges.end());
			gr2edges.insert(i.first);
		}
		if (gr2edges.size() == 0 && ewrtbg2 > 0.05) cerr << pgr1->gid << "(ale2) is empty after refining but has non-empty AS weight" << endl;
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
	if(DEBUG_MODE_ON && print_phaser_detail)
	{
		cout << "sc hs edge before phaser: ";
		auto && edges = sc.hs.edges;
		for(int i = 0; i < edges.size(); i++)
		{
			printf("hyper-edge (edges) %d: ( ", i);
			printv(edges[i]);
			printf(")\n");
		}
	}

	// make vertex pos map
	map<pair<int, int> map<int, genotype> > pos2vertices_gt;
	for (int i = 0; i < gr.num_vertices(); i++)
	{
		const vertex_info& vi = gr.get_vertex_info(i);
		if(vi.gt != ALLELE1 && vi.gt != ALLELE2) continue;
		
		pair<int, int> pp {vi.lpos.p32, vi.rpos.p32};
		if(pos2vertices_gt.find(pp) == pos2vertices_gt.end())
		{
			pos2vertices_gt.insert({pp, {i, vi.gt} });
		}
		else
		{
			pos2vertices_gt.find(pp)->insert({i, vi.gt};
		}
	}

	for (int allele_index = 0; allele_index < 2; allele_index++)
	{
		// only two potential alleles 
		assert (allele_index == 0 || allele_index == 1); 
		hyper_set*    phs      = (allele_index == 0)? phs1  : phs2;
		MED&          ewrt_cur = (allele_index == 0)? ewrt1 : ewrt2;
		genotype	  gt 	   = (allele_index == 0)? ALLELE1 : ALLELE2;
		// copy hs0 to hs1/hs2; remove undesired edges
		MVII edges_w_count;
		for (int j = 0; j < sc.hs.edges.size(); j++)
		{
			const  vector<int>& edge_idx_list = sc.hs.edges[j];
			int    c                          = sc.hs.ecnts[j];
			double bottleneck                 = c;
			bool   use_this                   = true;
			for(int edge_idx : edge_idx_list)
			{
				if(edge_idx == -1) continue;

				assert(edge_idx >= 0 && edge_idx < sc.i2e.size());
				edge_descriptor e = sc.i2e[edge_idx];

				const vertex_info& vis = gr.get_vertex_info(e->source());
				int ss = -1; // alternative s or t
				bool b = false;
				if(gt_conflict(vis.gt, gt))
				{
					b = true;
					auto it1 = pos2vertices_gt.find({vis.lpos.p32, vis.rpos.p32});
					if(it1 != pos2vertices_gt.end())
					{
						for(auto it2: it1->second)
						{
							genotype gg = it2.second;
							if(gg == gt) ss = it2.first;
						}
					}
				}

				int tt = -1;
				const vertex_info& vit = gr.get_vertex_info(e->target());
				if(gt_conflict(vit.gt, gt))
				{
					b = true;
					auto it1 = pos2vertices_gt.find({vit.lpos.p32, vit.rpos.p32});
					if(it1 != pos2vertices_gt.end())
					{
						for(auto it2: it1->second)
						{
							genotype gg = it2.second;
							if(gg == gt) tt = it2.first;
						}
					}
				}
				//FIXME: assert tt to ss has only one edge
				if(b == true && (ss < 0 || tt < 0)) continue;
				if(b == true)
				{
					PED ped = gr.edge(ss, tt);
					assert(ped.second == true);
					assert(ped.first != null_edge);
					e = ped.first;
				}

				if(e == null_edge) continue;
				
				if(ewrt_cur.find(e) == ewrt_cur.end())
				{
					use_this = false;
					break;
				}
				double w = ewrt_cur[e];
				assert(w >= 0);
				if(w < bottleneck) bottleneck = w;
				
				if(int(bottleneck < 0.999)) 
				{
					use_this = false;
					break;
				}
			}
			//TODO: edit edge_idx_list
			// add hyper_edge if all edges have AS weight > 1 (hs will be transformed)
			if (use_this && int(bottleneck) >= 1)
			{
				int allelic_c = int(bottleneck);
				auto it = edges_w_count.find(edge_idx_list);
				if (it == edges_w_count.end())
				{
					edges_w_count.insert({edge_idx_list, allelic_c});
				}
				else		// it may happen if one edge is a subset of another 
				{
					it->second = (it->second > allelic_c)? it->second: allelic_c;
				}
			}
		}
		phs->clear();
		phs->add_edge_list(edges_w_count);
		if(DEBUG_MODE_ON) assert(phs->edges_to_transform.size() == edges_w_count.size());
	}

	if(DEBUG_MODE_ON && print_phaser_detail)
	{
		cout << "hs0.size=" << sc.hs.edges.size() << endl;
		for (auto phs : {phs1, phs2})
		{
			cout << "phs_" << phs << "\t";
			cout << "edges.size=" << phs->edges.size() << "\t";
			cout << "edges2tf.size" << phs->edges_to_transform.size() << endl;
			for(const auto& es: phs->edges_to_transform) {printv(es); cout << endl;}
			for(const auto& es: phs->edges_to_transform) {for(int i: es) cout << sc.i2e[i] <<" "; cout <<endl;}
		}
	}
	return 0;
}

/*
** when there is no variants phased in sc0
** 1. do not split
** 2. assemble sc0 and collect transcripts to both container
*/
int phaser::assemble_scallop0(scallop& sc)
{
	if (DEBUG_MODE_ON) sc.gr.edge_integrity_examine();

	sc.paths.clear();
	sc.trsts.clear();
	sc.non_full_trsts.clear();
	sc.assemble_continue(is_allelic);

	//FIXME: for specific transcripts, use a probabilistic model to judge whether it is true

	trsts1 = sc.trsts;
	trsts2 = sc.trsts;
	non_full_trsts1 = sc.non_full_trsts;
	non_full_trsts2 = sc.non_full_trsts;

	if (DEBUG_MODE_ON)
	{
		for(const transcript& t: trsts1) assert(t.gt == UNPHASED || t.gt == NONSPECIFIC);
		for(const transcript& t: trsts2) assert(t.gt == UNPHASED || t.gt == NONSPECIFIC);
		for(const transcript& t: non_full_trsts1) assert(t.gt == UNPHASED || t.gt == NONSPECIFIC);
		for(const transcript& t: non_full_trsts2) assert(t.gt == UNPHASED || t.gt == NONSPECIFIC);
	}
	
	if(verbose >= 2)
	{
		printf("Collected %ld transcripts from non-specific splice graph with no variants%s\n", 
				trsts1.size(), sc.gr.gid.c_str());
	}

	return 0;
}

/*
** populate & build & assemble sc1, sc2; transform hs1, hs2;
** at the end, sc1, sc2 are ready to assemble
*/
int phaser::assemble_allelic_scallop()
{
	pgr1->gid = pgr1->gid + ".allele1";
	pgr2->gid = pgr2->gid + ".allele2";
	scallop sc1(pgr1,  *phs1, sc, true, false);
	scallop sc2(pgr2,  *phs2, sc, true, false);	
	allelic_transform(sc1, pgr1, x2y_1);
	allelic_transform(sc2, pgr2, x2y_2);
	sc1.assemble(is_allelic);  
	sc2.assemble(is_allelic);  

	trsts1 = sc1.trsts;
	trsts2 = sc2.trsts;
	non_full_trsts1 = sc1.non_full_trsts;
	non_full_trsts2 = sc2.non_full_trsts;

	if(verbose >= 2)
	{
		printf("Collected  %ld transcripts from allele1, %ld transcripts from allel2, of splice graph %s\n", 
				trsts1.size(), trsts2.size(), sc.gr.gid.c_str());
	}

	return 0;
}

/*
** transforms edge_descriptor and other pointers from sc0/hs0 to new pointers, using x2y
** objects transformed: sc, hs
*/
int phaser::allelic_transform(scallop& sc1, splice_graph* pgr, MEE& x2y)
{	
	scallop* psc = &sc1;
	if(DEBUG_MODE_ON && print_phaser_detail)
	{	
		set<edge_descriptor> sc_edges;
		set<edge_descriptor> gr_edges;
		set<edge_descriptor> mev_edges;
		PEEI sc_peei = psc->gr.edges();
		PEEI gr_peei = pgr->edges();

		cout << "DEBUG phaser::allelic_transform" << endl;
		cout << "pgr addr-" << pgr << endl;
		cout << "x2y size=" << x2y.size() << " print" << endl;
		for(auto xypair: x2y)
		{
			cout << "\t" << xypair.first << "\t" << xypair.second << endl;
		}
		cout << "finished printing x2y" << endl;

		for (auto i = sc_peei.first; i != sc_peei.second; ++i) sc_edges.insert(*i);
		for (auto j = gr_peei.first; j != gr_peei.second; ++j) gr_edges.insert(*j);
		for (pair<edge_descriptor, vector<int> > ev: psc->mev) mev_edges.insert(ev.first);

		assert(sc_edges == gr_edges);
		
		// hs_edges is a subset of sc_edges
		// sc_edges is a subset of mev_edges. The latter one contains edges removed in allelic_graph_refine
		if(print_phaser_detail)
		{
			cout << "mev edges before transform" << endl;
			for(auto e: mev_edges) cout << e << endl;
			cout << "sc edges before transform" << endl;
			for(auto e: sc_edges) cout << e << endl;
		}
	}

	psc->transform(pgr, sc.i2e, x2y);  // hs.transform called in sc

	if(DEBUG_MODE_ON)
	{
		psc->gr.edge_integrity_examine();

		set<edge_descriptor> sc_edges;
		set<edge_descriptor> gr_edges;
		set<edge_descriptor> mev_edges;
		PEEI sc_peei = psc->gr.edges();
		PEEI gr_peei = pgr->edges();
		for (auto i = sc_peei.first; i != sc_peei.second; ++i) sc_edges.insert(*i);
		for (auto j = gr_peei.first; j != gr_peei.second; ++j) gr_edges.insert(*j);
		for (pair<edge_descriptor, vector<int> > ev: psc->mev) mev_edges.insert(ev.first);

		assert(sc_edges == gr_edges);
		
		// hs_edges is a subset of sc_edges
		for (auto es : psc->hs.e2s) 
		{
			edge_descriptor hs_edge = psc->i2e[es.first];
			assert(sc_edges.find(hs_edge) != sc_edges.end());
		}

		// sc_edges is a subset of mev_edges. The latter one contains edges removed in allelic_graph_refine
		if(print_phaser_detail)
		{
			cout << "mev edges after transfrom" << endl;
			for(auto e: mev_edges) cout << e << endl;
			cout << "sc edges after transfrom" << endl;
			for(auto e: sc_edges) cout << e << endl;
		}

		for (edge_descriptor e : sc_edges) assert (mev_edges.find(e) != mev_edges.end()); 
	
		cout << "DEBUG: phaser::allelic_transform is completed and all edge_descriptor are properly transformed" << endl;
	}

	return 0;
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
	
	if(x+y<= 0) return {-1,-1};  // neighbors not splitted yet

	double z = (x + epsilon) / (x + y + 2 * epsilon);
	assert(z > 0 && z < 1);
	return {z, 1.0 - z};
}

/*
*	Transcripts will be assigned gt w.r.t. allelic scallop where whey were assembled
*	Nonspecific transcripts will be picked in transcript_set/filter
*/
int phaser::assign_transcripts_gt()
{
	for(transcript& t: trsts1) t.assign_gt(ALLELE1);
	for(transcript& t: trsts2) t.assign_gt(ALLELE2);
	for(transcript& t: non_full_trsts1) t.assign_gt(ALLELE1);
	for(transcript& t: non_full_trsts2) t.assign_gt(ALLELE2);
	if(verbose >= 2)
	{
		printf("Assigned genotypes for the above collected transcripts.\n");
	}
	return 0;
}