/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cassert>
#include <cstdio>
#include <map>
#include <set>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>  

#include "bundle.h"
#include "bundle_bridge.h"
#include "region.h"
#include "config.h"
#include "util.h"
#include "undirected_graph.h"
#include "as_pos.hpp"
#include "as_pos32.hpp"
#include "interval_map.h"

using namespace std;

bundle::bundle(bundle_base &b)
	: bb(b), br(b)
{
	br.build();
	prepare();
}

bundle::~bundle()
{}


int bundle::prepare()
{
	compute_strand();
	build_intervals();
	build_partial_exons();
	build_pos_pids_map();
	pexon_jset(jset);
	return 0;
}

int bundle::build(int mode, bool revise)
{
	build_splice_graph(mode);

	if(revise && to_revise_splice_graph)  revise_splice_graph();
	else gr.refine_splice_graph();
	
	build_hyper_set();
	return 0;
}

int bundle::compute_strand()
{
	if(library_type != UNSTRANDED) assert(bb.strand != '.');
	if(library_type != UNSTRANDED) return 0;

	int n0 = 0, np = 0, nq = 0;
	for(int i = 0; i < bb.hits.size(); i++)
	{
		if(bb.hits[i].xs == '.') n0++;
		if(bb.hits[i].xs == '+') np++;
		if(bb.hits[i].xs == '-') nq++;
	}

	if(np > nq) bb.strand = '+';
	else if(np < nq) bb.strand = '-';
	else bb.strand = '.';

	return 0;
}

int bundle::build_intervals()
{
	fmap.clear();
	set<hit*> added_hit;
	for(int i = 0; i < br.fragments.size(); i++)
	{
		fragment &fr = br.fragments[i];
		if(fr.paths.size() != 1 || fr.paths[0].type != 1) continue;
		const vector<as_pos32>& vv = br.get_aligned_intervals(fr);
		if(vv.size() <= 0) continue;
		assert(vv.size() % 2 == 0);

		/*
		if (DEBUG_MODE_ON && verbose >= 10)
		{
			for(auto h: {fr.h1, fr.h2})
			{
				hit ht = *h;
				cout << ht.qname << "itv size bridged=" << ht.itv_align.size() << " =" ;
				ht.print(true);
			}
		}
		*/

		for(int k = 0; k < vv.size() / 2; k++)
		{
			int32_t p = vv[2 * k + 0];
			int32_t q = vv[2 * k + 1];
			fmap += make_pair(ROI(p, q), 1);
			// if (DEBUG_MODE_ON && verbose >= 10) cout <<"itv added" << p << "-" << q << endl;
		}
		added_hit.insert(fr.h1);
		added_hit.insert(fr.h2);
	}

	for(int i = 0; i < bb.hits.size(); i++)
	{
		hit &ht = bb.hits[i];
		if((ht.flag & 0x100) >= 1 && !use_second_alignment) continue;
		if(added_hit.find(&ht) != added_hit.end()) continue;
		// if(ht.bridged == true) continue;
		// if(br.breads.find(ht.qname) != br.breads.end()) continue;

		for(int k = 0; k < ht.itv_align.size(); k++)
		{
			int32_t s = high32(ht.itv_align[k]);
			int32_t t = low32(ht.itv_align[k]);
			fmap += make_pair(ROI(s, t), 1);
		}
		// cout << ht.qname << "unbridged itv size=" << ht.itv_align.size() << endl;
		// ht.print();
	}
	return 0;
}

int bundle::build_partial_exons()
{
	pexons.clear();
	regional.clear();

	set<int32_t> m1, m2; // junction site
	for (auto&& j: br.junctions)
	{
		m1.insert(j.lpos.p32);
		m2.insert(j.rpos.p32);
	}

	vector<region>& regions = br.regions;
	// add non-AS pexons
	for (int i = 0; i < regions.size(); i++)
	{
		region& r =  regions[i];
		if(!r.is_allelic()) 
		{
			r.rebuild(&fmap); 
			for(int k = 0; k < r.pexons.size(); k++)
			{
				partial_exon& rpe = r.pexons[k];
				partial_exon pe (rpe);
				rpe.rid = i;
				rpe.rid2 = k;
				pe.rid = i;
				pe.rid2 = k;
				pexons.push_back(pe);
			}
		}
	}

	// add AS pexons directly
	for (int i = 0; i < regions.size(); i++)
	{
		region& r =  regions[i];
		if(r.is_allelic()) 
		{
			assert(r.pexons.size() == 0);
			int ltype = -1, rtype = -1;

			// left side is not junction, not var, & (empty vertex or no pexon)  => ltype += START_BOUNDARY;
			if (m1.find(r.lpos.p32) != m1.end()) ltype = r.ltype;
			else if (i >= 1 && regions[i-1].is_allelic()) ltype = r.ltype;
			else if (i >= 1 && regions[i-1].pexons.size() == 0) ltype = START_BOUNDARY;
			else if (i >= 1 && regions[i-1].pexons[regions[i-1].pexons.size() - 1].type != EMPTY_VERTEX) ltype = r.ltype;
			else ltype = START_BOUNDARY;

			// right side is not junction, not var, & empty => rtype += END_BOUNDARY;
			if (m2.find(r.rpos.p32) != m2.end()) rtype = r.rtype;
			else if (i < regions.size() - 1 && regions[i+1].is_allelic()) rtype = r.rtype;
			else if (i < regions.size() - 1 && regions[i+1].pexons.size() == 0) rtype = END_BOUNDARY;
			else if (i < regions.size() - 1 && regions[i+1].pexons[0].type != EMPTY_VERTEX) rtype = r.rtype;
			else rtype = END_BOUNDARY;

			assert(ltype != -1);
			assert(rtype != -1);
			assert(r.ave != 0);
			// assert(r.gt != UNPHASED);  // assuming all var phased //TODO: not always true, need to handle potential seq errors

			partial_exon pe(r.lpos, r.rpos, ltype, rtype, r.gt);
			pe.assign_as_cov(r.ave, r.max, r.dev);
			pe.rid = i;
			pe.rid2 = 0;
			pe.type = 0;  // assert not EMPTY_VERTEX
			r.pexons.push_back(pe);
			assert(r.pexons.size() == 1);
			pexons.push_back(pe);
		}
	}

	// sort, make pe.pid & regional
	// pexons and regions.pexons are different, but have same pid
	sort(pexons.begin(), pexons.end());
	for (int i = 0; i < pexons.size(); i ++)
	{
		partial_exon& pe = pexons[i];
		pe.pid = i;
		if((pe.lpos != bb.lpos || pe.rpos != bb.rpos) && (pe.ltype & START_BOUNDARY) && (pe.rtype & END_BOUNDARY)) 
			regional.push_back(true);
		else regional.push_back(false);		

		// region.pexons and bundle.pexons should have the same rid, rid2, pid
		assert(pe.rid >= 0 && pe.rid < regions.size());
		assert(pe.rid2 >= 0 && pe.rid2 < regions[pe.rid].pexons.size());
		partial_exon& rpe = regions[pe.rid].pexons[pe.rid2];
		assert(pe.lpos == rpe.lpos);
		assert(pe.rpos == rpe.rpos);
		assert(rpe.pid == -1);
		assert(rpe.rid == pe.rid);
		assert(rpe.rid2 == pe.rid2);
		if(i >= 1) assert(pe.lpos.p32 >= pexons[i-1].lpos);
		rpe.pid = i;
	}

	if(DEBUG_MODE_ON)
	{
		for(const region& r: regions)
			for(const partial_exon& pe: r.pexons)
				assert(pe.pid >= 0 && pe.pid <= pexons.size());
	}

	return 0;
}

int bundle::build_pos_pids_map()
{
	pos_pids.clear();

	for(const partial_exon& pe: pexons) 
	{
		pair<int32_t, int32_t> pospair {pe.lpos.p32, pe.rpos.p32};
		auto it = pos_pids.find(pospair);
		if (it == pos_pids.end()) pos_pids.insert({pospair, {pe.pid}});
		else it->second.push_back(pe.pid);
	}

	if(DEBUG_MODE_ON)
	{
		int c = 0;
		int rpos1 = -1;
		int pid1 = -1;
		for(const auto& pp: pos_pids)
		{
			int32_t lpos2 = pp.first.first;
			int32_t rpos2 = pp.first.second;
			const vector<int>& pids = pp.second;

			assert(rpos1 <= lpos2);
			rpos1 = rpos2;
			
			int max_pid = -1;
			for (int i = 0; i < pids.size(); i++)
			{	
				assert(pids[i] > pid1);
				if (pids[i] > max_pid) max_pid = pids[i];
			}
			pid1 = max_pid;
			c += pids.size();
			assert(pid1 >= 0);
			assert(c == pid1 + 1);
		}
	}

	return 0;
}

/*
**	equivalent to `junctions` and `link_partial_exons`, but it links pid-pid, not rid-rid
** 	i.e. adjacent pid(s) are connected if threaded by reads
**	jset_to_fill := map < (in-pid, out-pid) , count>
**	computed from fragments' paths/ hits' vlist, which are indices of regions (i.e. jset in bridger)
**	convert this region index jset to pexon index jset
 */
int bundle::pexon_jset(map<pair<int, int>, int >& pexon_jset)
{	
	const vector<region>& regions = br.regions;
	pexon_jset.clear(); 						// to be filled

	// bridged fragments
	map<pair<int, int>, vector<hit*>> m;		// < <rid1, rid2>, (*hit)s >
	for(int i = 0; i < br.fragments.size(); i++)
	{
		fragment &fr = br.fragments[i];
		if(fr.paths.size() != 1 || fr.paths[0].type != 1) continue;
		const vector<int32_t>& vv = br.get_splices_region_index(fr);
		if(vv.size() <= 0) continue;

		for(int k = 0; k < vv.size() - 1; k++)
		{
			pair<int, int> xy {vv[k], vv[k+1]};
			if(m.find(xy) == m.end())
			{
				vector<hit*> hv;
				hv.push_back(fr.h1);
				m.insert({xy, hv});
			}
			else
			{
				m[xy].push_back(fr.h1);
			}
		}
	}

	// unbridged hits
	for(int i = 0; i < bb.hits.size(); i++)
	{
		if(bb.hits[i].bridged == true) continue;
		if((bb.hits[i].flag & 0x100) >= 1) continue;
		if(br.breads.find(bb.hits[i].qname) != br.breads.end()) continue;

		vector<int> v = decode_vlist(bb.hits[i].vlist);
		if(v.size() == 0) continue;

		for(int k = 0; k < v.size() - 1; k++)
		{
			pair<int, int> xy {v[k], v[k+1]};
			if(m.find(xy) == m.end())
			{
				vector<hit*> hv;
				hv.push_back(&(bb.hits[i]));
				m.insert({xy, hv});
			}
			else
			{
				m[xy].push_back(&(bb.hits[i]));
			}
		}
	}

	/*
	if (DEBUG_MODE_ON && print_bundle_detail)
	{
		for(auto it = m.begin(); it != m.end(); it++)
		{
			int rid1 = it->first.first;
			int rid2 = it->first.second;
			int c = it->second.size();
			cout << "jset m: " << rid1 << "--" << rid2 << ", counts = " << c << endl;
		}
	}
	*/

	// populate jset, between regions
	map<pair<int, int>, vector<hit*>>::iterator it;
	for(it = m.begin(); it != m.end(); it++)
	{
		vector<hit*> &v = it->second;

		int rid1 = it->first.first;
		int rid2 = it->first.second;
		assert(rid1 >= 0 && rid1 < regions.size());
		assert(rid2 >= 0 && rid2 < regions.size());
		assert(rid1 < rid2);

		int pid1 = -1;
		int pid2 = -1;
		const vector<partial_exon>& pexons1 = regions[rid1].pexons;
		const vector<partial_exon>& pexons2 = regions[rid2].pexons;
		
		// rid to pid
		// assuming an edge always connect region1's last pexon to region2's first exon
		if (pexons1.size() >= 1 && pexons2.size() >= 1)
		{				
			const partial_exon& pe1 = pexons1[pexons1.size() - 1];
			const partial_exon& pe2 = pexons2[0];
			pid1 = pe1.pid;
			pid2 = pe2.pid;
			assert(pid1 < pid2);
			if(!pexons[pid1].rpos.samepos(regions[rid1].rpos)) pid1 = -1;
			if(!pexons[pid2].lpos.samepos(regions[rid2].lpos)) pid2 = -1;
		}
		if (pid1 < 0 || pid2 < 0 )	continue;

		assert(pexon_jset.find({pid1, pid2}) == pexon_jset.end());
		pexon_jset.insert({{pid1, pid2}, v.size()});
	}

	// populate jset, within regions 
	for(it = m.begin(); it != m.end(); it++)
	{	
		vector<hit*> &v = it->second;

		int rid1 = it->first.first;
		assert(rid1 >= 0 && rid1 < regions.size());
		
		const vector<partial_exon>& pexons1 = regions[rid1].pexons;
		for (int i = 0; i < pexons1.size() - 1; i++)
		{
			const partial_exon& pe1 = pexons1[i];
			const partial_exon& pe2 = pexons1[i+1];
			int p1 = pe1.pid;
			int p2 = pe2.pid;
			assert(p1 < p2);
			if (! pe1.rpos.samepos(pe2.rpos)) continue;
			
			auto j = pexon_jset.find({p1, p2}) ;
			if (j == pexon_jset.end())
			{
				pexon_jset.insert({{p1, p2}, v.size()});
			}
			else
			{
				j->second = j->second + v.size();
			}
		}

	}

	return 0;
}

vector<int> bundle::align_hit(hit &h)
{
	bool b = true;
	vector<int> sp2;
	vector<int> v = decode_vlist(h.vlist);

	if(DEBUG_MODE_ON && print_bundle_detail)
	{
		cout << "align_hit, decode list v " ;
		for(int i: v) cout << i <<  " ";
		cout << endl;
	}

	if(v.size() == 0) return sp2;

	for(int k = 0; k < v.size(); k++)
	{
		const region& r = br.regions[v[k]];
		if(r.pexons.size() == 0) 
		{
			b = false;
			break;
		}
		for(const partial_exon& pe: r.pexons) sp2.push_back(pe.pid);
	}

	//TODO: min_flank_length filter

	vector<int> e;
	if(b == false) return e;
	else return sp2;
}

vector<int> bundle::align_fragment(fragment &fr)
{
	bool b = true;
	vector<int> sp2;
	vector<int> v = br.get_splices_region_index(fr);

	if(DEBUG_MODE_ON && print_bundle_detail)
	{
		cout << "align_fragment, decode list v " ;
		for(int i: v) cout << i <<  " ";
		cout << endl;
	}


	if(v.size() == 0) return sp2;

	for(int k = 0; k < v.size(); k++)
	{
		const region& r = br.regions[v[k]];
		if(r.pexons.size() == 0) 
		{
			b = false;
			break;
		}
		for(const partial_exon& pe: r.pexons) sp2.push_back(pe.pid);
	}

	//TODO: min_flank_length filter

	if(DEBUG_MODE_ON && print_bundle_detail)
	{
		cout << "align_fragment, aligned pid list  " ;
		for(int i: sp2) cout << i <<  " ";
		cout << endl;
	}

	vector<int> e;
	if(b == false) return e;
	else return sp2;
}

int bundle::build_splice_graph(int mode)
{
	//FIXME: only two as vertices at the same sites. no "n" no "N", must be in reference, at most two.
	
	// build graph
	gr.clear();
	if (verbose >= 3) 
		cout << "splice graph build for bundle " << bb.chrm << ":" << bb.lpos << "-" << bb.rpos << " " <<bb.strand << " strand" << endl;
	
	// vertices: start, each region, end
	gr.add_vertex();
	vertex_info vi0;
	vi0.lpos = bb.lpos;
	vi0.rpos = bb.lpos;
	vi0.as_type = START_OR_SINK;
	gr.set_vertex_weight(0, 0);
	gr.set_vertex_info(0, vi0);

	int n_as = 0;

	for(int i = 0; i < pexons.size(); i++) // vertices for each (partial) exon
	{
		const partial_exon &r = pexons[i];
		int length = r.rpos.p32 - r.lpos.p32;
		assert(length >= 1);
		gr.add_vertex();
		if(mode == 1) gr.set_vertex_weight(i + 1, r.max < min_guaranteed_edge_weight ? min_guaranteed_edge_weight : r.max);
		if(mode == 2) gr.set_vertex_weight(i + 1, r.ave < min_guaranteed_edge_weight ? min_guaranteed_edge_weight : r.ave);
		vertex_info vi;
		vi.lpos = r.lpos;
		vi.rpos = r.rpos;
		vi.length = length;
		vi.gt = r.gt;
		
		// FIXME: not complete enumeration
		// UHPHASED_MONOVAR vs NS_NONVAR
		if (gt_as(r.gt))// FIXME: keep all var
		{
			vi.as_type = AS_DIPLOIDVAR;
			n_as += 1;
		} 
		else if (r.is_allelic() && r.gt == UNPHASED)
		{
			vi.as_type = AS_DIPLOIDVAR;
			n_as += 1;
		}
		else vi.as_type = NS_NONVAR;

		vi.stddev = r.dev;// < 1.0 ? 1.0 : r.dev;
		vi.regional = regional[i];
		vi.type = pexons[i].type;
		gr.set_vertex_info(i + 1, vi);
	}

	gr.add_vertex();
	vertex_info vin;
	vin.lpos = bb.rpos;
	vin.rpos = bb.rpos;
	vin.as_type = START_OR_SINK;
	gr.set_vertex_weight(pexons.size() + 1, 0);
	gr.set_vertex_info(pexons.size() + 1, vin);

	if(verbose >= 3) cout << "splice graph build jset edges\n";

	// edges: each jset_pair => and e2w;  including adjacent pexons
	// vertics: assign as_type
	for(const auto& jset_item: jset)
	{
		int  lpid   = jset_item.first.first;
		int  rpid   = jset_item.first.second;
		int  c      = jset_item.second;

		if(lpid< 0 || rpid < 0) continue;

		edge_descriptor p = gr.add_edge(lpid + 1, rpid + 1);

		assert(c >= 1);
		edge_info ei;
		ei.weight = c;
		gr.set_edge_info(p, ei);
		gr.set_edge_weight(p, c);

		// assign as_type
		// FIXME: not complete enumeration
		// UHPHASED_MONOVAR vs NS_NONVAR
		if (!decompose_as_neighor)  // internal explore only, default false
		{
			vertex_info& vx = gr.vinf[lpid+1];
			vertex_info& vy = gr.vinf[rpid+1];
			if(vx.is_as_vertex())
			{
				if (!vy.is_as_vertex()) vy.as_type = AJ_NONVAR;
			}
			else if (vy.is_as_vertex())
			{
				if (!vx.is_as_vertex()) vx.as_type = AJ_NONVAR;
			}
		}
	}

	// edges: adjacent pexons which are not included in jset
	for(auto i = pos_pids.begin(), b = prev(pos_pids.end(), 1); i != b ; ++i)
	{
		auto j = next(i, 1);
		int32_t rpos1 = i->first.second;
		int32_t lpos2 = j->first.first;
		if (rpos1 != lpos2) continue;

		const vector<int>& pid1s = i->second;
		const vector<int>& pid2s = j->second;
		for(int lpid: pid1s)
		{
			for(int rpid: pid2s)
			{
				assert(0 <= lpid);
				assert(lpid < rpid);
				if(jset.find({lpid, rpid}) != jset.end()) continue;  // threaded by reads, already added
				
				const partial_exon& pe1 = pexons[lpid];
				const partial_exon& pe2 = pexons[rpid];
				if(gt_conflict(pe1.gt, pe2.gt)) continue;

				assert(gr.edge(lpid + 1, rpid + 1).second == false); // edge not added
				edge_descriptor p = gr.add_edge(lpid + 1, rpid + 1);
				edge_info ei;
				double c = min_guaranteed_edge_weight;
				ei.weight = c;
				gr.set_edge_info(p, ei);
				gr.set_edge_weight(p, c);
			}
		}
	}

	// edges: connecting start/end and pexons
	int ss = 0;
	int tt = pexons.size() + 1;
	for(int i = 0; i < pexons.size(); i++)
	{
		const partial_exon &r = pexons[i];

		if(r.ltype == START_BOUNDARY)
		{
			edge_descriptor p = gr.add_edge(ss, i + 1);
			// double w = r.ave;
			// if(i >= 1 && pexons[i - 1].rpos == r.lpos) w -= pexons[i - 1].ave; 
			// if(w < 1.0) w = 1.0;
			double w = min_guaranteed_edge_weight;
			if(mode == 1) w = r.max;
			if(mode == 2) w = r.ave;
			if(mode == 1 && i >= 1 && pexons[i - 1].rpos.p32 == r.lpos.p32) w -= pexons[i - 1].max;
			if(mode == 2 && i >= 1 && pexons[i - 1].rpos.p32 == r.lpos.p32) w -= pexons[i - 1].ave;
			if(w < min_guaranteed_edge_weight) w = min_guaranteed_edge_weight;

			gr.set_edge_weight(p, w);
			edge_info ei;
			ei.weight = w;
			gr.set_edge_info(p, ei);
		}

		if(r.rtype == END_BOUNDARY) 
		{
			edge_descriptor p = gr.add_edge(i + 1, tt);
			// double w = r.ave;
			// if(i < pexons.size() - 1 && pexons[i + 1].lpos == r.rpos) w -= pexons[i + 1].ave; 
			// if(w < 1.0) w = 1.0;
			double w = min_guaranteed_edge_weight;
			if(mode == 1) w = r.max;
			if(mode == 2) w = r.ave;
			if(mode == 1 && i < pexons.size() - 1 && pexons[i + 1].lpos.p32 == r.rpos.p32) w -= pexons[i + 1].max;
			if(mode == 2 && i < pexons.size() - 1 && pexons[i + 1].lpos.p32 == r.rpos.p32) w -= pexons[i + 1].ave;
			if(w < min_guaranteed_edge_weight) w = min_guaranteed_edge_weight;
			gr.set_edge_weight(p, w);
			edge_info ei;
			ei.weight = w;
			gr.set_edge_info(p, ei);
		}
	}


	gr.strand = bb.strand;
	gr.chrm = bb.chrm;

	return 0;
}

int bundle::revise_splice_graph()
{
	bool b = false;
	
	if(DEBUG_MODE_ON && print_bundle_detail && output_graphviz_files) 
		gr.graphviz(bb.chrm + "." + to_string(bb.lpos) + "." + to_string(bb.rpos) + ".bef_revise.dot");

	bool gr_not_intact = gr.refine_splice_graph();
	if (DEBUG_MODE_ON) if(gr_not_intact) gr.graphviz("DEBUG_graph_not_intact." +bb.chrm + "." + to_string(bb.lpos) + "." + to_string(bb.rpos) + ".bef_revise.dot");
	
	while(true)
	{
		b = tackle_false_boundaries();
		if(b == true) continue;

		b = remove_false_boundaries();
		if(b == true) continue;

		b = gr.remove_inner_boundaries();
		if(b == true) continue;

		b = gr.remove_small_exons();
		if(b == true) continue;

		b = gr.remove_intron_contamination();
		if(b == true) continue;

		b = gr.remove_small_junctions();
		if(b == true) gr.refine_splice_graph();
		if(b == true) continue;

		b = gr.extend_start_boundaries();
		if(b == true) continue;

		b = gr.extend_end_boundaries();
		if(b == true) continue;

		b = gr.extend_boundaries();
		if(b == true) gr.refine_splice_graph();
		if(b == true) continue;

		b = gr.keep_surviving_edges();
		if(b == true) gr.refine_splice_graph();
		if(b == true) continue;

		break;
	}

	gr.edge_integrity_examine();
	gr.refine_splice_graph();
	gr.edge_integrity_examine();

	if(DEBUG_MODE_ON && print_bundle_detail && output_graphviz_files) 
		gr.graphviz(bb.chrm + "." + to_string(bb.lpos) + "." + to_string(bb.rpos) + ".aft_revise.dot");

	return 0;
}

// add by Mingfu -- to use paired-end reads to remove false boundaries
bool bundle::remove_false_boundaries()
{
	map<int, int> fb1;		// end
	map<int, int> fb2;		// start
	for(int i = 0; i < br.fragments.size(); i++)
	{
		fragment &fr = br.fragments[i];
		if(fr.paths.size() == 1 && fr.paths[0].type == 1) continue;
		//if(fr.h1->bridged == true || fr.h2->bridged == true) continue;

		// only use uniquely aligned reads
		//if(fr.h1->nh >= 2 || fr.h2->nh >= 2) continue;
		if(br.breads.find(fr.h1->qname) != br.breads.end()) continue;

		// calculate actual length
		vector<int> v = align_fragment(fr);
		
		if(v.size() <= 1) continue;

		int32_t tlen = 0;
		int32_t offset1 = (fr.lpos - pexons[v.front()].lpos);
		int32_t offset2 = (pexons[v.back()].rpos - fr.rpos);
		for(int i = 0; i < v.size(); i++)
		{
			int32_t l = pexons[v[i]].rpos - pexons[v[i]].lpos;
			tlen += l;
		}
		tlen -= offset1;
		tlen -= offset2;

		// int u1 = gr.locate_vertex(fr.h1->rpos - 1);
		// int u2 = gr.locate_vertex(fr.h2->pos);
		const vector<int>& u1v = align_hit((*fr.h1));
		int u1 = u1v.size() > 0? u1v.back(): -1;
		const vector<int>& u2v = align_hit((*fr.h2));
		int u2 = u2v.size() > 0? u2v[0]: -1;

		if(u1 < 0 || u2 < 0) continue;
		if(u1 >= u2) continue;

		vertex_info v1 = gr.get_vertex_info(u1);
		vertex_info v2 = gr.get_vertex_info(u2);

		int types = 0;
		int32_t lengths = 0;
		for(int k = 0; k < fr.paths.size(); k++) types += fr.paths[k].type;
		for(int k = 0; k < fr.paths.size(); k++) lengths += fr.paths[k].length;

		bool use = true;
		if(fr.paths.size() == 1 && types == 2 && tlen > 10000) use = false;
		//if(fr.paths.size() == 1 && types == 2 && lengths <= 1.5 * insertsize_high) use = false;
		//if(fr.paths.size() == 1 && types == 2 && tlen <= 1.5 * insertsize_high) use = false;
		//if(fr.paths.size() == 1 && types == 2 && lengths <= 2 * tlen) use = false;

		if(verbose >= 2) printf("%s: u1 = %d, %d%s-%d%s, u2 = %d, %d%s-%d%s, h1.rpos = %d, h2.lpos = %d, #bridging = %lu, types = %d, lengths = %d, tlen = %d, use = %c\n", 
				fr.h1->qname.c_str(), u1, v1.lpos.p32, v1.lpos.ale.c_str(), v1.rpos.p32, v1.rpos.ale.c_str(), u2, v2.lpos.p32, v2.lpos.ale.c_str(), v2.rpos.p32, v2.rpos.ale.c_str(), fr.h1->rpos, fr.h2->pos, fr.paths.size(), types, lengths, tlen, use ? 'T' : 'F');

		if(use == false) continue;

		//if(gr.get_vertex_info(u1).rpos == fr.h1->rpos)
		{
			if(fb1.find(u1) != fb1.end()) fb1[u1]++;
			else fb1.insert(make_pair(u1, 1));
		}

		//if(gr.get_vertex_info(u2).lpos == fr.h2->pos)
		{
			if(fb2.find(u2) != fb2.end()) fb2[u2]++;
			else fb2.insert(make_pair(u2, 1));
		}
	}

	bool b = false;
	for(auto &x : fb1)
	{
		PEB p = gr.edge(x.first, gr.num_vertices() - 1);
		vertex_info vi = gr.get_vertex_info(x.first);
		if(vi.type == EMPTY_VERTEX) continue;
		if(p.second == false) continue;
		double w = gr.get_vertex_weight(x.first);
		double z = log(1 + w) / log(1 + x.second);
		double s = log(1 + w) - log(1 + x.second);
		if(s > 1.5) continue;
		if(verbose >= 2) printf("detect false end boundary %d with %d reads, vertex = %d, w = %.2lf, type = %d, z = %.2lf, s = %.2lf\n", vi.rpos.p32, x.second, x.first, w, vi.type, z, s); 
		//gr.remove_edge(p.first);
		vi.type = EMPTY_VERTEX;
		gr.set_vertex_info(x.first, vi);
		b = true;
	}

	for(auto &x : fb2)
	{
		PEB p = gr.edge(0, x.first);
		vertex_info vi = gr.get_vertex_info(x.first);
		if(vi.type == EMPTY_VERTEX) continue;
		if(p.second == false) continue;
		double w = gr.get_vertex_weight(x.first);
		double z = log(1 + w) / log(1 + x.second);
		double s = log(1 + w) - log(1 + x.second);
		if(s > 1.5) continue;
		if(verbose >= 2) printf("detect false start boundary %d with %d reads, vertex = %d, w = %.2lf, type = %d, z = %.2lf, s = %.2lf\n", vi.lpos.p32, x.second, x.first, w, vi.type, z, s); 
		//gr.remove_edge(p.first);
		vi.type = EMPTY_VERTEX;
		gr.set_vertex_info(x.first, vi);
		b = true;
	}
	return b;
}

bool bundle::tackle_false_boundaries()
{
	bool b = false;
	vector<int> points(pexons.size(), 0);
	for(int k = 0; k < br.fragments.size(); k++)
	{
		fragment &fr = br.fragments[k];

		if(fr.paths.size() != 1) continue;
		if(fr.paths[0].type != 2) continue;
		if(br.breads.find(fr.h1->qname) != br.breads.end()) continue;

		vector<int> v = align_fragment(fr);
		if(v.size() <= 1) continue;

		int32_t offset1 = (fr.lpos - pexons[v.front()].lpos);
		int32_t offset2 = (pexons[v.back()].rpos - fr.rpos);

		int32_t tlen = 0;
		for(int i = 0; i < v.size(); i++)
		{
			int32_t l = pexons[v[i]].rpos - pexons[v[i]].lpos;
			tlen += l;
		}
		tlen -= offset1;
		tlen -= offset2;

		// print
		//fr.print(99);
		if(verbose >= 2) printf("break fragment %s: total-length = %d, bridge-length = %d\n", fr.h1->qname.c_str(), tlen, fr.paths[0].length);
		/*
		for(int i = 0; i < v.size(); i++)
		{
			int32_t l = pexons[v[i]].rpos - pexons[v[i]].lpos;
			if(i == 0) l -= offset1;
			if(i == v.size() - 1) l -= offset2;
			printf(" vertex %d: length = %d, region = %d-%d -> %d\n", v[i], l, pexons[v[i]].lpos, pexons[v[i]].rpos, pexons[v[i]].rpos - pexons[v[i]].lpos);
		}
		*/

		if(tlen < insertsize_low / 2.0) continue;
		if(tlen > insertsize_high * 2.0) continue;
		if(tlen >= fr.paths[0].length) continue;

		for(int i = 0; i < v.size() - 1; i++)
		{
			partial_exon &px = pexons[v[i + 0]];
			partial_exon &py = pexons[v[i + 1]];
			if(px.rtype == END_BOUNDARY) 
			{
				if(verbose >= 2) printf("break ending vertex %d, pos = %d%s\n", v[i], px.rpos.p32, px.rpos.ale.c_str());
				points[v[i + 0]] += 1;
			}
			if(py.ltype == START_BOUNDARY) 
			{
				if(verbose >= 2) printf("break starting vertex %d, pos = %d%s\n", v[i + 1], py.lpos.p32, py.lpos.ale.c_str());
				points[v[i + 1]] += 1;
			}
		}
	}

	for(int k = 0; k < points.size(); k++)
	{
		if(points[k] <= 0) continue;
		vertex_info vi = gr.get_vertex_info(k + 1);
		if(vi.type == EMPTY_VERTEX) continue;
		PEB p = gr.edge(k + 1, gr.num_vertices() - 1);
		if(p.second == false) continue;
		double w = gr.get_vertex_weight(k + 1);
		double z = log(1 + w) / log(1 + points[k]);
		double s = log(1 + w) - log(1 + points[k]);
		if(verbose >= 2) printf("tackle false end boundary %d with %d reads, vertex = %d, w = %.2lf, z = %.2lf, s = %.2lf\n", pexons[k].rpos.p32, points[k], k + 1, w, z, s);
		if(s > 1.5) continue;
		vi.type = EMPTY_VERTEX;
		gr.set_vertex_info(k + 1, vi);
		b = true;
	}

	for(int k = 0; k < points.size(); k++)
	{
		if(points[k] <= 0) continue;
		vertex_info vi = gr.get_vertex_info(k + 1);
		if(vi.type == EMPTY_VERTEX) continue;
		PEB p = gr.edge(0, k + 1);
		if(p.second == false) continue;
		double w = gr.get_vertex_weight(k + 1);
		double z = log(1 + w) / log(1 + points[k]);
		double s = log(1 + w) - log(1 + points[k]);
		if(verbose >= 2) printf("tackle false start boundary %d with %d reads, vertex = %d, w = %.2lf, z = %.2lf, s = %.2lf\n", pexons[k].lpos.p32, points[k], k + 1, w, z, s);
		if(s > 1.5) continue;
		vi.type = EMPTY_VERTEX;
		gr.set_vertex_info(k + 1, vi);
		b = true;
	}

	return b;
}

int bundle::print(int index)
{
	printf("Bundle %d: ", index);

	// statistic xs
	int n0 = 0, np = 0, nq = 0;
	for(int i = 0; i < bb.hits.size(); i++)
	{
		if(bb.hits[i].xs == '.') n0++;
		if(bb.hits[i].xs == '+') np++;
		if(bb.hits[i].xs == '-') nq++;
	}

	printf("tid = %d, #hits = %lu, #partial-exons = %lu, range = %s:%d-%d, orient = %c (%d, %d, %d)\n",
			bb.tid, bb.hits.size(), pexons.size(), bb.chrm.c_str(), bb.lpos, bb.rpos, bb.strand, n0, np, nq);

	if(verbose <= 1) return 0;

	// print hits
	for(int i = 0; i < bb.hits.size(); i++) bb.hits[i].print();

	// print fmap
	/*
	for(JIMI it = fmap.begin(); it != fmap.end(); it++)
	{
		printf("bundle.fmap %d: jmap [%d%s, %d%s) -> %d\n", 
			index, lower(it->first).p32, lower(it->first).ale.c_str(), upper(it->first).p32, upper(it->first).ale.c_str(), it->second);
	}
	*/

	// print regions
	const vector<region>& regions = br.regions;
	for(int i = 0; i < regions.size(); i++)
	{
		regions[i].print(i);	
	}

	// print partial exons
	for(int i = 0; i < pexons.size(); i++)
	{
		pexons[i].print(i);
	}

	// print jset 
	for(auto i: jset)
	{
		int pid1 = i.first.first;
		int pid2 = i.first.second;
		int c = i.second;
		cout << "jset: " << pid1 << "-" << pid2 << " " << " counts = " << c << endl;
	}

	printf("\n");

	return 0;
}

int bundle::build_hyper_set()
{
	map<vector<int>, int> m;

	for(int k = 0; k < br.fragments.size(); k++)
	{
		fragment &fr = br.fragments[k];
	
		if(fr.type != 0) continue;	// note by Qimin, skip if not paired-end fragments

		if(fr.h1->paired != true) printf("error type: %d\n", fr.type);
		assert(fr.h1->paired == true);
		assert(fr.h2->paired == true);

		if(fr.paths.size() != 1) continue;
		if(fr.paths[0].type != 1) continue;

		//if(fr.h1->bridged == false) continue;
		//if(fr.h2->bridged == false) continue;

		vector<int> v = align_fragment(fr);
		
		if(m.find(v) == m.end()) m.insert(pair<vector<int>, int>(v, fr.cnt));
		else m[v] += fr.cnt;
	}
	
	// note by Qimin, bridge umi-linked fragments into one single long path  //TODO:
	for(int k = 0; k < br.umiLink.size(); k++)
	{
		vector<int> v;
		v.clear();

		int cnt = 0;

		// if only one fr, no need to bridge into longer one
		if(br.umiLink[k].size() == 1)
		{
			fragment &fr = br.fragments[(br.umiLink[k][0])];

			if(fr.paths.size() != 1) continue;

			// TODO: "bridged" may not be correct
			if(fr.h1->bridged == false) continue;
			if(fr.h2->bridged == false) continue;

			v = align_fragment(fr);
			if(fr.paths.size() != 1 || fr.paths[0].type != 1) v.clear();

			if(m.find(v) == m.end()) m.insert(pair<vector<int>, int>(v, fr.cnt));
			else m[v] += fr.cnt;

			continue;
		}

		// if multiple fr in umi-link, bridge into one single long path
		for(int kk = 0; kk < br.umiLink[k].size(); kk++)
		{
			fragment &fr = br.fragments[(br.umiLink[k][kk])];

			// if unbridge, then trucate and add to m
			if(fr.paths.size() != 1 || fr.h1->bridged == false || fr.h2->bridged == false)
			{
				if(v.size() > 0)
				{
					if(m.find(v) == m.end()) m.insert(pair<vector<int>, int>(v, cnt));
					else m[v] += cnt;
				}

				v.clear();
				cnt = 0;

				continue;
			}

			// otherwise, add and merge cur_v to v
			vector<int> cur_v = align_fragment(fr);
			if(fr.paths.size() != 1 || fr.paths[0].type != 1) cur_v.clear();

			if(cur_v.size()==0)
			{
				if(v.size() > 0)
				{
					if(m.find(v) == m.end()) m.insert(pair<vector<int>, int>(v, cnt));
					else m[v] += cnt;
				}

				v.clear();
				cnt = 0;

				continue;

			}
			cnt += fr.cnt;

			v.insert(v.end(), cur_v.begin(), cur_v.end());
			sort(v.begin(), v.end());
			vector<int>::iterator iter = unique(v.begin(),v.end());
			v.erase(iter,v.end());
		}

		if(v.size() > 0)
		{
			/*
			printf("v = ");
			for(int ii = 0; ii < v.size(); ii++)
			{
				printf("%d ", v[ii]);
			}
			printf("\n");
			*/

			if(m.find(v) == m.end()) m.insert(pair<vector<int>, int>(v, cnt));
			else m[v] += cnt;
		}
	}
	
	for(int k = 0; k < bb.hits.size(); k++)
	{
		hit &h = bb.hits[k];

		// bridged used here, but maybe okay
		if(h.bridged == true) continue;

		vector<int> v = align_hit(h);
		
		if(m.find(v) == m.end()) m.insert(pair<vector<int>, int>(v, 1));
		else m[v] += 1;
	}

	/*
	if(DEBUG_MODE_ON && print_bundle_detail) 
	{
		cout << "bundle::build_hyper_set() get path from br.fragments; stage 2; size = " << m.size() << endl;
		for(const auto & mvii:m)
		{
			const vector<int>& vi = mvii.first;
			int i = mvii.second;
			cout << "\t";
			for(int j: vi)
			{
				cout << j << ", ";
			}
			cout << ": " << i << ";" <<  endl;
		}
	}
	*/

	hs.clear();
	for(map<vector<int>, int>::iterator it = m.begin(); it != m.end(); it++)
	{
		const vector<int> &v = it->first;
		int c = it->second;
		if(v.size() >= 2) hs.add_node_list(v, c);
	}

	if(DEBUG_MODE_ON && print_bundle_detail) {cout << "build_hyper_set completed. print hs." << endl; hs.print();}

	return 0;
}
