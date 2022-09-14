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
	// FIXME: remove accordingly count of AS pos
	// remove insufficiently supported as pos
	build_junctions();
	build_regions();
	build_partial_exons();

	build_partial_exon_map();
	link_partial_exons();
	return 0;
}

int bundle::build(int mode, bool revise)
{
	build_splice_graph(mode);
	// TODO: phasing
	// TODO: future AS bridge

	if(revise == true) revise_splice_graph(); // TODO: this might result an error for allelic regions
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
	for(int i = 0; i < br.fragments.size(); i++)
	{
		fragment &fr = br.fragments[i];
		if(fr.paths.size() != 1 || fr.paths[0].type != 1) continue;
		vector<as_pos32> vv = br.get_aligned_intervals(fr);
		if(vv.size() <= 0) continue;
		assert(vv.size() % 2 == 0);

		for(int k = 0; k < vv.size() / 2; k++)
		{
			int32_t p = vv[2 * k + 0];
			int32_t q = vv[2 * k + 1];
			fmap += make_pair(ROI(p, q), 1);
		}
	}

	for(int i = 0; i < bb.hits.size(); i++)
	{
		hit &ht = bb.hits[i];
		if(ht.bridged == true) continue;
		if((ht.flag & 0x100) >= 1) continue;
		if(br.breads.find(ht.qname) != br.breads.end()) continue;
		for(int k = 0; k < ht.itvm.size(); k++)
		{
			int32_t s = high32(ht.itvm[k]);
			int32_t t = low32(ht.itvm[k]);
			fmap += make_pair(ROI(s, t), 1);
		}
	}
	return 0;
}

int bundle::build_junctions()
{
	junctions.clear();

	map<as_pos, vector<hit*>> m;		// bridged fragments
	for(int i = 0; i < br.fragments.size(); i++)
	{
		fragment &fr = br.fragments[i];
		if(fr.paths.size() != 1 || fr.paths[0].type != 1) continue;
		vector<as_pos32> vv = br.get_splices(fr);
		if(vv.size() <= 0) continue;
		assert(vv.size() % 2 == 0);

		for(int k = 0; k < vv.size() / 2; k++)
		{
			as_pos p = as_pos(pack(vv[2 * k + 0], vv[2 * k + 1]), "$" );

			if(m.find(p) == m.end())
			{
				vector<hit*> hv;
				hv.push_back(fr.h1);
				//hv.push_back(fr.h2);
				m.insert(pair< as_pos, vector<hit*> >(p, hv));
			}
			else
			{
				m[p].push_back(fr.h1);
				//m[p].push_back(fr.h2);
			}
		}
	}


	// vector<as_pos> m;
	// vector<vector<int> > n;
	// for(int i = 0; i < hits.size(); i++)
	// {
	// 	vector<as_pos> v = hits[i].spos;
	// 	if(v.size() == 0) continue;

	// 	//hits[i].print();
	// 	for(int k = 0; k < v.size(); k++)
	// 	{
	// 		as_pos p = v[k];
	// 		vector<as_pos>::iterator it = find(m.begin(), m.end(), p);
	// 		if(it == m.end())
	// 		{
	// 			vector<int> hv;
	// 			hv.push_back(i);
	// 			m.push_back(p);
	// 			n.push_back(hv);
	// 			assert(p.ale == "$");
	// 		}
	// 		else
	// 		{
	// 			n[it - m.begin()].push_back(i);
	// 		}
	// 	}
	// }
	// assert(m.size() == n.size());

	for(int i = 0; i < bb.hits.size(); i++)
	{
		if(bb.hits[i].bridged == true) continue;
		if((bb.hits[i].flag & 0x100) >= 1) continue;
		if(br.breads.find(bb.hits[i].qname) != br.breads.end()) continue;

		vector<as_pos> v = bb.hits[i].spos;
		if(v.size() == 0) continue;

		for(int k = 0; k < v.size(); k++)
		{
			as_pos p = v[k];
			if(m.find(p) == m.end())
			{
				vector<hit*> hv;
				hv.push_back(&(bb.hits[i]));
				m.insert(pair< as_pos, vector<hit*> >(p, hv));
			}
			else
			{
				m[p].push_back(&(bb.hits[i]));
			}
		}
	}
	
	// for(int it = 0; it < m.size(); it++)
	// {
	// 	vector<int> &v = n[it];
	// 	if(v.size() < min_splice_boundary_hits) continue;

	// 	int s0 = 0;
	// 	int s1 = 0;
	// 	int s2 = 0;
	// 	int nm = 0;
	// 	// check strandness of junction
	// 	for(int k = 0; k < v.size(); k++)
	// 	{
	// 		hit &h = hits[v[k]];
	// 		if(h.xs == '.') s0++;
	// 		if(h.xs == '+') s1++;
	// 		if(h.xs == '-') s2++;
	// 	}

	// 	// printf("junction: %s:%d%s-%d%s (%d, %d, %d) %d %d\n", chrm.c_str(), p1.p32, p1.ale.c_str(), p2.p32, p2.ale.c_str(), s0, s1, s2, s1 < s2 ? s1 : s2, nm);

	// 	junction jc(m[it], v.size());

	// 	if(s1 == 0 && s2 == 0) jc.strand = '.';
	// 	else if(s1 >= 1 && s2 >= 1) jc.strand = '.';
	// 	else if(s1 > s2) jc.strand = '+';
	// 	else jc.strand = '-';
	// 	junctions.push_back(jc);

	// }

	// map<int64_t, vector<hit*>>::iterator it;
	for(auto it = m.begin(); it != m.end(); it++)
	{
		vector<hit*> &v = it->second;
		if(v.size() < min_splice_boundary_hits) continue;

		as_pos32 p1 = high32(it->first);
		as_pos32 p2 = low32(it->first);

		int s0 = 0;
		int s1 = 0;
		int s2 = 0;
		int nm = 0;
		for(int k = 0; k < v.size(); k++)
		{
			nm += v[k]->nm;
			if(v[k]->xs == '.') s0++;
			if(v[k]->xs == '+') s1++;
			if(v[k]->xs == '-') s2++;
		}

		//printf("junction: %s:%d-%d (%d, %d, %d) %d\n", bb.chrm.c_str(), p1, p2, s0, s1, s2, s1 < s2 ? s1 : s2);

		junction jc(it->first, v.size());
		// jc.nm = nm;
		if(s1 == 0 && s2 == 0) jc.strand = '.';
		else if(s1 >= 1 && s2 >= 1) jc.strand = '.';
		else if(s1 > s2) jc.strand = '+';
		else jc.strand = '-';
		junctions.push_back(jc);
	}

	sort(junctions.begin(), junctions.end());

	if (verbose >= 3)
	{
		cout << "bundle build_junction DEBUG: \n junctions size = " << junctions.size() << endl; 
		for (int i = 0; i < junctions.size(); i ++) junctions[i].print(bb.chrm, i);		
	}
	return 0;
}


int bundle::build_allelic_junctions()
{
	allelic_junctions.clear();
	vector<as_pos> m;
	vector<vector<int> > n;
	map<pair<as_pos32, as_pos32>, vector<int> > consecutive_al_junction;
	for(int i = 0; i < bb.hits.size(); i++)
	{
		vector<as_pos> v = bb.hits[i].apos;
		if(v.size() == 0) continue;
		for(int k = 0; k < v.size(); k++)
		{
			as_pos p = v[k];
			vector<as_pos>::iterator it = find(m.begin(), m.end(), p);
			if(it == m.end())
			{
				vector<int> hv;
				hv.push_back(i);
				m.push_back(p);
				n.push_back(hv);
				assert(p.ale != "$");  // assert being allele
			}
			else
			{
				n[it - m.begin()].push_back(i);
			}
			// consecutive junctions
			/*
			if(k < v.size() -1)
			{
				as_pos32 p1 = low32(p);
				as_pos32 p2 = high32(v[k+1]);
				if(p1.samepos(p2)) 
				{
					pair<as_pos32, as_pos32> con_al_j = make_pair(p1, p2);
					if (consecutive_al_junction.find(con_al_j) == consecutive_al_junction.end())
					{
						vector<int> hv;
						hv.push_back(i);
						consecutive_al_junction.insert(make_pair(con_al_j, hv));
					}
					else
					{
						consecutive_al_junction.find(con_al_j)->second.push_back(i);
					}
					
				}
			}
			*/
		}
	}

	assert(m.size() == n.size());
	allelic_itv.clear();
	for (int i = 0; i < m.size(); i ++) allelic_itv.insert(make_pair(m[i], n[i]));

	// sort m, n
	vector<pair< as_pos, vector<int> > > mn (allelic_itv.begin(), allelic_itv.end());
	sort(mn.begin(), mn.end());
	m.clear();
	n.clear();
	for(int i = 0; i < mn.size(); i++) {m.push_back(mn[i].first); n.push_back(mn[i].second);}

	if (verbose >= 3) 
	{
		printf("bundle al junctions m.size() = %lu, n.size() = %lu\n", m.size(), n.size());
		for (auto&& i: m) cout << "\tapos vector " << high32(i).p32 << high32(i).ale << "-" << low32(i).p32 << low32(i).ale << endl;
		cout << "printf apos finished\n";
	}

	for(int it = 0; it < m.size(); it++)
	{
		vector<int> &v = n[it];
		if(v.size() < min_allele_overlap) continue;

		as_pos32 p1 = high32(m[it]);
		as_pos32 p2 = low32(m[it]);

		int s0 = 0;
		int s1 = 0;
		int s2 = 0;
		int nm = 0;
		// check strandness of junction
		for(int k = 0; k < v.size(); k++)
		{
			hit &h = bb.hits[v[k]];
			if(h.xs == '.') s0++;
			if(h.xs == '+') s1++;
			if(h.xs == '-') s2++;
		}

		as_pos32 p1_prev(p1.p32, "$");
		as_pos32 p2_next(p2.p32, "$");

		// Add left al junction 
		// not beginning bundle & not consecutive al site
		int itt = it - 1;
		while(itt >= 0 && m[itt].sameasitv(m[it])) --itt;  //m[itt] = previous allele
		if((it == 0 && p1.p32 > bb.lpos) || (itt <= 0 && it >= 1 && p1.p32 > bb.lpos) 
			|| (it >= 1 && itt >= 0 && low32(m[itt]).leftto(p1)))
		{
			junction jc1(p1_prev, p1, v.size());
			if(s1 == 0 && s2 == 0) jc1.strand = '.';
			else if(s1 >= 1 && s2 >= 1) jc1.strand = '.';
			else if(s1 > s2) jc1.strand = '+';
			else jc1.strand = '-';
			allelic_junctions.push_back(jc1);
		}
		
		// Add right al junction
		// Not end of bundle
		itt = it + 1;
		while(itt < m.size() && m[itt].sameasitv(m[it])) ++itt;  //m[itt] = next allele
		if (((it == m.size() - 1 && p2.p32 < bb.rpos) || (it < m.size() - 1))
			&& !(it < m.size() - 1 && itt < m.size() && high32(m[itt]).samepos(p2))) 
		{
			junction jc2(p2, p2_next, v.size());
			if(s1 == 0 && s2 == 0) jc2.strand = '.';
			else if(s1 >= 1 && s2 >= 1) jc2.strand = '.';
			else if(s1 > s2) jc2.strand = '+';
			else jc2.strand = '-';
			allelic_junctions.push_back(jc2);
		}
	}

	// add consecutive junctions
	for(auto it = consecutive_al_junction.begin(); it != consecutive_al_junction.end(); ++it)
	{
		vector<int> &v = it->second;
		if(v.size() < min_allele_overlap) continue;

		as_pos32 p1 = it->first.first;
		as_pos32 p2 = it->first.second;

		int s0 = 0;
		int s1 = 0;
		int s2 = 0;
		int nm = 0;
		// check strandness of junction
		for(int k = 0; k < v.size(); k++)
		{
			hit &h = bb.hits[v[k]];
			if(h.xs == '.') s0++;
			if(h.xs == '+') s1++;
			if(h.xs == '-') s2++;
		}

		junction jc(p1, p2, v.size());
		if(s1 == 0 && s2 == 0) jc.strand = '.';
		else if(s1 >= 1 && s2 >= 1) jc.strand = '.';
		else if(s1 > s2) jc.strand = '+';
		else jc.strand = '-';
		allelic_junctions.push_back(jc);
	}

	sort(allelic_junctions.begin(), allelic_junctions.end());
	if (verbose >= 3)
	{
		cout << "bundle build_allelic_junctions DEBUG: \n al_junctions size = " << allelic_junctions.size() << endl; 
		for (int i = 0; i < allelic_junctions.size(); i++) allelic_junctions[i].print(bb.chrm, i);
	}

	return 0;
}

int bundle::build_regions()
{
	// TODO: Fix splice types in vector as bundle_bridge::regions
	MPI s1;
	s1.insert(PI(as_pos32(bb.lpos, "$"), START_BOUNDARY));
	s1.insert(PI(as_pos32(bb.rpos, "$"), END_BOUNDARY));
	for(int i = 0; i < junctions.size(); i++)
	{
		junction &jc = junctions[i];

		double ave, dev, max;
		evaluate_rectangle(fmap, jc.lpos, jc.rpos, ave, dev, max);
		
		as_pos32 l = jc.lpos;
		as_pos32 r = jc.rpos;
		assert(l.ale == "$");
		assert(r.ale == "$");
		if(s1.find(l) == s1.end()) s1.insert(make_pair(l, LEFT_SPLICE));
		// else if(s1[l] == RIGHT_SPLICE) s1[l] = LEFT_RIGHT_SPLICE; // TODO: Fix splice types in vector as bundle_bridge::regions

		if(s1.find(r) == s1.end()) s1.insert(make_pair(r, RIGHT_SPLICE));
		// else if(s1[r] == LEFT_SPLICE) s1[r] = LEFT_RIGHT_SPLICE;
	}

	vector<PPI> v(s1.begin(), s1.end());
	sort(v.begin(), v.end());

	regions.clear();
	// add non allelic regions
	// those regions may overlap w/ allelic pexons
	
	for(int k = 0; k < v.size() - 1; k++)
	{
		as_pos32 l = v[k].first;
		as_pos32 r = v[k + 1].first;
		assert(l.ale == "$");
		assert(r.ale == "$");
		int ltype = v[k].second; 
		int rtype = v[k + 1].second;
		// if(ltype == LEFT_RIGHT_SPLICE) ltype = RIGHT_SPLICE; // TODO: Fix splice types in vector as bundle_bridge::regions
		// if(rtype == LEFT_RIGHT_SPLICE) rtype = LEFT_SPLICE;
		// if(ltype == ALLELIC_LEFT_RIGHT_SPLICE) ltype = ALLELIC_LEFT_SPLICE;
		// if(rtype == ALLELIC_LEFT_RIGHT_SPLICE) rtype = ALLELIC_RIGHT_SPLICE;
		// regions.push_back(region(l, r, ltype, rtype, &nammap, &imap));
		regions.push_back(region(l, r, ltype, rtype, &fmap, &(bb.imap))); // TODO:
	}

	if(verbose >= 3) {for(auto it = regions.begin(); it != regions.end(); ++it) {it->print(0);}}
	return 0;
}

int bundle::build_partial_exons()
{
	pexons.clear();
	regional.clear();
	// add non allelic pexons from region
	for(int i = 0; i < regions.size(); i++)
	{
		region &r = regions[i];
		for(int k = 0; k < r.pexons.size(); k++)
		{
			partial_exon &pe = r.pexons[k];
			pe.rid = i;
			pe.pid = pexons.size();
			vector<partial_exon> pev = partition_allelic_partial_exons(pe);
			for(auto& partitioned_pe: pev)
			{
				// if((partitioned_pe.lpos != lpos || partitioned_pe.rpos != rpos) && partitioned_pe.ltype == START_BOUNDARY && partitioned_pe.rtype == END_BOUNDARY) regional.push_back(true);
				if((partitioned_pe.lpos != bb.lpos || partitioned_pe.rpos != bb.rpos) && partitioned_pe.ltype == START_BOUNDARY && partitioned_pe.rtype == END_BOUNDARY) regional.push_back(true);
				else regional.push_back(false);
			}
			pexons.insert(pexons.end(), pev.begin(), pev.end());
			if(verbose >= 3) printf("regional = %s, ", regional.back() ? "TRUE" : "FALSE"); pe.print(k);
		}
	}

	// add allelic pexons directly
	for (auto it = allelic_itv.begin(); it != allelic_itv.end(); ++it)
	{
		as_pos32 l = high32(it->first);
		as_pos32 r = low32(it->first);
		int ltype, rtype;
		if (l.p32 > bb.lpos) ltype = ALLELIC_LEFT_SPLICE; 
		else if (l.p32 == bb.lpos) ltype = START_BOUNDARY;
		else assert(0 && "Allelic sites out of left boundary to bundle");

		if (r.p32 < bb.rpos) rtype = ALLELIC_RIGHT_SPLICE;
		else if (r.p32 == bb.rpos) rtype = END_BOUNDARY;
		else assert(0 && "Allelic sites out of right boundary of bundle");

		partial_exon pe(l, r, ltype, rtype);
		pe.is_allelic = true;
		assert(it->second.size() > 0);
		pe.ave = double(it->second.size()) / double(it->first.ale.length());
		pe.dev = 0.0;
		if (pe.ave < min_allele_overlap) continue;
		pexons.push_back(pe);  													//TODO: allelic regions can only be intact currently
		regional.push_back(false);
	}

	assert(regional.size() == pexons.size());
	
	// sort
	vector<pair<partial_exon, bool> > tmp;
	for(int i = 0; i < pexons.size(); ++i) tmp.push_back(make_pair(pexons[i], regional[i]));
	sort(tmp.begin(), tmp.end());
	pexons.clear();
	regional.clear();
	for(int i = 0; i < tmp.size(); ++i) {pexons.push_back(tmp[i].first); regional.push_back(tmp[i].second);}
	
	if (verbose >= 3) 
	{
		cout << "print pexons:\n";
		for (int i = 0; i< pexons.size(); i++) {cout <<"regional?" <<regional[i] << "\t"; pexons[i].print(i); }
	}
	return 0;
}

vector<partial_exon> bundle::partition_allelic_partial_exons(const partial_exon& pe)
{
	MPI s2;  // only positions of al branching sites w/o seq
	for(int i = 0; i < allelic_junctions.size(); i++)
	{
		junction &jc = allelic_junctions[i];
		as_pos32 l = jc.lpos;
		as_pos32 r = jc.rpos;
		assert(l.ale != "$" || r.ale!= "$");
		if (l.ale == "$")
		{
			if(s2.find(l) == s2.end()) s2.insert(make_pair(l, ALLELIC_LEFT_SPLICE));
			// else if(s2[l] == ALLELIC_RIGHT_SPLICE) s2[l] = ALLELIC_LEFT_RIGHT_SPLICE;  // TODO: Fix splice types in vector as bundle_bridge::regions
		}
		if(r.ale == "$")
		{
			if(s2.find(r) == s2.end()) s2.insert(make_pair(r, ALLELIC_RIGHT_SPLICE));
			// else if(s2[r] == ALLELIC_LEFT_SPLICE) s2[r] = ALLELIC_LEFT_RIGHT_SPLICE;  // TODO: Fix splice types in vector as bundle_bridge::regions
		}
	}

	vector<partial_exon> vpe;
	as_pos32 l = pe.lpos;
	as_pos32 r = pe.rpos;
	int ltype = pe.ltype;
	int rtype = pe.rtype;

	MPI::iterator al_it = s2.upper_bound(l);

	if(al_it == s2.end() || al_it->first.rightsameto(r))
	{
		partial_exon al_pe(l, r, ltype, rtype);		
		evaluate_rectangle(bb.nammap, al_pe.lpos, al_pe.rpos, al_pe.ave, al_pe.dev, al_pe.max);
		al_pe.is_allelic = false;
		if(al_pe.ave > 0) vpe.push_back(al_pe);
		return vpe;
	}

	while (al_it != s2.end() && al_it->first.leftto(r))
	{
		partial_exon al_pe(l, al_it->first, ltype, al_it->second);		
		evaluate_rectangle(bb.nammap, al_pe.lpos, al_pe.rpos, al_pe.ave, al_pe.dev, al_pe.max);
		al_pe.is_allelic = false;
		if(al_pe.ave > 0) vpe.push_back(al_pe);		// filter out non allelic pexons in allelic regions
		
		l = al_it->first;
		ltype = al_it->second;
		al_it = s2.upper_bound(l);
	}
	
	assert(al_it != s2.begin());
	al_it = --al_it;
	assert(al_it->first.leftto(r));
	if (al_it->first.leftto(r))
	{
		partial_exon al_pe(al_it->first, r, al_it->second, rtype);		
		evaluate_rectangle(bb.nammap, al_pe.lpos, al_pe.rpos, al_pe.ave, al_pe.dev, al_pe.max);
		al_pe.is_allelic = false;
		if(al_pe.ave > 0) vpe.push_back(al_pe);
	}

	return vpe;
}

int bundle::build_partial_exon_map()
{
	pmap.clear();
	pmap_na.clear();
	pmap_a.clear();
	for(int i = 0; i < pexons.size(); i++)
	{
		partial_exon &p = pexons[i];
		pair<as_pos32, as_pos32> a = make_pair(p.lpos, p.rpos);

		assert(pmap.find(a) == pmap.end());
		pmap.insert(make_pair(a, i + 1));

		if (p.lpos.ale == "$" && p.rpos.ale == "$")
		{
			assert(pmap_na.find(a) == pmap_na.end());
			pmap_na.insert(make_pair(a, i + 1));
		}
		else
		{
			assert(pmap_a.find(a) == pmap_a.end());
			pmap_a.insert(make_pair(a, i + 1));
		}
	}

	// pmap_na and pmap_a is not overlapped 
	as_pos32 a(0);
	for(auto it = pmap_na.begin(); it != pmap_na.end(); ++it)
	{
		as_pos32 l = it->first.first;
		as_pos32 r = it->first.second;
		assert(a.leftsameto(l));
		assert(l.leftto(r));
		a = r;
	}

	// print pmap
	if(verbose >=3)
	{
		cout << "Print pexons map (all, non-allelic, allelic)\n";
		printf("pexons.size()=%lu, pmap.size()=%lu, pmap_na.size()=%lu, pmap_a.size()=%lu\n", pexons.size(), pmap.size(), pmap_na.size(), pmap_a.size());
		
		cout << "Print pmap\n";
		for(auto it = pmap.begin(); it != pmap.end(); ++it)
		{
			as_pos32 l = it->first.first;
			as_pos32 r = it->first.second;
			printf("pexon %d%s-%d%s :%d\n", l.p32, l.ale.c_str(), r.p32, r.ale.c_str(), it->second);
		}
		cout << "\nPrint pmap_na\n";
		for(auto it = pmap_na.begin(); it != pmap_na.end(); ++it)
		{
			as_pos32 l = it->first.first;
			as_pos32 r = it->first.second;
			printf("pexon %d%s-%d%s :%d\n", l.p32, l.ale.c_str(), r.p32, r.ale.c_str(), it->second);
		}
		cout << "\nPrint pmap_a\n";
		for(auto it = pmap_a.begin(); it != pmap_a.end(); ++it)
		{
			as_pos32 l = it->first.first;
			as_pos32 r = it->first.second;
			printf("pexon %d%s-%d%s :%d\n", l.p32, l.ale.c_str(), r.p32, r.ale.c_str(), it->second);
		}
	}

	return 0;
}

int bundle::locate_left_partial_exon(as_pos32 x)
{
	assert(x.ale == "$");
	auto it = pmap_na.upper_bound(make_pair(x,x)); // p1>x or p1=x p2>x
	
	if (it == pmap_na.begin() ) 
	{
		if (it->first.first.p32 != x.p32) return -1;
	}
	else if (it == pmap_na.end()) it = prev(it);
	else if (it->first.first.p32 != x.p32) it = prev(it);
	
	assert(it->second >= 1);
	assert(it->second <= pexons.size());
	int k = it->second - 1;
	
	as_pos32 p1 = it->first.first;
	as_pos32 p2 = it->first.second;
	if(p2.leftsameto(x)) return -1;
	assert(p2.rightto(x));
	assert(p1.leftsameto(x));

	// if(x - p1.p32 > min_flank_length && p2.p32 - x < min_flank_length) k++; //TODO: min_flank_length in pmap
	return k;
}

int bundle::locate_right_partial_exon(as_pos32 x)
{
	assert(x.ale == "$");
	auto it = pmap_na.upper_bound(make_pair(x,x));
	if (it == pmap_na.begin()) return -1;
	it = prev(it);

	assert(it->second >= 1);
	assert(it->second <= pexons.size());
	int k = it->second - 1;
	
	as_pos32 p1 = it->first.first;
	as_pos32 p2 = it->first.second;
	if(p2.leftto(x)) return -1;
	assert(p2.rightsameto(x));
	assert(p1.leftto(x));

	// if(x - p1.p32 > min_flank_length && p2.p32 - x < min_flank_length) k++; //TODO: min_flank_length in pmap
	return k;
}

// int bundle::build_hyper_edges2()
// {
// 	if(verbose >= 3) cout<< "Build hyper edges\n";
	
// 	sort(hits.begin(), hits.end());

// 	/*
// 	printf("----------------------\n");
// 	for(int k = 9; k < hits.size(); k++) hits[k].print();
// 	printf("======================\n");
// 	*/

// 	hs.clear();

// 	string qname;
// 	int hi = -2;
// 	vector<int> sp1;
// 	for(int i = 0; i < hits.size(); i++)
// 	{
// 		hit &h = hits[i];
		
// 		/*
// 		printf("sp1 = ( ");
// 		printv(sp1);
// 		printf(")\n");
// 		h.print();
// 		*/

// 		if(h.qname != qname || h.hi != hi)
// 		{
// 			set<int> s(sp1.begin(), sp1.end());
// 			if(s.size() >= 2) hs.add_node_list(s);
// 			sp1.clear();
// 		}

// 		qname = h.qname;
// 		hi = h.hi;

// 		if((h.flag & 0x4) >= 1) continue;

// 		vector<int> sp2;
// 		for(int k = 0; k < h.itvna.size(); ++k)   // find interval match in pexon 
// 		{
// 			as_pos32 p1 = high32(h.itvna[k]);
// 			as_pos32 p2 = low32(h.itvna[k]);

// 			int k1 = locate_left_partial_exon(p1); 
// 			int k2 = locate_right_partial_exon(p2);
// 			if(k1 < 0 || k2 < 0) continue;
				
// 			for(int j = k1; j <= k2; j++) sp2.push_back(j); // b/c each region in itvna is contiguous
// 		}
// 		for(int k = 0; k < h.apos.size(); ++k)
// 		{
// 			as_pos32 p1 = high32(h.apos[k]);
// 			as_pos32 p2 = low32(h.apos[k]);
// 			auto it = pmap_a.find(make_pair(p1, p2));
			
// 			// assert(it != pmap_a.end()); //FIXME: should be right
// 			if(it == pmap_a.end()) continue; //FIXME:

// 			sp2.push_back(it->second - 1);
// 		}
// 		sort(sp2.begin(), sp2.end());
		
// 		if(sp1.size() <= 0 || sp2.size() <= 0)
// 		{
// 			sp1.insert(sp1.end(), sp2.begin(), sp2.end());
// 			continue;
// 		}

// 		int x1 = -1, x2 = -1;
// 		if(h.isize < 0) 
// 		{
// 			x1 = sp1[max_element(sp1)];
// 			x2 = sp2[min_element(sp2)];
// 		}
// 		else
// 		{
// 			x1 = sp2[max_element(sp2)];
// 			x2 = sp1[min_element(sp1)];
// 		}

// 		vector<int> sp3;
// 		bool c = bridge_read(x1, x2, sp3);

// 		//printf("=========\n");

// 		if(c == false)
// 		{
// 			set<int> s(sp1.begin(), sp1.end());
// 			if(s.size() >= 2) hs.add_node_list(s);
// 			sp1 = sp2;
// 		}
// 		else
// 		{
// 			sp1.insert(sp1.end(), sp2.begin(), sp2.end());
// 			sp1.insert(sp1.end(), sp3.begin(), sp3.end());
// 		}
// 	}

// 	return 0;
// }

// bool bundle::bridge_read(int x, int y, vector<int> &v)
// {
// 	v.clear();
// 	if(x >= y) return true;

// 	PEB e = gr.edge(x + 1, y + 1);
// 	if(e.second == true) return true;
// 	//else return false;

// 	if(y - x >= 6) return false;

// 	long max = 9999999999;
// 	vector<long> table;
// 	vector<int> trace;
// 	int n = y - x + 1;
// 	table.resize(n, 0);
// 	trace.resize(n, -1);
// 	table[0] = 1;
// 	trace[0] = -1;
// 	for(int i = x + 1; i <= y; i++)
// 	{
// 		edge_iterator it1, it2;
// 		PEEI pei;
// 		for(pei = gr.in_edges(i + 1), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
// 		{
// 			int s = (*it1)->source() - 1;
// 			int t = (*it1)->target() - 1;
// 			assert(t == i);
// 			if(s < x) continue;
// 			if(table[s - x] <= 0) continue;
// 			table[t - x] += table[s - x];
// 			trace[t - x] = s - x;
// 			if(table[t - x] >= max) return false;
// 		}
// 	}

// 	//printf("x = %d, y = %d, num-paths = %ld\n", x, y, table[n - 1]);
// 	if(table[n - 1] != 1) return false;

// 	//printf("path = ");

// 	v.clear();
// 	int p = n - 1;
// 	while(p >= 0)
// 	{
// 		p = trace[p];
// 		if(p <= 0) break;
// 		v.push_back(p + x);
// 		//printf("%d ", p + x);
// 	}
// 	//printf("\n");
// 	//assert(v.size() >= 1);

// 	return true;
// }

vector<int> bundle::align_hit(hit &h)
{
	bool b = true;
	vector<int> sp2;
	vector<as_pos> v;
	h.get_aligned_intervals(v);
	if(v.size() == 0) return sp2;

	for(int k = 0; k < v.size(); k++)
	{
		as_pos32 p1 = high32(v[k]);
		as_pos32 p2 = low32(v[k]);

		int k1 = locate_left_partial_exon(p1);
		int k2 = locate_right_partial_exon(p2);
		if(k1 < 0 || k2 < 0) b = false;
		if(b == false) break;

		for(int j = k1; j <= k2; j++) sp2.push_back(j);
	}

	vector<int> e;
	if(b == false) return e;
	else return sp2;
}

vector<int> bundle::align_fragment(fragment &fr)
{
	bool b = true;
	vector<int> sp2;
	//if(fr.paths.size() != 1 || fr.paths[0].type != 1) return sp2;
	vector<as_pos32> v = br.get_aligned_intervals(fr);
	assert(v.size() % 2 == 0);
	if(v.size() == 0) return sp2;

	for(int k = 0; k < v.size() / 2; k++)
	{
		as_pos32 p1 = v[2 * k + 0];
		as_pos32 p2 = v[2 * k + 1];
		int k1 = locate_left_partial_exon(p1);
		int k2 = locate_right_partial_exon(p2);
		if(k1 < 0 || k2 < 0) b = false;
		if(b == false) break;
		for(int j = k1; j <= k2; j++) sp2.push_back(j);
	}

	vector<int> e;
	if(b == false) return e;
	else return sp2;
}


int bundle::link_partial_exons()
{
	if(pexons.size() == 0) return 0;

	MPI lm;
	MPI rm;
	for(int i = 0; i < pexons.size(); i++)
	{
		as_pos32 l = pexons[i].lpos;
		as_pos32 r = pexons[i].rpos;

		// assert(lm.find(l) == lm.end()); FIXME: 
		// assert(rm.find(r) == rm.end());
		if(lm.find(l) != lm.end() || rm.find(r) != rm.end()) continue; // FIXME: clean
		
		lm.insert(PPI(l, i));
		rm.insert(PPI(r, i));
	}

	// junctions connect pexons
	for(int i = 0; i < junctions.size(); i++)
	{
		junction &b = junctions[i];
		MPI::iterator li;
		MPI::iterator ri;
		li = rm.find(b.lpos);
		ri = lm.find(b.rpos);
		// assert(li != rm.end());//FIXME: 
		// assert(ri != lm.end());

		if(li != rm.end() && ri != lm.end())
		{
			b.lexon = li->second;
			b.rexon = ri->second;
		}
		else //FIXME: 
		{
			b.lexon = b.rexon = -1;
		}
	}

	// allelic junctions connect allelic pexons
	for(int i = 0; i < allelic_junctions.size(); i++)
	{
		junction &b = allelic_junctions[i];
		MPI::iterator li;
		MPI::iterator ri;
		li = rm.find(b.lpos);
		ri = lm.find(b.rpos);
		// assert(li != rm.end());
		// assert(ri != lm.end()); // FIXME:

		if(li != rm.end() && ri != lm.end())
		{
			b.lexon = li->second;
			b.rexon = ri->second;
		}
		else //FIXME: 
		{
			b.lexon = -1;
			b.rexon = -1;
		}
	}

	if (verbose >= 3) 
	{
		cout << "link_exons_junction DEBUG print: real junctions:\n"; 
		for (int i = 0; i < junctions.size(); ++i) junctions[i].print(bb.chrm, i); 
		cout << "link_exons_junction DEBUG print: allelic junctions:\n";
		for (int i = 0; i < allelic_junctions.size(); ++i) allelic_junctions[i].print(bb.chrm, i); 
		cout << "link_exons_junction DEBUG print finished. \n";
	}
	
	return 0;
}

int bundle::build_splice_graph(int mode)
{
	gr.clear();
	if (verbose >= 3) 
		cout << "splice graph build for bundle " << bb.chrm << bb.lpos << "-" << bb.rpos << endl;
	// vertices: start, each region, end
	gr.add_vertex();
	vertex_info vi0;
	vi0.lpos = bb.lpos;
	vi0.rpos = bb.lpos;
	gr.set_vertex_weight(0, 0);
	gr.set_vertex_info(0, vi0);
	for(int i = 0; i < pexons.size(); i++) // vertices for each (partial) exon
	{
		const partial_exon &r = pexons[i];
		int length = r.rpos.p32 - r.lpos.p32;
		assert(length >= 1);
		gr.add_vertex();
		// gr.set_vertex_weight(i + 1, r.ave < 1.0 ? 1.0 : r.ave);
		if(mode == 1) gr.set_vertex_weight(i + 1, r.max < min_guaranteed_edge_weight ? min_guaranteed_edge_weight : r.max);
		if(mode == 2) gr.set_vertex_weight(i + 1, r.ave < min_guaranteed_edge_weight ? min_guaranteed_edge_weight : r.ave);
		vertex_info vi;
		vi.lpos = r.lpos;
		vi.rpos = r.rpos;
		vi.length = length;
		vi.stddev = r.dev;// < 1.0 ? 1.0 : r.dev;
		vi.regional = regional[i];
		vi.type = pexons[i].type;
		gr.set_vertex_info(i + 1, vi);
	}

	gr.add_vertex();
	vertex_info vin;
	vin.lpos = bb.rpos;
	vin.rpos = bb.rpos;
	gr.set_vertex_weight(pexons.size() + 1, 0);
	gr.set_vertex_info(pexons.size() + 1, vin);
	if(verbose >= 3) cout << "splice graph build junction edges\n";
	// edges: each junction => and e2w
	set<pair<int, int> > edge_set;
	for(int i = 0; i < junctions.size(); i++)
	{
		const junction &b = junctions[i];

		if(b.lexon < 0 || b.rexon < 0) continue;

		const partial_exon &x = pexons[b.lexon];
		const partial_exon &y = pexons[b.rexon];

		edge_descriptor p = gr.add_edge(b.lexon + 1, b.rexon + 1);
		edge_set.insert(make_pair(b.lexon + 1, b.rexon + 1));

		assert(b.count >= 1);
		edge_info ei;
		ei.weight = b.count;
		ei.strand = b.strand;
		gr.set_edge_info(p, ei);
		gr.set_edge_weight(p, b.count);
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

	// edges: connecting adjacent pexons => e2w
	for(int i = 0; i < (int)(pexons.size()) - 1; i++)
	{
		const partial_exon &x = pexons[i];
		int k = 1;
		for (; k < pexons.size() - 1 - i; ++k)
		{
			const partial_exon &z = pexons[i + k];
			if (x.lpos.samepos(z.lpos)) continue;
			else break;
		}
		
		while (i+k < pexons.size() && x.rpos.samepos(pexons[i+k].lpos))
		{
			// TODO: find all k, might be more than one k
			const partial_exon &y = pexons[i+k]; // y is first non-allele after x
			assert(!x.lpos.samepos(y.lpos));
			assert(x.rpos.p32 == y.lpos.p32);
			
			int xd = gr.out_degree(i + 1);
			int yd = gr.in_degree(i + k + 1);
			// double wt = (xd < yd) ? x.ave : y.ave;
			double wt = min_guaranteed_edge_weight;
			if(mode == 1) wt = (xd < yd) ? x.max: y.max;
			if(mode == 2) wt = (xd < yd) ? x.ave: y.ave;
			//int32_t xr = compute_overlap(mmap, x.rpos - 1);						
			//int32_t yl = compute_overlap(mmap, y.lpos);
			//double wt = xr < yl ? xr : yl;
			if (edge_set.find(make_pair(i + 1, i + k + 1)) == edge_set.end() ) 		// is edge present in junction
			{
				edge_descriptor p = gr.add_edge(i + 1, i + k + 1);
				// double w = (wt < 1.0) ? 1.0 : wt;
				double w = (wt < min_guaranteed_edge_weight) ? min_guaranteed_edge_weight : wt;
				gr.set_edge_weight(p, w);
				edge_info ei;
				ei.weight = w;
				gr.set_edge_info(p, ei);
			}
			k++;
		}
	}

	//assert //TODO: all vertices should have s and t, connecting to sink or tail

	gr.strand = bb.strand;
	gr.chrm = bb.chrm;

	if(DEBUG_MODE_ON) {cout << "splice graph built\n"; gr.draw("./debug.example.gr.tex");}

	return 0;
}

int bundle::revise_splice_graph()
{
	bool b = false;
	while(true)
	{
		b = tackle_false_boundaries();
		if(b == true) continue;

		// b = extend_boundaries();
		b = remove_false_boundaries();
		if(b == true) continue;

		b = remove_inner_boundaries();
		if(b == true) continue;

		b = remove_small_exons();
		// if(b == true) refine_splice_graph();
		if(b == true) continue;

		b = remove_intron_contamination();
		if(b == true) continue;

		b = remove_small_junctions();
		if(b == true) refine_splice_graph();
		if(b == true) continue;

		// b = keep_surviving_edges();
		b = extend_start_boundaries();
		if(b == true) continue;

		b = extend_end_boundaries();
		if(b == true) continue;

		b = extend_boundaries();
		if(b == true) refine_splice_graph();
		if(b == true) continue;

		// b = remove_intron_contamination();
		b = keep_surviving_edges();
		if(b == true) refine_splice_graph();
		if(b == true) continue;

		break;
	}

	refine_splice_graph();
	return 0;
}

int bundle::refine_splice_graph()
{
	while(true)
	{
		bool b = false;
		for(int i = 1; i < gr.num_vertices() - 1; i++)
		{
			if(gr.degree(i) == 0) continue;
			if(gr.in_degree(i) >= 1 && gr.out_degree(i) >= 1) continue;
			gr.clear_vertex(i);
			b = true;
		}
		if(b == false) break;
	}
	return 0;
}

bool bundle::extend_start_boundaries()
{
	bool flag = false;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		PEB p = gr.edge(0, i);
		if(p.second == true) continue;

		double wv = gr.get_vertex_weight(i);
		double we = 0;
		PEEI pei = gr.in_edges(i);
		for(edge_iterator it = pei.first; it != pei.second; it++)
		{
			we += gr.get_edge_weight(*it);
		}

		if(wv < we || wv < 10 * we * we + 10) continue;

		edge_descriptor ee = gr.add_edge(0, i);
		gr.set_edge_weight(ee, wv - we);
		gr.set_edge_info(ee, edge_info());

		vertex_info vi = gr.get_vertex_info(i);
		if(verbose >= 2) printf("extend start boundary: vertex = %d, wv = %.2lf, we = %.2lf, pos = %d%s\n", i, wv, we, vi.lpos.p32, vi.lpos.ale.c_str());

		flag = true;
	}
	return flag;
}

bool bundle::extend_end_boundaries()
{
	bool flag = false;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		PEB p = gr.edge(i, gr.num_vertices() - 1);
		if(p.second == true) continue;

		double wv = gr.get_vertex_weight(i);
		double we = 0;
		PEEI pei = gr.out_edges(i);
		for(edge_iterator it = pei.first; it != pei.second; it++)
		{
			we += gr.get_edge_weight(*it);
		}

		if(wv < we || wv < 10 * we * we + 10) continue;

		edge_descriptor ee = gr.add_edge(i, gr.num_vertices() - 1);
		gr.set_edge_weight(ee, wv - we);
		gr.set_edge_info(ee, edge_info());

		vertex_info vi = gr.get_vertex_info(i);
		if(verbose >= 2) printf("extend end boundary: vertex = %d, wv = %.2lf, we = %.2lf, pos = %d%s\n", i, wv, we, vi.rpos.p32, vi.rpos.ale.c_str());

		flag = true;
	}
	return flag;
}

bool bundle::extend_boundaries()
{
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		int s = e->source();
		int t = e->target();
		int32_t p = gr.get_vertex_info(t).lpos - gr.get_vertex_info(s).rpos;
		double we = gr.get_edge_weight(e);
		double ws = gr.get_vertex_weight(s);
		double wt = gr.get_vertex_weight(t);

		if(p <= 0) continue;
		if(s == 0) continue;
		if(t == gr.num_vertices() - 1) continue;

		bool b = false;
		if(gr.out_degree(s) == 1 && ws >= 10.0 * we * we + 10.0) b = true;
		if(gr.in_degree(t) == 1 && wt >= 10.0 * we * we + 10.0) b = true;

		if(b == false) continue;

		if(gr.out_degree(s) == 1)
		{
			edge_descriptor ee = gr.add_edge(s, gr.num_vertices() - 1);
			gr.set_edge_weight(ee, ws);
			gr.set_edge_info(ee, edge_info());
		}
		if(gr.in_degree(t) == 1)
		{
			edge_descriptor ee = gr.add_edge(0, t);
			gr.set_edge_weight(ee, wt);
			gr.set_edge_info(ee, edge_info());
		}

		gr.remove_edge(e);

		return true;
	}

	return false;
}

VE bundle::compute_maximal_edges()
{
	typedef pair<double, edge_descriptor> PDE;
	vector<PDE> ve;

	undirected_graph ug;
	edge_iterator it1, it2;
	PEEI pei;
	for(int i = 0; i < gr.num_vertices(); i++) ug.add_vertex();
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		double w = gr.get_edge_weight(e);
		int s = e->source();
		int t = e->target();
		if(s == 0) continue;
		if(t == gr.num_vertices() - 1) continue;
		ug.add_edge(s, t);
		ve.push_back(PDE(w, e));
	}

	vector<int> vv = ug.assign_connected_components();

	sort(ve.begin(), ve.end());

	for(int i = 1; i < ve.size(); i++) assert(ve[i - 1].first <= ve[i].first);

	VE x;
	set<int> sc;
	for(int i = ve.size() - 1; i >= 0; i--)
	{
		edge_descriptor e = ve[i].second;
		double w = gr.get_edge_weight(e);
		if(w < 1.5) break;
		int s = e->source();
		int t = e->target();
		if(s == 0) continue;
		if(t == gr.num_vertices()) continue;
		int c1 = vv[s];
		int c2 = vv[t];
		assert(c1 == c2);
		if(sc.find(c1) != sc.end()) continue;
		x.push_back(e);
		sc.insert(c1);
	}
	return x;
}

bool bundle::keep_surviving_edges()
{
	set<int> sv1;
	set<int> sv2;
	SE se;
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		double w = gr.get_edge_weight(*it1);
		int32_t p1 = gr.get_vertex_info(s).rpos;
		int32_t p2 = gr.get_vertex_info(t).lpos;
		if(w < min_surviving_edge_weight) continue;
		se.insert(*it1);
		sv1.insert(t);
		sv2.insert(s);
	}

	VE me = compute_maximal_edges();
	for(int i = 0; i < me.size(); i++)
	{
		edge_descriptor ee = me[i];
		se.insert(ee);
		sv1.insert(ee->target());
		sv2.insert(ee->source());
	}

	while(true)
	{
		bool b = false;
		for(SE::iterator it = se.begin(); it != se.end(); it++)
		{
			edge_descriptor e = (*it);
			int s = e->source(); 
			int t = e->target();
			if(sv1.find(s) == sv1.end() && s != 0)
			{
				edge_descriptor ee = gr.max_in_edge(s);
				assert(ee != null_edge);
				// if (ee == null_edge) throw BundleError(); 
				assert(se.find(ee) == se.end());
				se.insert(ee);
				sv1.insert(s);
				sv2.insert(ee->source());
				b = true;
			}
			if(sv2.find(t) == sv2.end() && t != gr.num_vertices() - 1)
			{
				edge_descriptor ee = gr.max_out_edge(t);
				assert(ee != null_edge);
				// if (ee == null_edge) throw BundleError(); 
				assert(se.find(ee) == se.end());
				se.insert(ee);
				sv1.insert(ee->target());
				sv2.insert(t);
				b = true;
			}
			if(b == true) break;
		}
		if(b == false) break;
	}

	VE ve;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		if(se.find(*it1) != se.end()) continue;
		ve.push_back(*it1);
	}

	for(int i = 0; i < ve.size(); i++)
	{
		if(verbose >= 2) printf("remove edge (%d, %d), weight = %.2lf\n", ve[i]->source(), ve[i]->target(), gr.get_edge_weight(ve[i]));
		gr.remove_edge(ve[i]);
	}

	if(ve.size() >= 1) return true;
	else return false;
}

bool bundle::remove_small_exons()
{
	bool flag = false;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		if(vi.type == EMPTY_VERTEX) continue;

		bool b = true;
		edge_iterator it1, it2;
		PEEI pei;
		int32_t p1 = gr.get_vertex_info(i).lpos;
		int32_t p2 = gr.get_vertex_info(i).rpos;

		if(p2 - p1 >= min_exon_length) continue;
		if(gr.degree(i) <= 0) continue;

		for(pei = gr.in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int s = e->source();
			//if(gr.out_degree(s) <= 1) b = false;
			if(s != 0 && gr.get_vertex_info(s).rpos == p1) b = false;
			if(b == false) break;
		}
		for(pei = gr.out_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int t = e->target();
			//if(gr.in_degree(t) <= 1) b = false;
			if(t != gr.num_vertices() - 1 && gr.get_vertex_info(t).lpos == p2) b = false;
			if(b == false) break;
		}

		if(b == false) continue;

		// only consider boundary small exons
		if(gr.edge(0, i).second == false && gr.edge(i, gr.num_vertices() - 1).second == false) continue;

		//gr.clear_vertex(i);
		if(verbose >= 2) printf("remove small exon: length = %d, pos = %d-%d\n", p2 - p1, p1, p2);
		vi.type = EMPTY_VERTEX;
		gr.set_vertex_info(i, vi);

		flag = true;
	}
	return flag;
}

bool bundle::remove_small_junctions()
{
	SE se;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.degree(i) <= 0) continue;

		bool b = true;
		edge_iterator it1, it2;
		PEEI pei;
		int32_t p1 = gr.get_vertex_info(i).lpos;
		int32_t p2 = gr.get_vertex_info(i).rpos;
		double wi = gr.get_vertex_weight(i);

		// compute max in-adjacent edge
		double ws = 0;
		for(pei = gr.in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int s = e->source();
			double w = gr.get_vertex_weight(s);
			if(s == 0) continue;
			if(gr.get_vertex_info(s).rpos != p1) continue;
			if(w < ws) continue;
			ws = w;
		}

		// remove small in-junction
		for(pei = gr.in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int s = e->source();
			double w = gr.get_edge_weight(e);
			if(s == 0) continue;
			if(gr.get_vertex_info(s).rpos == p1) continue;
			if(ws < 2.0 * w * w + 18.0) continue;
			if(wi < 2.0 * w * w + 18.0) continue;

			se.insert(e);
		}

		// compute max out-adjacent edge
		double wt = 0;
		for(pei = gr.out_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int t = e->target();
			double w = gr.get_vertex_weight(t);
			if(t == gr.num_vertices() - 1) continue;
			if(gr.get_vertex_info(t).lpos != p2) continue;
			if(w < wt) continue;
			wt = w;
		}

		// remove small in-junction
		for(pei = gr.out_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			double w = gr.get_edge_weight(e);
			int t = e->target();
			if(t == gr.num_vertices() - 1) continue;
			if(gr.get_vertex_info(t).lpos == p2) continue;
			if(ws < 2.0 * w * w + 18.0) continue;
			if(wi < 2.0 * w * w + 18.0) continue;

			se.insert(e);
		}

	}

	if(se.size() <= 0) return false;

	for(SE::iterator it = se.begin(); it != se.end(); it++)
	{
		edge_descriptor e = (*it);
		if(verbose >= 2) 
		{
			vertex_info v1 = gr.get_vertex_info(e->source());
			vertex_info v2 = gr.get_vertex_info(e->target());
			printf("remove small junction: length = %d, pos = %d%s-%d%s\n", v2.lpos - v1.rpos, v2.lpos.p32, v2.lpos.ale.c_str(), v1.rpos.p32, v1.rpos.ale.c_str());
		}
		gr.remove_edge(e);
	}

	return true;
}

bool bundle::remove_inner_boundaries()
{
	bool flag = false;
	int n = gr.num_vertices() - 1;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		if(vi.type == EMPTY_VERTEX) continue;
		
		if(gr.in_degree(i) != 1) continue;
		if(gr.out_degree(i) != 1) continue;

		PEEI pei = gr.in_edges(i);
		edge_iterator it1 = pei.first, it2 = pei.second;
		edge_descriptor e1 = (*it1);

		pei = gr.out_edges(i);
		it1 = pei.first;
		it2 = pei.second;
		edge_descriptor e2 = (*it1);

		int s = e1->source();
		int t = e2->target();

		if(s != 0 && t != n) continue;
		if(s != 0 && gr.out_degree(s) == 1) continue;
		if(t != n && gr.in_degree(t) == 1) continue;

		if(vi.stddev >= 0.01) continue;

		if(verbose >= 2) printf("remove inner boundary: vertex = %d, weight = %.2lf, length = %d, pos = %d-%d\n",
				i, gr.get_vertex_weight(i), vi.length, vi.lpos.p32, vi.rpos.p32);

		// gr.clear_vertex(i);
		vi.type = EMPTY_VERTEX;
		gr.set_vertex_info(i, vi);
		flag = true;
	}
	return flag;
}

bool bundle::remove_intron_contamination()
{
	bool flag = false;
	for(int i = 1; i < gr.num_vertices(); i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		if(vi.type == EMPTY_VERTEX) continue;
		
		if(gr.in_degree(i) != 1) continue;
		if(gr.out_degree(i) != 1) continue;

		edge_iterator it1, it2;
		PEEI pei = gr.in_edges(i);
		it1 = pei.first;
		edge_descriptor e1 = (*it1);
		pei = gr.out_edges(i);
		it1 = pei.first;
		edge_descriptor e2 = (*it1);
		int s = e1->source();
		int t = e2->target();
		double wv = gr.get_vertex_weight(i);

		if(s == 0) continue;
		if(t == gr.num_vertices() - 1) continue;
		if(gr.get_vertex_info(s).rpos != vi.lpos) continue;
		if(gr.get_vertex_info(t).lpos != vi.rpos) continue;

		PEB p = gr.edge(s, t);
		if(p.second == false) continue;

		edge_descriptor ee = p.first;
		double we = gr.get_edge_weight(ee);

		if(wv > we) continue;
		if(wv > max_intron_contamination_coverage) continue;

		if(verbose >= 2) printf("clear intron contamination %d, weight = %.2lf, length = %d, edge weight = %.2lf\n", i, wv, vi.length, we);

		// gr.clear_vertex(i);
		vi.type = EMPTY_VERTEX;
		gr.set_vertex_info(i, vi);

		flag = true;
	}
	return flag;
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

		int u1 = gr.locate_vertex(fr.h1->rpos - 1);
		int u2 = gr.locate_vertex(fr.h2->pos);

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
				if(verbose >= 2) printf("break starting vertex %d, pos = %d\n", v[i + 1], py.lpos.p32, py.lpos.ale.c_str());
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

// int bundle::find_contamination_chain()
// {
// 	int min_vertices = 5;
// 	double max_coverage = 4.0;
// 	int32_t max_distance = 2000;

// 	vector<int> chain;
// 	vector<string> types;
// 	for(int i = 1; i < pexons.size() - 1; i++)
// 	{
// 		string type = "";
// 		partial_exon &pe = pexons[i];
// 		if(pe.max > max_coverage) continue;

// 		if(pe.ltype == START_BOUNDARY && pe.rtype == END_BOUNDARY) type = "island";
// 		if(pe.ltype == START_BOUNDARY && pe.rtype == RIGHT_SPLICE) type = "start";
// 		if(pe.ltype == LEFT_SPLICE && pe.rtype == RIGHT_SPLICE && gr.edge(i - 1, i + 1).second == true) type = "intron";
// 		if(pe.ltype == LEFT_SPLICE && pe.rtype == END_BOUNDARY) type = "end";

// 		if(type == "") continue;

// 		chain.push_back(i);
// 		types.push_back(type);
// 	}

// 	if(chain.size() == 0) return 0;

// 	int32_t pre = 0;
// 	for(int k = 0; k < chain.size(); k++)
// 	{
// 		partial_exon &pe = pexons[chain[k]];
// 		printf("chain %d, pexon = %d, type = %s, pos = %d%s-%d%s, len = %d, cov = %.2lf, dist = %d\n", k, chain[k], types[k].c_str(), pe.lpos.p32, pe.lpos.ale.c_str(), pe.rpos.p32, pe.rpos.ale.c_str(), pe.rpos - pe.lpos, pe.max, pe.lpos - pre);
// 		pre = pe.rpos;
// 	}

// 	int k1 = -1;
// 	pre = 0 - max_distance - 1;
// 	for(int k = 0; k < chain.size(); k++)
// 	{
// 		partial_exon &pe = pexons[chain[k]];
// 		int32_t dist = pe.lpos - pre;
// 		if(dist > max_distance)
// 		{
// 			if(k - k1 > min_vertices)
// 			{
// 				printf("delete chain: %d-%d\n", k1 + 1, k - 1);
// 				for(int i = k1 + 1; i < k; i++)
// 				{
// 					gr.clear_vertex(chain[k] + 1);
// 				}
// 			}
// 			k1 = k - 1;
// 		}
// 		pre = pe.rpos;
// 	}

// 	if((int)(chain.size()) > min_vertices + k1)
// 	{
// 		printf("delete chain: %d-%d\n", k1 + 1, (int)(chain.size()) - 1);
// 		for(int i = k1 + 1; i < chain.size(); i++)
// 		{
// 			gr.clear_vertex(chain[i] + 1);
// 		}
// 	}
// 	return 0;
// }


// int bundle::count_junctions() const
// {
// 	int x = 0;
// 	for(int i = 0; i < junctions.size(); i++)
// 	{
// 		x += junctions[i].count;
// 	}
// 	return x;
// }

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

	// print regions
	for(int i = 0; i < regions.size(); i++)
	{
		regions[i].print(i);
	}

	// print junctions 
	for(int i = 0; i < junctions.size(); i++)
	{
		junctions[i].print(bb.chrm, i);
	}

	// print partial exons
	for(int i = 0; i < pexons.size(); i++)
	{
		pexons[i].print(i);
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

	
	// note by Qimin, bridge umi-linked fragments into one single long path
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

	hs.clear();
	for(map<vector<int>, int>::iterator it = m.begin(); it != m.end(); it++)
	{
		const vector<int> &v = it->first;
		int c = it->second;
		if(v.size() >= 2) hs.add_node_list(v, c);
	}
	return 0;
}
