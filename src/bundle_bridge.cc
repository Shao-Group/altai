/*
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cassert>
#include <cstdio>
#include <map>
#include <iomanip>
#include <fstream>

#include "bundle_bridge.h"
#include "region.h"
#include "config.h"
#include "util.h"
#include "undirected_graph.h"
#include "bridger.h"

bundle_bridge::bundle_bridge(bundle_base &b)
	: bb(b)
{
}

bundle_bridge::~bundle_bridge()
{}

int bundle_bridge::build()
{
	build_junctions();
	// build_allelic_junctions(); // not implemented, not used
	extend_junctions();
	build_regions();
	align_hits_transcripts();
	index_references();

	build_fragments();		// pair MEPPS

	if (verbose >= 3) print(1);

	bridger bdg(this); 		// bridge
	bdg.bridge();

	return 0;
}

int bundle_bridge::build_junctions()
{
	junctions.clear();
	map< as_pos, vector<int> > m;
	for(int i = 0; i < bb.hits.size(); i++)
	{
		vector<as_pos> v = bb.hits[i].spos;
		if(v.size() == 0) continue;

		for(int k = 0; k < v.size(); k++)
		{
			as_pos p = v[k];

			if(m.find(p) == m.end())
			{
				vector<int> hv;
				hv.push_back(i);
				m.insert(pair< as_pos, vector<int> >(p, hv));
			}
			else
			{
				m[p].push_back(i);
			}
		}
	}

	map< as_pos, vector<int> >::iterator it;
	for(it = m.begin(); it != m.end(); it++)
	{
		vector<int> &v = it->second;
		if(v.size() < min_splice_boundary_hits) continue;

		as_pos32 p1 = high32(it->first);
		as_pos32 p2 = low32(it->first);

		int s0 = 0;
		int s1 = 0;
		int s2 = 0;
		int nm = 0;
		for(int k = 0; k < v.size(); k++)
		{
			hit &h = bb.hits[v[k]];
			nm += h.nm;
			if(h.xs == '.') s0++;
			if(h.xs == '+') s1++;
			if(h.xs == '-') s2++;
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
		cout << "bundle_bridge build_junction: \n junctions size = " << junctions.size() << endl; 
		for (int i = 0; i < junctions.size(); i ++) junctions[i].print("NA", i);		
	}

	return 0;
}

int bundle_bridge::build_allelic_junctions()
{
	/*
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
		cout << "bundle_bridge build_allelic_junctions DEBUG: \n al_junctions size = " << allelic_junctions.size() << endl; 
		for (int i = 0; i < allelic_junctions.size(); i++) allelic_junctions[i].print(bb.chrm, i);
	}
	*/
	throw runtime_error("bundle_bridge::build_allelic_junctions() not implemented as it won't be used.");
	return 0;	
}

int bundle_bridge::extend_junctions()  // not used w/o ref
{
	map< as_pos, vector<int> > m;
	for(int i = 0; i < ref_trsts.size(); i++)
	{
		vector<PI32> v = ref_trsts[i].get_intron_chain();
		for(int k = 0; k < v.size(); k++)
		{
			assert(v[k].first < v[k].second);
			// TODO, TODO
			if(v[k].first <= bb.lpos) continue;
			if(v[k].second >= bb.rpos) continue;
			as_pos p = as_pos(pack(v[k].first.p32, v[k].second.p32), v[k].first.ale);

			if(m.find(p) == m.end())
			{
				vector<int> hv;
				hv.push_back(i);
				m.insert(pair< as_pos, vector<int> >(p, hv));
			}
			else
			{
				m[p].push_back(i);
			}
		}
	}

	map< as_pos, vector<int> >::iterator it;
	for(it = m.begin(); it != m.end(); it++)
	{
		vector<int> &v = it->second;

		int s0 = 0;
		int s1 = 0;
		int s2 = 0;
		for(int k = 0; k < v.size(); k++)
		{
			char c = ref_trsts[v[k]].strand;
			if(c == '.') s0++;
			if(c == '+') s1++;
			if(c == '-') s2++;
		}

		junction jc(it->first, 0 - v.size());

		if(s1 == 0 && s2 == 0) jc.strand = '.';
		else if(s1 >= 1 && s2 >= 1) jc.strand = '.';
		else if(s1 > s2) jc.strand = '+';
		else jc.strand = '-';
		junctions.push_back(jc);
	}
	return 0;
}

int bundle_bridge::build_regions()
{
	/*
		integral positions of splice-sites are identified based on h.spos.
		....................  pseudo-splice sites ...............  h.apos.
		Splice types of such positions are stored in `pos_splicetypes`.
		ALLELIC_LEFT_SPLICE is starting position of a variant (inclusive)
		ALLELIC_RIGHT_SPLICE is ending position of a variant (exclusive)
		set<splice_types> is converted to int constant by `splicetype_set_to_int`
	*/
	
	map<int, set<int> > pos_splicetypes;  // < position_int, {splice_types} >

	// add non-allelic pos to pos_splicetypes
	pos_splicetypes.insert({bb.lpos, {START_BOUNDARY}});
	pos_splicetypes.insert({bb.rpos, {END_BOUNDARY}});
	for(int i = 0; i < junctions.size(); i++)
	{
		junction &jc = junctions[i];
		int32_t l = jc.lpos.p32;
		int32_t r = jc.rpos.p32;
		pos_splicetypes[l].insert(LEFT_SPLICE);
		pos_splicetypes[r].insert(RIGHT_SPLICE);
	}
	
	map<pair<int, int>, set<string> > poses_seqs;  // < (pos, pos), {allele_seqs} >
	for (const hit& h: bb.hits)
	{
		for (const as_pos& p: h.apos)
		{
			pair<int, int> p_int {high32(p), low32(p)};
			poses_seqs[p_int].insert(p.ale);
		}
	}

	// add allelic pos (int) to pos_splicetypes
	for (auto&& a: poses_seqs)
	{
		auto poses = a.first;
		int32_t l = poses.first;
		int32_t r = poses.second;
		pos_splicetypes[l].insert(ALLELIC_LEFT_SPLICE);
		pos_splicetypes[r].insert(ALLELIC_RIGHT_SPLICE);
	}

	if (verbose >= 3)
	{
		cout << "bundle_bridge build regions" << endl;
		for (auto && p: pos_splicetypes)
		{
			cout << "pos_splicetypes " << p.first << ": {" ;
			for (auto && ii : p.second) cout << ", " << ii;
			cout << "}" << endl;
		}
		for(auto && p: poses_seqs)
		{
			cout << "poses_seqs (" << p.first.first << ", " << p.first.second << "): {";
			for (auto && ii: p.second) cout << ii << ", ";
			cout << "}" << endl;
		}
	}

	regions.clear();
	auto i1 = pos_splicetypes.begin();
	auto i2 = poses_seqs.begin();
	while (i1 != pos_splicetypes.end() && i2 != poses_seqs.end())
	{
		int32_t l1 = i1->first;
		int32_t r1 = std::next(i1, 1)->first;
		int32_t l2 = i2->first.first;
		int32_t r2 = i2->first.second;
		auto&& ltypes = i1->second;
		auto&& rtypes = std::next(i1)->second;

		assert(l2 >= l1);

		as_pos32 l, r;
		int ltype, rtype;
		if (l2 >= r1) // non-AS region
		{
			l = l1;
			r = r1;
			ltype = splicetype_set_to_int(ltypes);
			rtype = splicetype_set_to_int(rtypes);
			i1 ++;
			region rr(l, r, ltype, rtype);
			evaluate_rectangle(bb.mmap, l, r, rr.ave, rr.dev, rr.max);
			regions.push_back(rr);
		}
		else  // AS region, build all variants at same position
		{
			assert (l1 == l2);
			assert (r1 == r2);
			for (auto&& a: i2->second)
			{
				l = as_pos32(l2, a);
				r = as_pos32(r2, a);
				ltype = ALLELIC_LEFT_SPLICE; 
				rtype = ALLELIC_RIGHT_SPLICE;
				region rr(l, r, ltype, rtype);
				// evaluate_rectangle(bb.mmap, l, r, rr.ave, rr.dev, rr.max);
				// FIXME:  
				// rr.assign_as_cov(); 
				regions.push_back(rr);
			}
			i2 ++;
			i1 ++;
		}
	}
	assert (i2 == poses_seqs.end());
	while (i1 != prev(pos_splicetypes.end()) )  // remaining non-AS region 
	{
		int32_t l1 = i1->first;
		int32_t r1 = std::next(i1)->first;
		auto&& ltypes = i1->second;
		auto&& rtypes = std::next(i1)->second;

		as_pos32 l, r;
		int ltype, rtype;
		
		l = l1;
		r = r1;
		ltype = splicetype_set_to_int(ltypes);
		rtype = splicetype_set_to_int(rtypes);
		i1 ++;
		region rr(l, r, ltype, rtype);
		evaluate_rectangle(bb.mmap, l, r, rr.ave, rr.dev, rr.max);
		regions.push_back(rr);
		
	}
	sort(regions.begin(), regions.end());

	// print & assert
	if (verbose >= 3)
	{
		for (auto&& r: regions) r.print(123);
	}
	if (DEBUG_MODE_ON)
	{
		for(int k = 0; k < regions.size(); k++)
		{
			if(k >= 1) 
			{
				bool _continuous = regions[k - 1].rpos.samepos(regions[k].lpos);
				bool _same = regions[k - 1].lpos.samepos(regions[k].lpos) && regions[k - 1].rpos.samepos(regions[k].rpos);
				assert( _continuous || _same );
			}
			
		}
	}

	return 0;
}


int bundle_bridge::splicetype_set_to_int(set<int>& s)
{
	int ss = 0;
	for (int i: s)
	{
		ss += i;
	}
	return ss;
}


int bundle_bridge::align_hits_transcripts()
{
	map<as_pos32, int> m1;
	map<as_pos32, int> m2;
	for(int k = 0; k < regions.size(); k++)
	{
		m1.insert(pair<as_pos32, int>(regions[k].lpos, k));	
		m2.insert(pair<as_pos32, int>(regions[k].rpos, k));	
	}	
	
	if (DEBUG_MODE_ON)
	{
		for (auto&& i: m1)
		{
			as_pos32 pp = i.first;
			as_pos32 cc = i.second;
			cout << "bundle_bridge::align_hits_transcripts() m1::(region.lpos, idx) = " << pp.aspos32string() << " " << cc << endl;
		}
		for (auto&& i: m2)
		{
			as_pos32 pp = i.first;
			as_pos32 cc = i.second;
			cout << "bundle_bridge::align_hits_transcripts() m2::(region.rpos, idx) = " << pp.aspos32string() << " " << cc << endl;
		}
	}

	for(int i = 0; i < bb.hits.size(); i++)
	{
		align_hit(m1, m2, bb.hits[i], bb.hits[i].vlist);
		bb.hits[i].vlist = encode_vlist(bb.hits[i].vlist);
	}

	ref_phase.resize(ref_trsts.size());
	for(int i = 0; i < ref_trsts.size(); i++)
	{
		align_transcript(m1, ref_trsts[i], ref_phase[i]);
	}

	if (verbose >= 3)
	{
		for (auto&& h : bb.hits)
		{
			h.print();
			cout << h.qname << " vlist: [";
			printv(decode_vlist(h.vlist));
			cout << "]" << endl;
		}
		for (auto && r: regions) r.print(123);
	}

	return 0;
}

int bundle_bridge::align_hit(const map<as_pos32, int> &m1, const map<as_pos32, int> &m2, const hit &h, vector<int> &vv)
{
	vv.clear();
	vector<as_pos> v;
	h.get_aligned_intervals(v);
	if(v.size() == 0 && !h.has_variant() ) return 0;

	vector<PI> sp;
	sp.resize(v.size());

	as_pos32 p1 = high32(v.front());
	as_pos32 p2 = low32(v.back());

	sp[0].first = locate_region_left(m1, p1);
	for(int k = 1; k < v.size(); k++)
	{
		p1 = high32(v[k]);
		auto it = m1.find(p1);
		cout << p1.aspos32string() << "124" << endl;
		assert(it != m1.end());
		sp[k].first = it->second;
	}

	sp[sp.size() - 1].second = locate_region_right(m2, p2);
	for(int k = 0; k < v.size() - 1; k++)
	{
		p2 = low32(v[k]);
		auto it = m2.find(p2);
		assert(it != m2.end());
		sp[k].second = it->second; 
	}

	for(int k = 0; k < sp.size(); k++)
	{
		assert(sp[k].first <= sp[k].second);
		if(k > 0) assert(sp[k - 1].second < sp[k].first);
		for(int j = sp[k].first; j <= sp[k].second; j++) 
		{
			vv.push_back(j);
			if (regions[j].is_allelic()) assert(sp[k].first == sp[k].second);
		}
	}

	return 0;
}

int bundle_bridge::align_transcript(const map<as_pos32, int> &m, const transcript &t, vector<int> &vv)
{
	throw runtime_error("bundle_bridge::align_transcript() not used & not implemented");
	vv.clear();
	int k1 = -1;
	int k2 = -1;
	for(int k = 0; k < t.exons.size(); k++)
	{
		if(t.exons[k].second > bb.lpos)
		{
			k1 = k;
			break;
		}
	}
	for(int k = t.exons.size() - 1; k >= 0; k--)
	{
		if(t.exons[k].first < bb.rpos)
		{
			k2 = k;
			break;
		}
	}

	if(k1 > k2) return 0;
	if(k1 == -1 || k2 == -1) return 0;

	vector<PI> sp;
	sp.resize(k2 + 1);

	as_pos32 p1 = t.exons[k1].first.p32 > bb.lpos ? t.exons[k1].first : as_pos32(bb.lpos);
	as_pos32 p2 = t.exons[k2].second.p32 < bb.rpos ? t.exons[k2].second : as_pos32(bb.rpos);

	sp[k1].first = locate_region(p1);
	for(int k = k1 + 1; k <= k2; k++)
	{
		p1 = t.exons[k].first;
		map<as_pos32, int>::const_iterator it = m.find(p1);
		assert(it != m.end());
		sp[k].first = it->second;
	}

	sp[k2].second = locate_region(p2 - 1);
	for(int k = k1; k < k2; k++)
	{
		p2 = t.exons[k].second;
		map<as_pos32, int>::const_iterator it = m.find(p2);
		assert(it != m.end());
		sp[k].second = it->second - 1; 
	}

	for(int k = k1; k <= k2; k++)
	{
		assert(sp[k].first <= sp[k].second);
		if(k > k1) assert(sp[k - 1].second < sp[k].first);
		for(int j = sp[k].first; j <= sp[k].second; j++) vv.push_back(j);
	}

	return 0;
}

int bundle_bridge::index_references()
{
	ref_index.clear();
	ref_index.resize(regions.size());
	for(int k = 0; k < ref_phase.size(); k++)
	{
		vector<int> &v = ref_phase[k];
		for(int j = 0; j < v.size(); j++)
		{
			int x = v[j];
			ref_index[x].push_back(PI(k, j));
		}
	}
	return 0;
}

int bundle_bridge::locate_region_left(const map<as_pos32, int> &m, as_pos32 x)
{
	if(regions.size() == 0) return -1;

	if (x.ale != "$") 
	{
		auto it = m.find(x);
		assert(it != m.end());
		return it->second;
	}

	return locate_region(x);
}

int bundle_bridge::locate_region_right(const map<as_pos32, int> &m, as_pos32 x)
{
	if(regions.size() == 0) return -1;

	if (x.ale != "$") 
	{
		auto it = m.find(x);
		assert(it != m.end());
		return it->second;
	}

	return locate_region(x-1);
}

int bundle_bridge::locate_region(as_pos32 x)
{
	if(regions.size() == 0) return -1;
	assert (x.ale == "$");
	int k1 = 0;
	int k2 = regions.size();
	while(k1 < k2)
	{
		int m = (k1 + k2) / 2;
		region &r = regions[m];
		if(x.rightsameto(r.lpos) && x.leftsameto(r.rpos)) return m;
		else if(x < r.lpos) k2 = m;
		else k1 = m;
	}
	return -1;
}

int bundle_bridge::build_fragments()
{

	int ctp = 0;// count fragments number from paired-end reads
	int ctu = 0;// count fragments number from UMI
	int ctb = 0;// count fragments number for both

	// TODO: parameters
	int32_t max_misalignment1 = 20;
	int32_t max_misalignment2 = 10;

	fragments.clear();
	if(bb.hits.size() == 0) return 0;

	int max_index = bb.hits.size() + 1;
	if(max_index > 1000000) max_index = 1000000;

	vector< vector<int> > vv;
	vv.resize(max_index);

	// first build index
	for(int i = 0; i < bb.hits.size(); i++)
	{
		hit &h = bb.hits[i];
		if(h.isize >= 0) continue;
		if(h.vlist.size() == 0) continue;

		// do not use hi; as long as qname, pos and isize are identical
		int k = (h.qhash % max_index + h.pos % max_index + (0 - h.isize) % max_index) % max_index;
		/*
		SI si(h.qname, h.hi);
		MSI &m = vv[k];
		assert(m.find(si) == m.end());
		m.insert(PSI(si, i));
		*/
		vv[k].push_back(i);
	}

	for(int i = 0; i < bb.hits.size(); i++)
	{
		hit &h = bb.hits[i];
		if(h.paired == true) continue;
		if(h.isize <= 0) continue;
		if(h.vlist.size() == 0) continue;

		int k = (h.qhash % max_index + h.mpos % max_index + h.isize % max_index) % max_index;

		/*
		h.print();
		for(int j = 0; j < vv[k].size(); j++)
		{
			hit &z = bb.hits[vv[k][j]];
			printf(" ");
			z.print();
		}
		*/

		int x = -1;
		for(int j = 0; j < vv[k].size(); j++)
		{
			hit &z = bb.hits[vv[k][j]];
			//if(z.hi != h.hi) continue;
			if(z.paired == true) continue;
			if(z.pos != h.mpos) continue;
			if(z.isize + h.isize != 0) continue;
			if(z.qhash != h.qhash) continue;
			if(z.qname != h.qname) continue;
			x = vv[k][j];
			break;
		}

		/*
		SI si(h.qname, h.hi);
		MSI::iterator it = vv[k].find(si);
		if(it == vv[k].end()) continue;
		int x = it->second;
		*/

		//printf("HIT: i = %d, x = %d, bb.hits[i].vlist = %lu | ", i, x, bb.hits[i].vlist.size(), bb.hits[i].qname.c_str()); bb.hits[i].print();

		if(x == -1) continue;
		if(bb.hits[x].vlist.size() == 0) continue;

		fragment fr(&bb.hits[i], &bb.hits[x]);

		// ===============================
		// TODO: dit for UMI
		bb.hits[i].pi = x;
		bb.hits[x].pi = i;
		bb.hits[i].fidx = fragments.size();
		bb.hits[x].fidx = fragments.size();
		ctp += 1;
		fr.type = 0; 
		// ================================
		fr.lpos = h.pos;
		fr.rpos = bb.hits[x].rpos;

		vector<int> v1 = decode_vlist(bb.hits[i].vlist);
		vector<int> v2 = decode_vlist(bb.hits[x].vlist);
		fr.k1l = fr.h1->pos - regions[v1.front()].lpos;
		fr.k1r = regions[v1.back()].rpos - fr.h1->rpos;
		fr.k2l = fr.h2->pos - regions[v2.front()].lpos;
		fr.k2r = regions[v2.back()].rpos - fr.h2->rpos;

		fr.b1 = true;
		if(v1.size() <= 1) 
		{
			fr.b1 = false;
		}
		else if(v1.size() >= 2 && v1[v1.size() - 2] == v1.back() - 1)
		{
			if(fr.h1->rpos - regions[v1.back()].lpos > max_misalignment1 + fr.h1->nm) fr.b1 = false;
		}
		else if(v1.size() >= 2 && v1[v1.size() - 2] != v1.back() - 1)
		{
			if(fr.h1->rpos - regions[v1.back()].lpos > max_misalignment2 + fr.h1->nm) fr.b1 = false;
		}

		fr.b2 = true;
		if(v2.size() <= 1)
		{
			fr.b2 = false;
		}
		else if(v2.size() >= 2 || v2[1] == v2.front() + 1)
		{
			if(regions[v2.front()].rpos.p32 - fr.h2->pos > max_misalignment1 + fr.h2->nm) fr.b2 = false;
		}
		else if(v2.size() >= 2 || v2[1] != v2.front() + 1)
		{
			if(regions[v2.front()].rpos.p32 - fr.h2->pos > max_misalignment2 + fr.h2->nm) fr.b2 = false;
		}

		fragments.push_back(fr);

		bb.hits[i].paired = true;
		bb.hits[x].paired = true;

	}

	//printf("total bb.hits = %lu, total fragments = %lu\n", bb.hits.size(), fragments.size());
	

	// TODO
	// ===============================================
	// build fragments based on UMI
	// ===============================================
	
	vector<string> ub;
	vector< vector<int> > hlist;

	int sp = 0;

	// create ub and corresponding hlist
	for(int i = 0; i < bb.hits.size(); i++)
	{
		hit &h = bb.hits[i];
		assert(h.pos >= sp);
		sp = h.pos;

		if((h.flag & 0x4) >= 1 || h.umi == "") continue;

		bool new_umi = true;
		int ubidx = -1;
		for(int j = 0; j < ub.size(); j++)
		{
			if(h.umi == ub[j])
			{
				new_umi = false;
				ubidx = j;
				break;
			}
		}

		if(new_umi)
		{
			ub.push_back(h.umi);
			vector<int> head;
			head.push_back(i);
			hlist.push_back(head);
		}
		else hlist[ubidx].push_back(i);
	}


	/*	
	// print ub and hlist
        assert(ub.size() == hlist.size());
        for(int k = 0; k < ub.size(); k++)
        {
                printf("ub: %s, hlist # %d: (  ", ub[k].c_str(), k);
                for(int kk = 0; kk < hlist[k].size(); kk++)
                {
                        printf("hit %d: ", hlist[k][kk]);
			vector<int> v1 = decode_vlist(bb.hits[(hlist[k][kk])].vlist);
			for(int kkk = 0; kkk < v1.size(); kkk++)
			{
				printf("%d ", v1[kkk]);
			}
			printf(";");

                }
                printf(")\n");
        }
	*/

	

	/*
	//print ub
	for(int k = 0; k < ub.size(); k++)
        {
                printf("check repeat UMI %s\n", ub[k].c_str());
        }
	*/

	

	// build fragments based on hlist
	// every two consecutive hits in hlist stored as one fragment
	int hidx1 = -1;
	int hidx2 = -1;

	umiLink.clear();

	for(int i = 0; i < hlist.size(); i++)
	{
		assert(hlist[i].size() > 0);
		if(hlist[i].size() == 1) continue;

		vector<int> flist;
		flist.clear();

		for(int j = 0; j < hlist[i].size() - 1; j++)
		{
			hidx1 = hlist[i][j];
			hidx2 = hlist[i][j+1];

			// check whether has been paired
			if ( bb.hits[hidx1].pi == hidx2 && bb.hits[hidx2].pi == hidx1 && bb.hits[hidx1].paired == true && bb.hits[hidx2].paired == true)
			{
				//printf("Already exist paired-end fr: (%d, %d), fr# %d\n", hidx1, hidx2, bb.hits[hidx2].fidx);
				//assert(bb.hits[hidx1].fidx == bb.hits[hidx2].fidx);
				int fr_idx = bb.hits[hidx1].fidx;
				//printf("bothing....................fr.size() = %d, fidx = %d\n", fragments.size(), fr_idx);
				fragments[fr_idx].type = 2;
				ctb += 1;

				// TODO add fr index to umiLink
				flist.push_back(fr_idx);

				
				continue;
			}

			if(bb.hits[hidx1].vlist.size() == 0 || bb.hits[hidx2].vlist.size() == 0) continue;

			fragment fr(&bb.hits[hidx1], &bb.hits[hidx2]);
			fr.type = 1;
			ctu += 1;
			fr.lpos = bb.hits[hidx1].pos;
			fr.rpos = bb.hits[hidx2].rpos;

			vector<int> v1 = decode_vlist(bb.hits[hidx1].vlist);
			vector<int> v2 = decode_vlist(bb.hits[hidx2].vlist);
			fr.k1l = fr.h1->pos - regions[v1.front()].lpos;
			fr.k1r = regions[v1.back()].rpos - fr.h1->rpos;
			fr.k2l = fr.h2->pos - regions[v2.front()].lpos;
			fr.k2r = regions[v2.back()].rpos - fr.h2->rpos;

			fr.b1 = true;
			if(v1.size() <= 1)
			{
				fr.b1 = false;
			}
			else if(v1.size() >= 2 && v1[v1.size() - 2] == v1.back() - 1)
			{
				if(fr.h1->rpos - regions[v1.back()].lpos > max_misalignment1 + fr.h1->nm) fr.b1 = false;
			}
			else if(v1.size() >= 2 && v1[v1.size() - 2] != v1.back() - 1)
			{
				if(fr.h1->rpos - regions[v1.back()].lpos > max_misalignment2 + fr.h1->nm) fr.b1 = false;
			}

			fr.b2 = true;
			if(v2.size() <= 1)
			{
					fr.b2 = false;
			}
			else if(v2.size() >= 2 || v2[1] == v2.front() + 1)
			{
					if(regions[v2.front()].rpos.p32 - fr.h2->pos > max_misalignment1 + fr.h2->nm) fr.b2 = false;
			}
			else if(v2.size() >= 2 || v2[1] != v2.front() + 1)
			{
					if(regions[v2.front()].rpos.p32 - fr.h2->pos > max_misalignment2 + fr.h2->nm) fr.b2 = false;
			}

			fragments.push_back(fr);
			bb.hits[hidx1].paired = true;
			bb.hits[hidx2].paired = true;

			int cur_fidx = fragments.size() - 1;
			flist.push_back(cur_fidx);

		}

		umiLink.push_back(flist);

	}

	return 0;
}

int32_t bundle_bridge::compute_aligned_length(int32_t k1l, int32_t k2r, const vector<int>& v)
{
	if(v.size() == 0) return 0;
	int32_t flen = 0;
	for(int i = 0; i < v.size(); i++)
	{
		int k = v[i];
		flen += regions[k].rpos - regions[k].lpos;
	}
	return flen - k1l - k2r;
}

int bundle_bridge::print(int index)
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

	printf("tid = %d, #hits = %lu, #fragments = %lu, #ref-trsts = %lu, range = %s:%d-%d, orient = %c (%d, %d, %d)\n",
			bb.tid, bb.hits.size(), fragments.size(), ref_trsts.size(), bb.chrm.c_str(), bb.lpos, bb.rpos, bb.strand, n0, np, nq);

	if(verbose <= 1) return 0;

	// print junctions 
	for(int i = 0; i < junctions.size(); i++)
	{
		junctions[i].print(bb.chrm, i);
	}

	// print bb.hits
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

	printf("\n");

	return 0;
}

vector<int32_t> bundle_bridge::build_accumulate_length(const vector<int> &v)
{
	int32_t x = 0;
	vector<int32_t> acc;
	acc.resize(v.size());
	for(int i = 0; i < v.size(); i++)
	{
		x += regions[v[i]].rpos - regions[v[i]].lpos;
		acc[i] = x;
	}
	return acc;
}

vector<as_pos32> bundle_bridge::get_aligned_intervals(fragment &fr)
{
	vector<as_pos32> vv;
	//if(fr.h1->bridged == false) return vv;
	//if(fr.h2->bridged == false) return vv;
	if(fr.paths.size() != 1) return vv;
	assert(fr.paths[0].type == 1 || fr.paths[0].type == 2);
	//if(fr.paths[0].type != 1) return vv;

	vector<as_pos32> v = get_splices(fr);
	if(v.size() >= 1 && fr.h1->pos >= v.front()) return vv;
	if(v.size() >= 1 && fr.h2->rpos <= v.back()) return vv;

	v.insert(v.begin(), fr.h1->pos);
	v.push_back(fr.h2->rpos);
	return v;
}

vector<as_pos32> bundle_bridge::get_splices(fragment &fr)
{
	vector<as_pos32> vv;
	//if(fr.h1->bridged == false) return vv;
	//if(fr.h2->bridged == false) return vv;
	if(fr.paths.size() != 1) return vv;
	assert(fr.paths[0].type == 1 || fr.paths[0].type == 2);
	//if(fr.paths[0].type != 1) return vv;

	vector<int> v = decode_vlist(fr.paths[0].v);
	if(v.size() <= 0) return vv;

	for(int i = 0; i < v.size() - 1; i++)
	{
		as_pos32 pp = regions[v[i + 0]].rpos;
		as_pos32 qq = regions[v[i + 1]].lpos;
		if(pp >= qq) continue;
		vv.push_back(pp);
		vv.push_back(qq);
	}
	return vv;
}
