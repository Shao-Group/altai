/*
Part of Altai
(c) 2022 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include "phaser.h"

int phaser::phase()
{
	write_phasing_profile(bd.gr, bd.hs, bd);
	write_phasing_profile_classify(bd.gr, bd.hs, bd);
	return 0;
}

int phaser::write_phasing_profile(...) //TODO:
{
	MEI e2i;// edge map, from edge to index
    VE i2e;// edge map, from index to ed
	gr.get_edge_indices(i2e, e2i);
	IN_ILP_MODULE = true;

	hs.build(gr, e2i);

	map <PI32, map<vector<int>, int> > as_exon_phase_map; // AS exon to map <edge_list, count>

	PI32 last_pos = make_pair(as_pos32(-1),as_pos32(-1));
	for(int i = 0; i < hs.edges.size(); i++)
	{
		const vector<int> &v = hs.edges[i];
		int c = hs.ecnts[i];

		for (auto v_0: v)
		{
			edge_descriptor d= i2e[v_0];

			vertex_info v1 = gr.vinf[d->source()];
			vertex_info v2 = gr.vinf[d->target()];
			
			PI32 source_exon = make_pair(v1.lpos, v1.rpos);
			PI32 target_exon = make_pair(v2.lpos, v2.rpos);

			if (source_exon != last_pos)
			{
				if (v1.lpos.ale != "$" || v1.rpos.ale != "$") 
				{
					auto _p = as_exon_phase_map.find(source_exon);
					if (_p == as_exon_phase_map.end()) 
					{
						map<vector<int>, int> _f;
						_f.insert(make_pair(v, c));
						as_exon_phase_map.insert(make_pair(source_exon, _f));
					}
					else
					{
						auto _f0 = _p->second.find(v);
						if (_f0 == _p->second.end())
						{
							_p->second.insert(make_pair(v,c));
						}
						else
						{
							_f0->second += c;
						}
					}
				}
			}
			if (target_exon != source_exon)
			{
				if (v2.lpos.ale != "$" || v2.rpos.ale != "$") 
				{
					auto _p = as_exon_phase_map.find(target_exon);
					if (_p == as_exon_phase_map.end()) 
					{
						map<vector<int>, int> _f;
						_f.insert(make_pair(v, c));
						as_exon_phase_map.insert(make_pair(target_exon, _f));
					}
					else
					{
						auto _f0 = _p->second.find(v);
						if (_f0 == _p->second.end())
						{
							_p->second.insert(make_pair(v,c));
						}
						else
						{
							_f0->second += c;
						}
					}
				}
			}
			last_pos = target_exon;
		}
	}

	map <PI32, map<edge_descriptor, int> > as_exon_phase_map2; // AS exon to map <edge, count>
	for (auto as_ex_p: as_exon_phase_map)
	{
		PI32 a = as_ex_p.first;
		map<vector<int>, int> & edges_w_count_map = as_ex_p.second; 
		for (auto edges_w_count: edges_w_count_map)
		{
			int c = edges_w_count.second;
			vector<int> _j = edges_w_count.first;
			for (auto j : _j)
			{
				edge_descriptor d= i2e[j];

				if (as_exon_phase_map2.find(a) == as_exon_phase_map2.end())
				{
					map<edge_descriptor, int> f;
					f.insert(make_pair(d, c));
					as_exon_phase_map2.insert(make_pair(a, f));
				}
				else
				{
					map<edge_descriptor, int> &f = as_exon_phase_map2.find(a)->second;
					if (f.find(d) == f.end())
					{
						f.insert(make_pair(d, c));
					}
					else
					{
						int _c = f.find(d)->second + c;
						f.insert(make_pair(d, _c));
					}
				}
			}
		}		
	}


	// FIXME: incorrect
	/*
	map <PI32, map<edge_descriptor, int> > as_exon_phase_map3; // AS exon to map <unique edge, count>, unique edge is not phased to non-compatible variant
	vector <PI32> ky;
	for (auto _ky: as_exon_phase_map2) ky.push_back(_ky.first);
	
	if (ky.size() == 0)  return 0;
	sort(ky.begin(), ky.end()); // sort s.t. noncompatible var are adjacent
	PI32 last_exon = ky[0];
	map <edge_descriptor, int> & last_exon_edges = as_exon_phase_map2[last_exon];

	for (int i = 1; i < ky.size(); i++)
	{
		PI32 as_ex = ky[i];		
		map <edge_descriptor, int> & exon_edges = as_exon_phase_map2[as_ex];
		if (!(as_ex.first.samepos(last_exon.first) && as_ex.second.samepos(last_exon.second))) // if same exon pos, add diff edges; otherwise, no way to distinguish
		{
			map<edge_descriptor, int> _f;
			_f.clear();
			for (auto _e : (last_exon_edges)) // TODO: debug copy pointer, use index instead of edge_descriptor 
			{
				if (exon_edges->find(*(_e.first)) == exon_edges->end()) {_f.insert(_e);}
			}
			as_exon_phase_map3.insert(make_pair(last_exon, _f));

			map<edge_descriptor, int> _f2;
			_f2.clear();
			for (auto _e : (exon_edges))
			{
				if (last_exon_edges->find(*(_e.first)) == last_exon_edges->end()) {_f2.insert(_e);}
			}
			as_exon_phase_map3.insert(make_pair(as_ex, _f2));
			continue;
		}		
		last_exon = as_ex;
		last_exon_edges = exon_edges;	
	}
	*/

	// write
	ofstream phase_file;
	phase_file.open(phasing_profile.c_str(), fstream::app);
	phase_file << "Phasing profile of bundle: " << bd.chrm << ":" << bd.lpos << "-" << bd.rpos << endl;	
	gr.print(); // stdout
	phase_file << "\tsize as_exon_phase_map = " << as_exon_phase_map.size() <<":[";
	for (auto as_ex_p: as_exon_phase_map)
	{
		phase_file << as_ex_p.second.size()<<", ";
	}
	phase_file << "]" << endl;

	for (auto as_ex_p2: as_exon_phase_map2) //map <PI32, map<edge_descriptor, int> >
	{
		phase_file << "ASExon "<< as_ex_p2.first.first.p32 << as_ex_p2.first.first.ale << "-" << as_ex_p2.first.second.p32 << as_ex_p2.first.second.ale;
		phase_file << " PhasingEdgesList[" << as_ex_p2.second.size() << "] = {";
		for (auto p: as_ex_p2.second)
		{
			phase_file << p.first->source() <<"->" << p.first->target(); // edge index
			vertex_info v1 = gr.vinf[p.first->source()];
			vertex_info v2 = gr.vinf[p.first->target()];
			phase_file << "(" << v1.rpos.p32 << v1.rpos.ale <<"->" << v2.lpos.p32 << v2.lpos.ale << ")" ; // edge as_pos32
			phase_file << "[" << p.second << "]" << ", "; // counts

		}
		phase_file << "}" << endl;
	}
	phase_file.close();


	// // write2
	// ofstream phase_file2;
	// string phasing_profile2 = phasing_profile + "2";
	// phase_file2.open(phasing_profile2.c_str(), fstream::app);
	// phase_file2 << "Phasing profile2 of bundle: " << bd.chrm << ":" << bd.lpos << "-" << bd.rpos << endl;	
	// phase_file2 << "\tsize as_exon_phase_map3 = " << as_exon_phase_map3.size() <<":[";
	// for (auto as_ex_p: as_exon_phase_map3)
	// {
	// 	phase_file2 << as_ex_p.second.size()<<", ";
	// }
	// phase_file2 << "]" << endl;

	// for (auto as_ex_p3: as_exon_phase_map3) //map <PI32, map<edge_descriptor, int> >
	// {
	// 	phase_file2 << "ASExon "<< as_ex_p3.first.first.p32 << as_ex_p3.first.first.ale << "-" << as_ex_p3.first.second.p32 << as_ex_p3.first.second.ale;
	// 	phase_file2 << " AllelicPhasingEdgesList[" << as_ex_p3.second.size() << "] = {";
	// 	for (auto p3: as_ex_p3.second)
	// 	{
	// 		phase_file2 << p3.first->source() <<"->" << p3.first->target(); // edge index
	// 		// vertex_info v1 = gr.vinf[p.first->source()];
	// 		// vertex_info v2 = gr.vinf[p.first->target()];
	// 		// phase_file2 << "(" << v1.rpos.p32 << v1.rpos.ale <<"->" << v2.lpos.p32 << v2.lpos.ale << ")" ; // edge as_pos32
	// 		phase_file2 << "[" << p3.second << "]" << ", "; // counts

	// 	}
	// 	phase_file2 << "}" << endl;
	// }
	// phase_file2.close();

	return 0;
}

/*  
**  Three classes of profiling: 
**	1. directly phased 
**  2. AS junctions phased 
**  3a. not affected by decomposition -- assign by abundance
**  3b. affected by decomposition -- new algorithm 
*/
int assephasermbler::write_phasing_profile_classify(...) //TODO:
{
	MEI e2i;// edge map, from edge to index
    VE i2e;// edge map, from index to ed
	gr.get_edge_indices(i2e, e2i);
	IN_ILP_MODULE = true;

	hs.build(gr, e2i);

	map <int, map<vector<int>, int> > as_exon_phase_map; // AS_exon_index to map <edge_list, count>
	map <int, PI32> i2exon;  // node_index to AS exon
	map <PI32, int> exon2i;  // AS exon to node_index

	PI32 last_pos = make_pair(as_pos32(-1),as_pos32(-1));
	for(int i = 0; i < hs.edges.size(); i++)
	{
		const vector<int> &v = hs.edges[i];
		int c = hs.ecnts[i];

		for (int v_0: v)
		{
			edge_descriptor d= i2e[v_0];

			vertex_info v1 = gr.vinf[d->source()];
			vertex_info v2 = gr.vinf[d->target()];
			
			PI32 source_exon = make_pair(v1.lpos, v1.rpos);
			PI32 target_exon = make_pair(v2.lpos, v2.rpos);

			i2exon[d->source()] = source_exon;
			i2exon[d->target()] = target_exon;
			exon2i[source_exon] = d->source();
			exon2i[target_exon] = d->target();

			if (source_exon != last_pos)
			{
				if (v1.lpos.ale != "$" || v1.rpos.ale != "$") 
				{
					auto _p = as_exon_phase_map.find(d->source());
					if (_p == as_exon_phase_map.end()) 
					{
						map<vector<int>, int> _f;
						_f.insert(make_pair(v, c));
						as_exon_phase_map.insert(make_pair(d->source(), _f));
					}
					else
					{
						auto _f0 = _p->second.find(v);
						if (_f0 == _p->second.end())
						{
							_p->second.insert(make_pair(v,c));
						}
						else
						{
							_f0->second += c;
						}
					}
				}
			}
			if (target_exon != source_exon)
			{
				if (v2.lpos.ale != "$" || v2.rpos.ale != "$") 
				{
					auto _p = as_exon_phase_map.find(d->target());
					if (_p == as_exon_phase_map.end()) 
					{
						map<vector<int>, int> _f;
						_f.insert(make_pair(v, c));
						as_exon_phase_map.insert(make_pair(d->target(), _f));
					}
					else
					{
						auto _f0 = _p->second.find(v);
						if (_f0 == _p->second.end())
						{
							_p->second.insert(make_pair(v,c));
						}
						else
						{
							_f0->second += c;
						}
					}
				}
			}
			last_pos = target_exon;
		}
	}

	map <int, map<int, int> > as_exon_phase_map2; // AS_exon_index to map <edge_idx, count>
	for (auto as_ex_p: as_exon_phase_map)
	{
		int a = as_ex_p.first;
		map<vector<int>, int> & edges_w_count_map = as_ex_p.second; 
		for (auto edges_w_count: edges_w_count_map)
		{
			int c = edges_w_count.second;
			vector<int> _j = edges_w_count.first;
			for (auto j : _j)
			{
				int d= j;

				if (as_exon_phase_map2.find(a) == as_exon_phase_map2.end())
				{
					map<int, int> f;
					f.insert(make_pair(d, c));
					as_exon_phase_map2.insert(make_pair(a, f));
				}
				else
				{
					map<int, int> &f = as_exon_phase_map2.find(a)->second;
					if (f.find(d) == f.end())
					{
						f.insert(make_pair(d, c));
					}
					else
					{
						int _c = f.find(d)->second + c;
						f.insert(make_pair(d, _c));
					}
				}
			}
		}		
	}

	// get heterozygous AS exons
	vector<PI32> _AS_exons;
	set<PI32> hetero_AS_exons;
	vector<PI32> hetero_AS_exons_sorted;
	vector<int> hetero_AS_exons_sorted_index;
	set<int> hetero_AS_exons_index_set;
	_AS_exons.clear();
	hetero_AS_exons.clear();
	hetero_AS_exons_sorted.clear();
	hetero_AS_exons_sorted_index.clear();
	for (auto i: exon2i)  _AS_exons.push_back(i.first);
	sort(_AS_exons.begin(), _AS_exons.end());
	for (int i = 1; i < _AS_exons.size(); i++)
	{
		if (_AS_exons[i-1].first.samepos(_AS_exons[i].first) && _AS_exons[i-1].second.samepos(_AS_exons[i].second))
		{
			hetero_AS_exons.insert(_AS_exons[i-1]);
			hetero_AS_exons.insert(_AS_exons[i]);
		}
	}
	for (PI32 i: hetero_AS_exons) hetero_AS_exons_sorted.push_back(i);
	sort(hetero_AS_exons_sorted.begin(), hetero_AS_exons_sorted.end());
	for (PI32 i: hetero_AS_exons_sorted) hetero_AS_exons_sorted_index.push_back(exon2i[i]);
	hetero_AS_exons_index_set.insert(hetero_AS_exons_sorted_index.begin(), hetero_AS_exons_sorted_index.end());

	// get cat 1 AS exons
	map<int, set<int> > exon_index_to_index_set;
	set< set<int> > set_of_index_set; // directly phased AS exon indexes, each size >= 2
	set<int> union_all_index_set;
	for (int i = 0; i < hs.edges.size(); i++)
	{
		const vector<int> &v = hs.edges[i];
		set<int> n;
		for (int j: v) // hetero AS nodes only
		{
			int so = i2e[j]->source();
			int si = i2e[j]->target();
			if (hetero_AS_exons_index_set.find(so) != hetero_AS_exons_index_set.end()) n.insert(so);
			if (hetero_AS_exons_index_set.find(si) != hetero_AS_exons_index_set.end()) n.insert(si);
		}
		set<int> _s(n.begin(), n.end());

		vector< set<int> > sets_to_union;
		sets_to_union.push_back(_s);
		set<int> unioned_set;
		for (int j: v) sets_to_union.push_back(exon_index_to_index_set[j]);
		for (set<int> k: sets_to_union) unioned_set.insert(k.begin(), k.end());
		for (int j: unioned_set) exon_index_to_index_set[j] = unioned_set;
	}
	for (auto i: exon_index_to_index_set) 
	{
		if (i.second.size() >=2) 
		{
			set_of_index_set.insert(i.second);
			union_all_index_set.insert(i.second.begin(), i.second.end());
		}
	}
	// get cat 3a AS exons
	// set<int> to_phase;
	set <int> cat3;
	map<int, set<int> > cat3_index_to_index_set;
	for (int i = 0; i < hetero_AS_exons_sorted_index.size(); i++)
	{
		int idx = hetero_AS_exons_sorted_index[i];
		// fall into cat 1
		if (union_all_index_set.find(i) != union_all_index_set.end()) continue; 
		// still fall into cat 1 because preceding counter-allele falls into it
		if (i - 1 < hetero_AS_exons_sorted_index.size()
			&& union_all_index_set.find(hetero_AS_exons_sorted_index[i-1]) != union_all_index_set.end())
		{
			PI32 pre_exon = i2exon[hetero_AS_exons_sorted_index[i-1]];
			PI32 this_exon = i2exon[idx];
			if (pre_exon.first.samepos(this_exon.first) && pre_exon.second.samepos(this_exon.second))
			{
				union_all_index_set.insert(idx);
				continue;
			}
		}
		// still fall into cat 1 because next counter-allele falls into it
		if (i + 1 < hetero_AS_exons_sorted_index.size()
			&& union_all_index_set.find(hetero_AS_exons_sorted_index[i + 1]) != union_all_index_set.end())
		{
			PI32 next_exon = i2exon[hetero_AS_exons_sorted_index[i + 1]];
			PI32 this_exon = i2exon[idx];
			if (next_exon.first.samepos(this_exon.first) && next_exon.second.samepos(this_exon.second))
			{
				union_all_index_set.insert(idx);
				continue;
			}
		}
		// to_phase.insert(idx);

		// search a path between this idx & all other idx
		set<int> v;
		for (int j = 0; j < hetero_AS_exons_sorted_index.size(); j++)
		{
			bool all_path_threading = false;
			if (j == i) continue;
			int idx_may_phase = hetero_AS_exons_sorted_index[j];

			if (i2exon[idx].first >= i2exon[idx_may_phase].second)
				if (gr.check_all_path(idx_may_phase, idx)) all_path_threading = true;
			else if (i2exon[idx].second <= i2exon[idx_may_phase].first)
				if (gr.check_all_path(idx, idx_may_phase)) all_path_threading = true;
			
			if (all_path_threading) 
			{
				cat3.insert(idx);
				v.insert(idx);
				v.insert(idx_may_phase);
			}
		}
		vector< set<int> > sets_to_union;
		sets_to_union.push_back(v);
		set<int> unioned_set;
		for (int j: v) sets_to_union.push_back(cat3_index_to_index_set[j]);
		for (set<int> k: sets_to_union) unioned_set.insert(k.begin(), k.end());
		for (int j: unioned_set) cat3_index_to_index_set[j] = unioned_set;
	}


	// write
	ofstream phase_file;
	string phasing_profile_car = phasing_profile + ".cat";
	phase_file.open(phasing_profile_car.c_str(), fstream::app);
	phase_file << "Phasing profile of bundle: " << bd.chrm << ":" << bd.lpos << "-" << bd.rpos << endl;	
	gr.print(); // stdout
	
	int largest_size = 0;
	for (set<int> i: set_of_index_set) 
		if (i.size() > largest_size) largest_size = i.size();
	int largest_size_cat3 = 0;
	for (auto i: cat3_index_to_index_set) 
		if (i.second.size() > largest_size_cat3) largest_size_cat3 = i.second.size();
	phase_file << "# total hetero variant " << hetero_AS_exons_sorted_index.size() << endl ; // edge as_pos32
	if (hetero_AS_exons_sorted_index.size() != 0)
	{
		phase_file << "# cat 1 variant (largest phase_set) " << largest_size << " | " << float(largest_size)/hetero_AS_exons_sorted_index.size() * 100 << "%" << endl;
		phase_file << "# cat 1 variant (union) " << union_all_index_set.size() << " | " << float(union_all_index_set.size())/hetero_AS_exons_sorted_index.size() *100 << "%" << endl;
		phase_file << "# cat 3a variant (union & exclud cat1) " << cat3.size() << " | " << float(cat3.size())/ hetero_AS_exons_sorted_index.size() * 100 << "%" << endl;
		phase_file << "# cat 3a variant (largest & includ cat1) " << largest_size_cat3 << " | " << float(largest_size_cat3)/ hetero_AS_exons_sorted_index.size() * 100 << "%" << endl;
	}	
	phase_file << endl;

	phase_file.close();
	return 0;
}