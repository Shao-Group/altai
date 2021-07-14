/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include "filter.h"
#include "config.h"
#include "gurobi_c++.h"
#include <cassert>
#include <algorithm>
#include "as_pos32.hpp"

filter::filter(const vector<transcript> &v)
	:trs(v)
{}

int filter::simple_phase_set_by_coverage()
{
	int NUM_PHASE_SET = 2; // for human, number of phase set is 2 (diploid)

	vector<PI32> as_exons;
	for(int i = 0; i < trs.size(); ++i)
	{
		for (PI32 j: trs[i].exons)
		{
			if (j.first.ale !="$") 
			{
				as_exons.push_back(j);
			}
		}
	}
	sort(as_exons.begin(), as_exons.end());
	map<PI32, int> as_exons_enum;
	for(int i = 0; i < as_exons.size(); ++i) 
	{
		as_exons_enum.insert(make_pair(as_exons[i], i));
	}

	map<int, vector<int> > trs_as_exon_map;
	for(int i = 0; i < trs.size(); ++i)
	{
		vector<int> e;
		for (PI32 j: trs[i].exons)
		{
			if (j.first.ale !="$") 
			{
				e.push_back(as_exons_enum.find(j)->second);
			}
		}
		trs_as_exon_map.insert(make_pair(i, e));
	}

	
	try 
	{
		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);
		
		// Create variables, all variables are binary
		// x_{i,k} 		- whether i-th AS-exon is in k-th phase set 
		// y_m 			- whether m-th transcript is retained after filter
		// z_{m,k} 		- whether k-th phase set satisfies m-th transcript
		GRBVar x[as_exons.size()][NUM_PHASE_SET];
		GRBVar y[trs.size()];
		GRBVar z[trs.size()][NUM_PHASE_SET];
		GRBLinExpr obj = GRBLinExpr();
		for (int i = 0; i < as_exons.size(); ++i)
		{
			for (int k = 0; k < NUM_PHASE_SET; ++k)
			{
				x[i][k] = model.addVar(0, 1, 0, GRB_BINARY);
			}
		}
		for (int m = 0; m < trs.size(); ++m) 
		{
			y[m] = model.addVar(0, 1, 0, GRB_BINARY);
			obj += trs[m].coverage * y[m];
		}
		for (int m = 0; m < trs.size(); ++m) 
		{
			for (int k = 0; k < NUM_PHASE_SET; ++k)
			{
				z[m][k] = model.addVar(0, 1, 0, GRB_BINARY);
			}
		}
		if (DEBUG_MODE_ON) cout << "Variables Created successfully\n"; //CLEAN:

		// Set Objective 
		model.setObjective(obj, GRB_MAXIMIZE);
		if (DEBUG_MODE_ON) cout << "Obj set successfully\n"; //CLEAN:

		// Genome phase set constraint
		// Each PS has at most one AS vertex at each loci
		for (int i = 0; i < as_exons.size(); ++i)
		{
			as_pos32 l = as_exons[i].first;
			as_pos32 r = as_exons[i].second;
			for (int j = 1; j < as_exons.size() - i; ++j)
			{
				as_pos32 a = as_exons[i + j].first;
				as_pos32 b = as_exons[i + j].second;

				assert(l.leftsameto(a));
				if (a.rightsameto(r)) break;
				else 
				{
					// x[i] and x[i + j] overlap, add constrait
					for (int k = 0; k < NUM_PHASE_SET; ++k)
					{
						model.addConstr(x[i][k] + x[i + j][k] <= 1);
					}
				}
			}
		}
		if (DEBUG_MODE_ON) cout << "genome PS constrint successfully\n"; //CLEAN:
		// Transcript phase set constraint
		// A retained transcript must be satisfied by at least one phase set
		int m = 0;
		for (auto it = trs_as_exon_map.begin(); it != trs_as_exon_map.end(); ++it)
		{
			vector<int> &v = it->second;
			int n = v.size();
			if (n >= 1)  // otherwise y_m has no AS-exons and has no constraints
			{
				for (int k = 0; k < NUM_PHASE_SET; ++k)
				{
					GRBLinExpr expr_xz = GRBLinExpr();  
					for (int j = 0; j < n; ++j)
					{
						expr_xz += x[v[j]][k] / double(n); 
					}
					expr_xz -= z[m][k];
					model.addConstr(expr_xz, GRB_GREATER_EQUAL, 0);  // Sum_i(x[i][k]) / n - z[m][k] >= 0
				}

				GRBLinExpr expr_zy = GRBLinExpr();
				for (int k = 0; k < NUM_PHASE_SET; ++k)
				{
					expr_zy += z[m][k];
				}
				expr_zy -= y[m];
				model.addConstr(expr_zy, GRB_GREATER_EQUAL, 0);  // Sum_k(z[m][k]) - y[m] >= 0 
			}
			++m;
		}
		assert(m == trs_as_exon_map.size());
		cout << "transcript PS constrint successfully\n"; //CLEAN:

		// Optimize
		model.optimize();
		vector<transcript> v;
		for (int m = 0; m < trs.size(); ++m)
		{
			if (y[m].get(GRB_DoubleAttr_X) >= 1) v.push_back(trs[m]);
			cout << m << " y_m: " <<y[m].get(GRB_DoubleAttr_X) << endl;
			assert(y[m].get(GRB_DoubleAttr_X) == 0 || y[m].get(GRB_DoubleAttr_X) == 1); //CLEAN:
			for (int k = 0; k < NUM_PHASE_SET; ++k) 
			{
				cout << m << ' ' <<k << " z_mk: " <<z[m][k].get(GRB_DoubleAttr_X) << endl;
			}
		}
		/* // print x_ik values
		for (int k = 0; k < NUM_PHASE_SET; ++k) 
		{
			for (int i = 0; i < as_exons.size(); ++i)
			{
				cout << i << ' ' <<k << " x_ik: " <<x[i][k].get(GRB_DoubleAttr_X) << endl;
			}
		}
		*/
		trs = v;
	}
	catch (GRBException e)
	{
		cout << "Error code = " << e.getErrorCode() << endl;
    	cout << e.getMessage() << endl;
		throw BundleError();
	}

	return 0;
}

int filter::simple_phase_set_by_variant_number()
{
	int NUM_PHASE_SET = 2; // for human, number of phase set is 2 (diploid)

	vector<PI32> as_exons;
	for(int i = 0; i < trs.size(); ++i)
	{
		for (PI32 j: trs[i].exons)
		{
			if (j.first.ale !="$") 
			{
				as_exons.push_back(j);
			}
		}
	}
	sort(as_exons.begin(), as_exons.end());
	map<PI32, int> as_exons_enum;
	for(int i = 0; i < as_exons.size(); ++i) 
	{
		as_exons_enum.insert(make_pair(as_exons[i], i));
	}

	map<int, vector<int> > trs_as_exon_map;
	for(int i = 0; i < trs.size(); ++i)
	{
		vector<int> e;
		for (PI32 j: trs[i].exons)
		{
			if (j.first.ale !="$") 
			{
				e.push_back(as_exons_enum.find(j)->second);
			}
		}
		trs_as_exon_map.insert(make_pair(i, e));
	}

	map<int, vector<int> > as_exon_trs_map;
	for (auto it = trs_as_exon_map.begin(); it != trs_as_exon_map.end(); ++it)
	{
		vector<int> asexons = it->second;
		for (auto jt = asexons.begin(); jt != asexons.end(); ++jt)
		{
			if (as_exon_trs_map.find(*jt) == as_exon_trs_map.end())
			{
				vector<int> v;
				v.push_back(it->first);
				as_exon_trs_map.insert(make_pair(*jt, v));
			}
			else 
			{
				as_exon_trs_map.find(*jt)->second.push_back(it->first);
			}
		}
	}

	
	try 
	{
		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);
		
		// Create variables, all variables are binary
		// x_{i,k} 		- whether i-th AS-exon is in k-th phase set 
		// v_i	 		- whether i-th AS-exon is retained after filter
		// y_m 			- whether m-th transcript is retained after filter
		// z_{m,k} 		- whether k-th phase set satisfies m-th transcript
		GRBVar x[as_exons.size()][NUM_PHASE_SET];
		GRBVar vi[as_exons.size()];
		GRBVar y[trs.size()];
		GRBVar z[trs.size()][NUM_PHASE_SET];
		GRBLinExpr obj = GRBLinExpr();
		double total_cov;
		for (int m = 0; m < trs.size(); ++m) 
		{
			total_cov += trs[m].coverage;
		}
		for (int i = 0; i < as_exons.size(); ++i)
		{
			for (int k = 0; k < NUM_PHASE_SET; ++k)
			{
				x[i][k] = model.addVar(0, 1, 0, GRB_BINARY);
			}
		}
		for (int i = 0; i < as_exons.size(); ++i) 
		{
			vi[i] = model.addVar(0, 1, 0, GRB_BINARY);
			obj += vi[i];
		}
		for (int m = 0; m < trs.size(); ++m) 
		{
			y[m] = model.addVar(0, 1, 0, GRB_BINARY);
			obj += trs[m].coverage * y[m] / total_cov;
		}
		for (int m = 0; m < trs.size(); ++m) 
		{
			for (int k = 0; k < NUM_PHASE_SET; ++k)
			{
				z[m][k] = model.addVar(0, 1, 0, GRB_BINARY);
			}
		}

		// Set Objective y[m].coverage * y[m]
		model.setObjective(obj, GRB_MAXIMIZE);

		// Genome phase set constraint
		// Each PS has at most one AS vertex at each loci
		for (int i = 0; i < as_exons.size(); ++i)
		{
			as_pos32 l = as_exons[i].first;
			as_pos32 r = as_exons[i].second;
			for (int j = 1; j < as_exons.size() - i; ++j)
			{
				as_pos32 a = as_exons[i + j].first;
				as_pos32 b = as_exons[i + j].second;

				assert(l.leftsameto(a));
				if (a.rightsameto(r)) break;
				else 
				{
					// x[i] and x[i + j] overlap, add constrait
					for (int k = 0; k < NUM_PHASE_SET; ++k)
					{
						model.addConstr(x[i][k] + x[i + j][k] <= 1);
					}
				}
			}
		}
		if (DEBUG_MODE_ON) cout << "genome PS constrint successfully\n"; //CLEAN:
		// Transcript phase set constraint
		// A retained transcript must be satisfied by at least one phase set
		int m = 0;
		for (auto it = trs_as_exon_map.begin(); it != trs_as_exon_map.end(); ++it)
		{
			vector<int> &v = it->second;
			int n = v.size();
			if (n >= 1)  // otherwise y_m has no AS-exons and has no constraints
			{
				for (int k = 0; k < NUM_PHASE_SET; ++k)
				{
					GRBLinExpr expr_xz = GRBLinExpr();  
					for (int j = 0; j < n; ++j)
					{
						expr_xz += x[v[j]][k] / double(n); 
					}
					expr_xz -= z[m][k];
					model.addConstr(expr_xz, GRB_GREATER_EQUAL, 0);  // Sum_i(x[i][k]) / n - z[m][k] >= 0
				}

				GRBLinExpr expr_zy = GRBLinExpr();
				for (int k = 0; k < NUM_PHASE_SET; ++k)
				{
					expr_zy += z[m][k];
				}
				expr_zy -= y[m];
				model.addConstr(expr_zy, GRB_GREATER_EQUAL, 0);  // Sum_k(z[m][k]) - y[m] >= 0 
			}
			++m;
		}
		assert(m == trs_as_exon_map.size());
		if(DEBUG_MODE_ON) cout << "transcript PS constrint successfully\n"; //CLEAN:


		for (auto it = as_exon_trs_map.begin(); it != as_exon_trs_map.end(); ++it)
		{
			GRBLinExpr expr_v = GRBLinExpr(); 
			for (int ii = 0; ii < it->second.size(); ii ++)
			{
				expr_v += y[it->second[ii]];
			}
			model.addConstr(expr_v, GRB_GREATER_EQUAL, vi[it->first]);
		}
		

		// Optimize
		model.optimize();
		vector<transcript> v;
		for (int m = 0; m < trs.size(); ++m)
		{
			if (y[m].get(GRB_DoubleAttr_X) >= 1) v.push_back(trs[m]);
			cout << m << " y_m: " <<y[m].get(GRB_DoubleAttr_X) << endl;
			assert(y[m].get(GRB_DoubleAttr_X) == 0 || y[m].get(GRB_DoubleAttr_X) == 1); //CLEAN:
			for (int k = 0; k < NUM_PHASE_SET; ++k) 
			{
				cout << m << ' ' <<k << " z_mk: " <<z[m][k].get(GRB_DoubleAttr_X) << endl;
			}
		}
		/* // print x_ik values
		for (int k = 0; k < NUM_PHASE_SET; ++k) 
		{
			for (int i = 0; i < as_exons.size(); ++i)
			{
				cout << i << ' ' <<k << " x_ik: " <<x[i][k].get(GRB_DoubleAttr_X) << endl;
			}
		}
		*/
		trs = v;
	}
	catch (GRBException e)
	{
		cout << "Error code = " << e.getErrorCode() << endl;
    	cout << e.getMessage() << endl;
		throw BundleError();
	}

	return 0;
}

int filter::keep_as_transcripts_only()
{
	vector<transcript> v;
	for(int i = 0; i < trs.size(); i++ )
	{
		bool is_as_transcript = false;
		transcript t = trs[i];
		int e = t.exons.size();
		for (int k = 0; k < e; k++) 
		{	
			string a = t.exons[k].first.ale;
			if (a != "$") 
			{
				is_as_transcript = true;
				break;
			}
			string b = t.exons[k].second.ale;
			if (b != "$") 
			{
				is_as_transcript = true;
				break;
			}
		}

		if (is_as_transcript)
		{
			v.push_back(t);
		}
	}
}

int filter::filter_length_coverage()
{
	vector<transcript> v;
	for(int i = 0; i < trs.size(); i++)
	{
		int e = trs[i].exons.size();
		int minl = min_transcript_length_base + e * min_transcript_length_increase;
		if(trs[i].length() < minl) continue;
		if(e == 1 && trs[i].coverage < min_single_exon_coverage) continue;
		if(e >= 2 && trs[i].coverage < min_transcript_coverage) continue;
		v.push_back(trs[i]);
	}
	trs = v;
	return 0;
}

int filter::remove_nested_transcripts()	//TODO: not used, to make AS compatible
{
	set<int> s;
	for(int i = 0; i < trs.size(); i++)
	{
		vector<PI32> v = trs[i].exons;
		if(v.size() <= 1) continue;
		double w1 = trs[i].coverage;
		bool b = false;
		for(int k = 1; k < v.size(); k++)
		{
			as_pos32 p = v[k - 1].second;
			as_pos32 q = v[k - 0].first;

			for(int j = 0; j < trs.size(); j++)
			{
				if(trs[j].exons.size() <= 1) continue;
				PI32 pq = trs[j].get_bounds();
				double w2 = trs[j].coverage;

				if(w2 >= w1 && pq.first.rightto(p) && pq.second.leftto(q))
				{
					b = true;
					break;
				}
			}
			if(b == true) break;
		}
		if(b == true) s.insert(i);
	}

	vector<transcript> v;
	for(int i = 0; i < trs.size(); i++)
	{
		if(s.find(i) != s.end()) continue;
		v.push_back(trs[i]);
	}

	trs = v;
	return 0;
}

int filter::join_single_exon_transcripts()
{
	while(true)
	{
		bool b = join_transcripts();
		if(b == false) break;
	}
	return 0;
}

bool filter::join_transcripts()  //TODO:make AS compatible
{
	sort(trs.begin(), trs.end(), transcript_cmp);
	//print();

	int32_t mind = min_bundle_gap;
	int ki = -1, kj = -1;
	for(int i = 0; i < trs.size(); i++)
	{
		int j = locate_next_transcript(i);
		if(j == -1) continue;
		if(trs[i].exons.size() >= 2 && trs[j].exons.size() >= 2) continue;
		int32_t d = trs[j].get_bounds().first.p32 - trs[i].get_bounds().second.p32;
		if(d > mind) continue;
		mind = d;
		ki = i;
		kj = j;
	}
	if(ki == -1 || kj == -1) return false;
	if(mind > min_bundle_gap - 1) return false;

	if(verbose >= 3) printf("join transcript %d and %d\n", ki, kj);

	if(trs[ki].exons.size() >= 2)
	{
		assert(trs[kj].exons.size() == 1);
		as_pos32 p1 = trs[ki].get_bounds().second;
		as_pos32 p2 = trs[kj].get_bounds().second;
		trs[ki].add_exon(p1, p2); // TODO:make AS some seq before p2 may have allele
		trs[kj].sort();
		trs[ki].shrink();
		trs.erase(trs.begin() + kj);
		return true;
	}
	else if(trs[kj].exons.size() >= 2)
	{
		assert(trs[ki].exons.size() == 1);
		as_pos32 p1 = trs[ki].get_bounds().first;
		as_pos32 p2 = trs[kj].get_bounds().first;
		trs[kj].add_exon(p1, p2); // TODO:make AS some seq before p2 may have allele
		trs[kj].sort();
		trs[kj].shrink();
		trs.erase(trs.begin() + ki);
		return true;
	}
	else
	{
		assert(trs[ki].exons.size() == 1);
		assert(trs[kj].exons.size() == 1);
		as_pos32 p1 = trs[ki].get_bounds().first;
		as_pos32 p2 = trs[kj].get_bounds().first;
		trs[kj].add_exon(p1, p2); // TODO:make AS some seq before p2 may have allele
		trs[kj].sort();
		trs[kj].shrink();
		double cov = 0;
		cov += trs[ki].coverage * trs[ki].length();
		cov += trs[kj].coverage * trs[kj].length();
		cov /= (trs[ki].length() + trs[kj].length());
		trs[kj].coverage = cov;
		trs.erase(trs.begin() + ki);
		return true;
	}

	return true;
}

int filter::locate_next_transcript(int t)
{
	if(t < 0 || t >= trs.size()) return -1;
	PI32 p = trs[t].get_bounds();
	int a = 0;
	int b = trs.size() - 1;
	if(trs[b].get_bounds().first.leftto(p.second)) return -1;
	while(true)
	{
		assert(a <= b);
		if(a == b) return a;
		int k = (a + b) / 2;
		if(trs[k].get_bounds().first.samepos(p.second)) return k;				// TODO: AS compatible
		if(trs[k].get_bounds().first.leftsameto(p.second)) a = k + 1;			// TODO: AS compatible
		if(trs[k].get_bounds().first.rightto(p.second)) b = k;					// TODO: AS compatible
	}
	assert(false);
	return -1;
}

int filter::merge_single_exon_transcripts(vector<transcript> &trs0)
{
	typedef pair<PI32, int> PPI;
	vector<PPI> vv;
	for(int i = 0; i < trs0.size(); i++)
	{
		vector<PI32> v = trs0[i].exons;
		for(int k = 0; k < v.size(); k++)
		{
			vv.push_back(PPI(v[k], i));
		}
	}

	sort(vv.begin(), vv.end());

	set<int> fb;
	for(int i = 0; i < vv.size(); i++)
	{
		as_pos32 p1 = vv[i].first.first;
		as_pos32 q1 = vv[i].first.second;
		int k1 = vv[i].second;
		transcript &t1 = trs0[k1];
		if(t1.exons.size() != 1) continue;
		if(t1.strand != '.') continue;

		bool b = false;
		for(int k = i - 1; k >= 0 && k >= i - 10; k--)
		{
			as_pos32 p2 = vv[k].first.first;
			as_pos32 q2 = vv[k].first.second;
			int k2 = vv[k].second;
			if(fb.find(k2) != fb.end()) continue;
			transcript &t2 = trs0[k2];
			if(t2.seqname != t1.seqname) continue;

			assert(p1.rightsameto(p2));
			if(q2.leftto(q1)) continue;

			//if(b == true) printf("AAA insert k1 = %d (%d, %d) to fb with k2 = %d (%d, %d)\n", k1, p1, q1, k2, p2, q2);

			b = true;
			break;
		}

		if(b == true) fb.insert(k1);
		if(b == true) continue;

		for(int k = i + 1; k < vv.size(); k++)
		{
			as_pos32 p2 = vv[k].first.first;
			as_pos32 q2 = vv[k].first.second;
			int k2 = vv[k].second;
			if(fb.find(k2) != fb.end()) continue;
			transcript &t2 = trs0[k2];
			if(t2.seqname != t1.seqname) continue;

			if(p2.rightto(p1)) break;
			assert(p2 == p1);
			if(q2.leftto(q1)) continue;
			b = true;

			//if(b == true) printf("BBB insert k1 = %d (%d, %d) to fb with k2 = %d (%d, %d)\n", k1, p1, q1, k2, p2, q2);

			break;
		}
		if(b == true) fb.insert(k1);
	}

	vector<transcript> v;
	for(int i = 0; i < trs0.size(); i++)
	{
		if(fb.find(i) != fb.end()) continue;
		v.push_back(trs0[i]);
	}
	trs0 = v;
	return 0;
}

int filter::merge_single_exon_transcripts()
{
	typedef vector<transcript> VT;
	typedef pair<string, VT> PSVT;
	typedef map<string, VT> MSVT;

	MSVT msvt;
	for(int i = 0; i < trs.size(); i++)
	{
		transcript &t = trs[i];
		string s = t.seqname;
		if(msvt.find(s) == msvt.end())
		{
			VT v;
			v.push_back(t);
			msvt.insert(PSVT(s, v));
		}
		else
		{
			msvt[s].push_back(t);
		}
	}

	trs.clear();
	for(MSVT::iterator it = msvt.begin(); it != msvt.end(); it++)
	{
		merge_single_exon_transcripts(it->second);
		trs.insert(trs.end(), it->second.begin(), it->second.end());
	}

	return 0;
}

int filter::print()
{
	for(int i = 0; i < trs.size(); i++)
	{
		transcript &t = trs[i];
		printf("transcript %d: exons = %lu, pos = %d-%d\n",
				i, t.exons.size(), t.get_bounds().first.p32, t.get_bounds().second.p32);
	}
	return 0;
}

bool transcript_cmp(const transcript &x, const transcript &y)
{
	if(x.exons[0].first < y.exons[0].first) return true;
	else return false;
}
