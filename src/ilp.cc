#include "ilp.h"
#include "config.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <numeric>
#include <limits>

#define PHASING_DIFF 0
#define PHASING_WEIGHT_RATIO_50 1
#define PHASING_COUNT_RATIO_50 1
#define RATIO_THRESHOLD 30

//int PHASING_COUNT;

using namespace std;

ilp::ilp(const splice_graph &_gr, const hyper_set &_HS, GRBEnv *env):gr(_gr),HS(_HS)
{
    //env = new GRBEnv();
    model = new GRBModel(*env);
    //model->getEnv().set(GRB_IntParam_OutputFlag, 0);
    model->set("TimeLimit", "5.0");
    model->set("MIPGap", "0.01");
    
    gr.get_edge_indices(i2e, e2i);
    HS.edges.clear();
    HS.ecnts.clear();
    HS.e2s.clear();

    HS.build(gr, e2i);
    H = HS.edges;
    H_weight = HS.ecnts;
    map_e_to_phasing = HS.e2s;

    cout << "Size of HS: " << HS.nodes.size() << endl;
    cout << "Size of _HS: " << _HS.nodes.size() << endl;
    cout << "Size of HS edges: " << HS.edges.size() << endl;
    cout << "Size of _HS edges: " << _HS.edges.size() << endl;

    cout << "Size of H: " << H.size() << endl;
    //cout << H_weight.size() << endl;
    Hv.clear();
    //cout << M << endl;
    int cnt = 0;
    double ave_pweight = 0;
    double mid_pweight = 0;
    double ave_plen = 0;
    for(vector < vector <int> >::iterator it = H.begin(); it != H.end(); it++,cnt++)
    {
        vector<int> hv;
        //cout << cnt << "\tWeight: " << H_weight[cnt] << endl;
        //print_phasing_path(*it);
        for(vector<int>::iterator itv = it->begin(); itv != it->end(); itv++)
        {
            int s = i2e[*itv]->source();
            int t = i2e[*itv]->target();
            if(itv == it->begin())
            {
                hv.push_back(s);
                //cout << s << '\t';
            }
            hv.push_back(t);
            //cout << t << '\t';
        }
        //cout << endl;
        Hv.push_back(hv);

        ave_pweight += H_weight[cnt];
        ave_plen += hv.size();
        if(cnt == H.size()/2) mid_pweight = H_weight[cnt];
    }
    ave_pweight /= 1.0*(H.size());
    ave_plen /= 1.0*(H.size());
    
    /*ofstream fout("phasing-path.txt");
    cnt = 0;
    for(int i = 0; i < Hv.size(); i++, cnt++)
    {
        fout << '>' << cnt << '-' << H_weight[cnt] << endl;
        for(vector<int>::const_iterator it = Hv[cnt].begin(); it != Hv[cnt].end(); it++)
        {
            fout << *it << ' ';
        }
        fout << endl;

    }
    fout.close();*/
   
    //topological_balance();
    //remove_nest_phasing_path();
    add_edge_not_in_H();
    //sqrtdiff_trivial_vertex();

    //POSSIBLE_INTRON.clear();
    //sqrtdiff_trivial_vertex();
    //remove_possible_intron();
    //

    cout << "Edges: " << endl;
    for(int i = 0; i < gr.num_edges(); i++)
    {
        edge_descriptor e = i2e[i];
        printf("[%d->%d]\tWeight: %.2f\n", e->source(), e->target(), gr.get_edge_weight(e));
    }

    M = H.size();
    cout << "\n#Phasing path: " << M << endl;
    for(int i = 0; i < M; i++)
    {
        printf("%d: weight = %d, list = ",i, H_weight[i]);
        print_phasing_path(H[i]);
		printf("\n");

    }
    printf("\n");

    assert(M == Hv.size());
    assert(M == H_weight.size());
    //gr.print();
    
    //N = gr.num_edges()-gr.num_vertices()+2;
    AntiChain.clear();
    //int anti_num = max_antichain();
    POSSIBLE_INTRON.clear();
    sqrtdiff_trivial_vertex();
    //remove_possible_intron();
    edge_start_mid_end();


    PHASING_PATH_THRESHOLD = 5.0;
    //PHASING_PATH_THRESHOLD = 0.01*gr.compute_average_edge_weight();
    
    if(M > RATIO_THRESHOLD)
    {
        PHASING_WEIGHT_RATIO = PHASING_WEIGHT_RATIO_50;
        PHASING_COUNT_RATIO = PHASING_COUNT_RATIO_50;
    }
    else
    {
        PHASING_WEIGHT_RATIO = 1.0;
        PHASING_COUNT_RATIO = 1.0;
    }

    cout << "PHASING_PATH_THRESHOLD: " << PHASING_PATH_THRESHOLD << endl;
    FLOW_NUM = gr.num_edges()-gr.num_vertices()+2;

    /*cout << "Sequence of N: " << endl;
    for(int i = 1; i <= M; i++) 
    {
        PHASING_COUNT = i;
        int anti_num = max_antichain();
        AntiChain.clear();

        //cout << anti_num << '\t';
        cout << "N: " << anti_num << endl;
        cout << endl;
    }
    cout << endl;*/

    N1 = max_antichain();
    cout << "\nAnti-Chain1:\n";
    set<int>::const_iterator itac;
    for(itac = AntiChain.begin(); itac != AntiChain.end(); itac++)
    {
        printf("%d: weight = %d, list = ",*itac, H_weight[*itac]);
        print_phasing_path(H[*itac]);
        printf("\n");
    }

    if(N1 == 1)
    {
        N2 = N1;
    }
    else
    {
        H_nest.clear();
        N2  = max_antichain();
    }
    
    
    cout << "\nAnti-Chain2:\n";
    for(itac = AntiChain.begin(); itac != AntiChain.end(); itac++)
    {
        printf("%d: weight = %d, list = ",*itac, H_weight[*itac]);
        print_phasing_path(H[*itac]);
        printf("\n");
    }
    printf("\n");

    /*max_weight_antichain();
    cout << "\nAnti-Chain-max-weight1:\n";
    for(itac = AntiChain.begin(); itac != AntiChain.end(); itac++)
    {
        printf("%d: weight = %d, list = ",*itac, H_weight[*itac]);
        print_phasing_path(H[*itac]);
        printf("\n");
    }*/

    /*AntiChain.clear();
    N1 = greedy_antichain();
    N2 = N1;*/

    printf("\nN1 = %d; N2 = %d\n", N1, N2);
    if(N1 > 0)
    {
        ofstream stat_file;
        stat_file.open("stat.csv", fstream::app);
        stat_file.setf(ios::fixed, ios::floatfield);
        stat_file.precision(2);
        stat_file << ave_pweight << ',' << mid_pweight << ',' << ave_plen << ',' <<  N1 << ',' << N2 << ',';
        stat_file.close();
    }

    /*while(anti_num - anti_num_ori >1 && anti_num < flow_num)
    {
        anti_num_ori = anti_num;
        anti_num = max_antichain();
    }*/
    //int anti_num_second = max_antichain();


    //int flow_num = gr.num_edges()-gr.num_vertices()+2;
    //cout << anti_num << '\t' << anti_num_second << endl;
    cout << "\nV-E+2 = " << FLOW_NUM << endl;
    /*if(anti_num >= flow_num)
        N = anti_num;
    else
    {
        N = min((anti_num+1)*anti_num/2, flow_num);
        //N = min(N, anti_num+5);
    }*/
    N = AntiChain.size();
    //assert(N == N1);

    MAX_EDGE_WEIGHT = gr.get_edge_weight(gr.compute_maximum_edge_w());
    //cout << H.size() << endl;
    cout << *max_element(H_weight.begin(), H_weight.end()) << endl;
    MAX_EDGE_WEIGHT = H.size()*max(MAX_EDGE_WEIGHT, 1.0*(*max_element(H_weight.begin(), H_weight.end())));
    cout << MAX_EDGE_WEIGHT << endl;
    assert(MAX_EDGE_WEIGHT > 1);
    //assert(N>0);
    //topological_balance();

}

ilp::~ilp()
{
	if(model != NULL) delete model;
}

int ilp::solve(int suggestN)
{
    cout << "N: " << N << endl;
    if(N<=0)
    {
        cout << "No path predicted!" << endl;
        return 0;
    }
    // variables
    if(verbose >= 2)printf("\nAdd variables for supported edge ...\n");
    add_edge_support_variables();
    
    if(verbose >= 2) printf("Add edge abundance variables ...\n");
    add_edge_abundance_variables();
    
    if(verbose >= 2) printf("Add predicted edge coverage variables for phasing path...\n");
    add_predicted_edge_coverage_variables();
    
    if(verbose >= 2)printf("Add predicted vertex coverage variables for phasing path...\n");
    add_predicted_vertex_coverage_variables();

    if(verbose >= 2) printf("Add variables for supported phasing path ...\n");
    add_phasing_support_variables();
    
    add_path_abundance_variables();
    add_predicted_path_coverage_variables(); 
    add_phasing_cover_variables();

    add_unbalanced_junction_variables();
    //add_duplicate_path_variables();
    add_non_phasing_edge_variables();

    //constraints
    if(verbose >= 2) printf("Add valid path constraint...\n");
    add_valid_path_constraint_st();
    
    if(verbose >= 2) printf("Add x to y constraint...\n");
    add_x_to_y_constraint();
    
    if(verbose >= 2) printf("Add phasing path cover constraint...\n");
    add_phasing_cover_constraint();

    printf("Add edge conflict constraint...\n");
    add_edge_conflict_constraint(); 

    if(verbose >= 2) printf("Add anti-chain constraint...\n");
    add_one_spot_constraint();
    
    if(verbose >= 2) printf("Add conflict&nest constraint...\n");
    add_conflict_nest_constraint();
    
    printf("Add t0x0 constraint...\n");
    add_t0x0_constraint();
    
    printf("Add x to t constraint...\n");
    add_x_to_t_constraint();
    
    add_s0y0_constraint();
    //add_y_to_s_constraint();
    add_s_to_t_constraint();

    add_unbalanced_junction_constraint();
    //add_duplicate_path_constraint();
    //add_adjacent_edge_constraint();

    if(verbose >= 2) printf("Set objective...\n");
    set_objective();
    
    if(verbose >= 2) printf("Optimizing...\n");
    model->update();
    model->write("debug.lp");
    
    int err = 0;

    try
    {
        model->optimize();
    }
    catch(GRBException e)
    {
        cout << "Error code =" << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
        err = 1;
    }
    catch(...)
    {
        cout << "Exception during optimization." << endl;
        err = 1;
    }

    double gap = -1;
    ofstream stat_file;
    stat_file.open("stat.csv", fstream::app);
    stat_file.setf(ios::fixed, ios::floatfield);
    stat_file.precision(2);
    if(!err)
    {
        gap = fabs(model->get(GRB_DoubleAttr_ObjBound)-model->get(GRB_DoubleAttr_ObjVal))/fabs(model->get(GRB_DoubleAttr_ObjVal));
        cout << "Final gap value: " << gap << endl;
        if(gap>0.2)
        {
            cout << "Gap of solution is large." << endl;
            err = 1;
            
            stat_file << -1 << ',' << -1  << ',' << -1 <<  ',' << -1 << ',' << -1  << ',' << -1 << ',';
        }
        else
        {
            if(verbose >= 2) cout << "Obj: " << model->get(GRB_DoubleAttr_ObjVal) << endl;
            collect_results(suggestN);
            
            stat_file << gap << ',' << obj_UBJ.get(GRB_DoubleAttr_X) << ',' << obj_phasing_count.get(GRB_DoubleAttr_X) << ',' << obj_edge_weight.get(GRB_DoubleAttr_X) << ',' << obj_phasing_cover.get(GRB_DoubleAttr_X) << ',' <<  model->get(GRB_DoubleAttr_ObjVal) << ',';

        }
    }
    else
    {
        stat_file << -1 << ',' << -1  << ',' << -1 << ',' << -1 << ',' << -1  << ',' << -1 << ',';
    }
   
    stat_file.close();
    return err;
}

int ilp::edge_start_mid_end()
{
    set<int> start_eset;
    set<int> end_eset;
    set<int> mid_eset;
    set<int> mid_end;
    set<int> mid_start;
    for(int i = 0; i < H.size(); i++)
    {
        int start = H[i][0];
        start_eset.insert(start);
        //cout << "Start: " << start << endl;
        int end = H[i][H[i].size()-1];
        end_eset.insert(end);
        //cout << "End: " << end << endl;
        for(int j = 1; j < H[i].size()-1; j++)
        {
            int mid = H[i][j];
            mid_eset.insert(mid);
            //cout << "Mid: " << mid << endl;
            //
            mid_end.insert(mid);
            mid_start.insert(mid);
        }
        mid_end.insert(end);
        mid_start.insert(start);
    }

    set<int> non_ovlp_start;
    set<int> non_ovlp_end;
    set_difference(start_eset.begin(), start_eset.end(), mid_eset.begin(), mid_eset.end(), inserter(non_ovlp_start, non_ovlp_start.begin()));
    set_difference(end_eset.begin(), end_eset.end(), mid_eset.begin(), mid_eset.end(), inserter(non_ovlp_end, non_ovlp_end.begin()));
    set<int>::const_iterator it1, it2;
    //cout << "Non_ovlp_start: " << endl;
    for(it1 = non_ovlp_start.begin(); it1 != non_ovlp_start.end(); it1++)
    {
        //cout << *it1 << endl;
        edge_descriptor e = i2e[*it1];
        int s = e->source();
        if(gr.get_out_weights(s)/gr.get_edge_weight(e) > 50)
        {
            NON_OVLP_EDGE.insert(*it1);
        }
    }

    //cout << "Non_ovlp_end: " << endl;
    for(it2 = non_ovlp_end.begin(); it2 != non_ovlp_end.end(); it2++)
    {
        //cout << *it2 << endl;
        edge_descriptor e = i2e[*it2];
        int t = e->target();
        if(gr.get_in_weights(t)/gr.get_edge_weight(e) > 50)
        {
            NON_OVLP_EDGE.insert(*it2);
        }
    }

    return 0;
}


int ilp::sqrtdiff_trivial_vertex()
{
    cout << "sqrtdiff: " << endl;
    SQRTDIFF_E.clear();
    SQRTDIFF.clear();
    for(int i = 0; i < gr.num_edges(); i++)
    {
        SQRTDIFF_E.push_back(1.0);
    }
    for(int i = 0; i < gr.num_vertices(); i++)
    {
        SQRTDIFF.push_back(1.0);
    }
    for(int v = 0; v < gr.num_vertices(); v++)
    {
        if(gr.in_degree(v) == 1 && gr.out_degree(v) == 1)
        {
            double diff = 1.0;

            PEEI pei1 = gr.in_edges(v);
            int prev = (* pei1.first)->source();
            int eid1 = e2i[*pei1.first];
            double w1 = gr.get_vertex_weight(prev);

            double w2 = gr.get_vertex_weight(v);

            PEEI pei2 = gr.out_edges(v);
            int post = (* pei2.first)->target();
            int eid2 = e2i[*pei2.first];
            double w3 = gr.get_vertex_weight(post);

            
            /*if(prev == 0) w1 = w3; //s
            if(post == gr.num_vertices()-1) w3 = w1; //t
            //cout << w1 << "\t" << w2 << "\t" << w3 << endl;
            diff = (w1+w3)/(2.0*w2);
            //cout << v << ": " << diff << endl;*/
            
            int prev_right = gr.get_vertex_info(prev).rpos;
            int post_left = gr.get_vertex_info(post).lpos;
            int v_right = gr.get_vertex_info(v).rpos;
            int v_left = gr.get_vertex_info(v).lpos;
            bool prev_close = (abs(prev_right-v_left)<3);
            bool post_close = (abs(post_left-v_right)<3);
            if(!prev_close && !post_close) continue;

            double intron_thresh = 20.0;
            PEB peb = gr.edge(prev, post);
            if(peb.second)
            {
                diff = max(diff, gr.get_edge_weight(peb.first)/w2);
                if(diff>=intron_thresh) POSSIBLE_INTRON.insert(v);
            }
            
            /*if(prev == 0 && gr.in_degree(post)>1)
            {
                w1 = w3;
                //diff = max(1.0, max(diff, (w1+w3)/(2.0*w2))/10);
                diff = max(diff, w3/w2);
                if(gr.degree(prev) == 1) diff = 1;
                if(diff>=intron_thresh) POSSIBLE_INTRON.insert(v);
             }
            
            if(post == (gr.num_vertices()-1) && gr.out_degree(prev)>1) 
            {                    
                w3 = w1;
                //diff = max(1.0, max(diff, (w1+w3)/(2.0*w2))/10);
                diff = max(diff, w1/w2);
                if(gr.degree(post) == 1) diff = 1;
                if(diff>=intron_thresh) POSSIBLE_INTRON.insert(v);
            }*/
            
            if(gr.in_degree(post)==1 && gr.out_degree(post)==1)
            {
                PEEI pei3 = gr.out_edges(post);
                int post_post = (* pei3.first)->target();
                int eid3 = e2i[*pei3.first];
                double w3_post = gr.get_vertex_weight(post_post);

                PEB peb2 = gr.edge(prev, post_post);
                if(peb2.second)
                {
                    diff = max(diff, gr.get_edge_weight(peb2.first)/max(w2,w3));
                    if(diff>=intron_thresh)
                    {
                        POSSIBLE_INTRON.insert(v);
                        POSSIBLE_INTRON.insert(post);
                    }
                }
                
                /*if(prev == 0 && gr.in_degree(post_post)>1)
                {
                    diff = max(diff, w3_post/max(w2,w3));
                    if(gr.degree(prev) == 1) diff = 1;
                    if(diff>=intron_thresh)
                    {
                        POSSIBLE_INTRON.insert(v);
                        POSSIBLE_INTRON.insert(post);
                    }

                }*/

            }
            
            /*if(post == (gr.num_vertices()-1) && gr.in_degree(prev)==1 && gr.out_degree(prev)==1)
            {
                PEEI pei3 = gr.in_edges(prev);
                int prev_prev = (* pei3.first)->source();
                int eid3 = e2i[*pei3.first];
                double w1_prev = gr.get_vertex_weight(prev_prev);
                if(gr.out_degree(prev_prev)>1)
                {
                    diff = max(diff, w1_prev/max(w1, w2));
                    if(gr.degree(post) == 1) diff = 1;
                    if(diff>=intron_thresh)
                    {
                        POSSIBLE_INTRON.insert(v);
                        POSSIBLE_INTRON.insert(prev);
                    }

                }

            }*/

            if(diff > 1)
            {
                cout << v << ": " << sqrt(diff) << endl;
                SQRTDIFF[v] = sqrt(diff);
                //SQRTDIFF.push_back(sqrt(diff));
                //cout << v << ": " << diff << endl;
                //SQRTDIFF.push_back(diff); 
            }
        }

    }
    for(int v = 0; v < gr.num_vertices(); v++)
    {
        if(SQRTDIFF[v]>1)
        {
            PEEI pei1 = gr.in_edges(v);
            int eid1 = e2i[*pei1.first];
            int prev = (* pei1.first)->source();

            PEEI pei2 = gr.out_edges(v);
            int eid2 = e2i[*pei2.first];
            int post = (* pei2.first)->target();

            SQRTDIFF_E[eid1] *= SQRTDIFF[v];
            SQRTDIFF_E[eid2] *= SQRTDIFF[v];

            printf("[%d->%d->%d] sqrtdiff:%f\n", prev, v, post, SQRTDIFF[v]);
        }
    }

    return 0;
}

int ilp::add_edge_not_in_H()
{
    set<int> phasing_edges;
    int n = H.size();
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < H[i].size(); j++)
        {
            phasing_edges.insert(H[i][j]);
        }
    }
    if(phasing_edges.size() == gr.num_edges()) return 0;

    vector<int> edge_not_phasing;
    set<int> all_edges;
    for(int i = 0; i < gr.num_edges(); i++)
    {
        all_edges.insert(i);
    }
    set_difference(all_edges.begin(), all_edges.end(), phasing_edges.begin(), phasing_edges.end(), inserter(edge_not_phasing, edge_not_phasing.begin()));
    for(int i = 0; i < edge_not_phasing.size(); i++)
    {
        int eid = edge_not_phasing[i];
        //cout << "Edge not phased: " << eid << '\t';
        edge_descriptor e = i2e[eid];
        int pcnt = (int)gr.get_edge_weight(e);
        //if(pcnt < 2) continue;

        int s = e->source(), t = e->target();
        //printf("[%d -> %d]\n", s, t);
        if(s == 0 || t == gr.num_vertices()-1) continue;
        cout << "Edge not phased: " << eid << '\t';
        printf("[%d -> %d]\n", s, t);

        int pid = H.size();
        vector<int> h_new; 
        h_new.push_back(eid);
        H.push_back(h_new);

        vector<int> hv_new;
        hv_new.push_back(s);
        hv_new.push_back(t);
        Hv.push_back(hv_new);

        H_weight.push_back(pcnt);
        //cout << H_weight.back() << endl;

        set<int> p_new;
        p_new.insert(pid);
        map_e_to_phasing.insert(make_pair(eid, p_new));

    }
    return 0;
}

int ilp::remove_possible_intron()
{
    vector < vector <int> > Hv_ori = Hv;
    vector < vector <int> > H_ori = H;
    vector <int> H_weight_ori = H_weight;
    Hv.clear();
    H.clear();
    H_weight.clear();
   
    for(int i = 0; i < Hv_ori.size(); i++)
    {
        bool clear = false;
        for(set<int>::const_iterator itv=POSSIBLE_INTRON.begin(); itv != POSSIBLE_INTRON.end(); itv++)
        {
            //cout << "lala" << *itv << endl;
            if(find(Hv_ori[i].begin(), Hv_ori[i].end(), *itv) != Hv_ori[i].end())
            {
                clear = true;
                cout << "Delete " << i << "th path." << endl;
                break;
            }
        }

        if(!clear)
        {
            Hv.push_back(Hv_ori[i]);
            H.push_back(H_ori[i]);
            H_weight.push_back(H_weight_ori[i]);
        }

    }
    return 0;
}

int ilp::remove_nest_phasing_path()
{
    vector < vector <int> > Hv_ori = Hv;
    vector < vector <int> > H_ori = H;
    vector <int> H_weight_ori = H_weight;
    Hv.clear();
    H.clear();
    H_weight.clear();
   
    //cout << Hv_ori.size()<< endl;
    int cnt = 0;
    ofstream fout("phasing-path-wo-nest.txt");
    for(int i = 0; i < Hv_ori.size(); i++)
    {
        bool clear = false;
        for(int j = 0; j < Hv_ori.size(); j++)
        {
            //printf("phasing %d <-> %d\n", i, j);
            if(i==j) continue;
            
            int end1 = Hv_ori[i].size()-1;
            int end2 = Hv_ori[j].size()-1;
            if(Hv_ori[i][0]>=Hv_ori[j][0] && Hv_ori[i][end1]<=Hv_ori[j][end2])
            {
                if(whether_compatible(Hv_ori[j], Hv_ori[i]))
                {
                    //cout << "Clear: " << i << " by " << j << endl;
                    clear = true;
                    break;
                }
                
            }
        }
        if(!clear)
        {
            Hv.push_back(Hv_ori[i]);
            H.push_back(H_ori[i]);
            H_weight.push_back(H_weight_ori[i]);
            fout << '>' << cnt << '-' << H_weight[cnt] << endl;
            //fout << '>' << cnt << endl;
            for(vector<int>::const_iterator it = Hv[cnt].begin(); it != Hv[cnt].end(); it++)
            {
                fout << *it << ' ';
            }
            fout << endl;
            cnt++;
        }
    }
    fout.close();
    return 0;
}

int ilp::topological_balance()
{
    for(int v = 1; v < gr.num_vertices()-1; v++)
    {
        double inw = gr.get_in_weights(v);
        double outw = gr.get_out_weights(v);
        double ratio = inw/outw;
        PEEI pei = gr.out_edges(v);
        edge_iterator ei1;
        for(ei1 = pei.first; ei1 != pei.second; ei1++)
        {
            double out_ew = gr.get_edge_weight(*ei1);
            gr.set_edge_weight(*ei1, out_ew*ratio);
            printf("[%d->%d] %f -> %f\n", (*ei1)->source(), (*ei1)->target(), out_ew, out_ew*ratio);
        }
    }
    return 0;
}

bool ilp::sort_pair_decreasing(pair<int, int> p1, pair<int, int> p2)
{
    return (p1.first>p2.first);
}

int ilp::max_weight_antichain()
{
    printf("\nCalculating antichain with max weight:\n");

    set< pair<int, int> > cnt_vid;
    for(int i = 0; i < M; i++)
    {
        cnt_vid.insert(make_pair(H_weight[i], i));
    }
    
    set<int> all;
    for(int i = 0; i < Hv.size(); i++)
    {
        all.insert(i);
    }

    set< pair<int, int> >::reverse_iterator itp;
    set<int> minor_ori = minor_phasing_path;
    map<int, int> nest_ori = H_nest;
    int max_acc_cnt = 0;
    set<int> max_weight_chain;
    int sec_acc_cnt = 0;
    set<int> sec_weight_chain;
    for(itp = cnt_vid.rbegin(); itp != cnt_vid.rend(); itp++)
    {
        AntiChain.clear();
        int vid = itp->second;
        minor_phasing_path.clear();
        H_nest.clear();

        /*printf("\nCandidate path id: %d conflict with ", vid);
        for(int i = 0; i < H_conflict[vid].size(); i++)
        {
            printf("%d\t", H_conflict[vid][i]);
        }
        printf("\n");*/

        set_difference(all.begin(), all.end(), H_conflict[vid].begin(), H_conflict[vid].end(), inserter(minor_phasing_path, minor_phasing_path.begin()));
        minor_phasing_path.erase(vid);
        int tempN = max_antichain();
        
        int acc_cnt = 0;
        for(set<int>::const_iterator it = AntiChain.begin();it != AntiChain.end(); it++)
        {
            acc_cnt += H_weight[*it];
        
        }
        
        printf("Antichain for phasing path %d: N = %d, chain's weight = %d, path: ", vid, tempN, acc_cnt);
        print_phasing_path(H[vid]);
        printf("\n");

        /*if( tempN == N1 && acc_cnt>max_acc_cnt)
        {
            sec_acc_cnt = max_acc_cnt;
            sec_weight_chain = max_weight_chain;

            max_acc_cnt = acc_cnt;
            max_weight_chain = AntiChain;
        }*/
        if(tempN == N1)
        {
            max_acc_cnt = acc_cnt;
            max_weight_chain.insert(AntiChain.begin(), AntiChain.end());
        }
    }
    AntiChain = max_weight_chain;
    //AntiChain.insert(sec_weight_chain.begin(), sec_weight_chain.end());
    minor_phasing_path = minor_ori;
    H_nest = nest_ori;

    return AntiChain.size();
}
int ilp::greedy_antichain()
{
    set<int> vec;
    for(int i = 0; i < Hv.size(); i++)
    {
        vec.insert(i);
    }
    set<int> all = vec;
    set<int> greedy_chain;
    set<int> minor_pp = minor_phasing_path;
    set<int> minor_temp;
    minor_phasing_path.clear();
    int greedyN = 0;
    while(vec.size()>0)
    {
        //if(greedyN == N1) break;
        
        int cnt = 0;
        int cand = 0;
        for(set<int>::const_iterator it = vec.begin(); it != vec.end(); it++)
        {
            if(H_weight[*it] > cnt)
            {
                cand = *it;
                cnt = H_weight[*it];
            }
        }

        //check whether cand is in max-antichain
        // vec.erase(cand);
        minor_temp = minor_phasing_path;
        set_difference(all.begin(), all.end(), H_conflict[cand].begin(), H_conflict[cand].end(), inserter(minor_phasing_path, minor_phasing_path.begin()));

        int tempN = max_antichain();
        printf("Check phasing path (%d): N = %d\n", cand, tempN);
        if(tempN == N1-greedyN-1)
        {
            greedyN++;
            greedy_chain.insert(cand);
            printf("Candidate path id: %d\n", cand);
            /*printf("Candidate path id: %d conflict with ", cand);
            for(int i = 0; i < H_conflict[cand].size(); i++)
            {
                printf("%d\t", H_conflict[cand][i]);
            }
            printf("\n");*/
            
            set<int> conflict(H_conflict[cand].begin(), H_conflict[cand].end());
            set<int> new_vec;
            set_intersection(vec.begin(), vec.end(), conflict.begin(), conflict.end(), inserter(new_vec, new_vec.begin()));
            vec = new_vec;
            
            set_difference(all.begin(), all.end(), vec.begin(), vec.end(), inserter(minor_phasing_path, minor_phasing_path.begin()));

        }
        else
        {
            minor_phasing_path = minor_temp;
            vec.erase(cand);
        }
        AntiChain.clear();
    }

    cout << "\nGreedy anti-chain:\n";
    for(set<int>::const_iterator it = greedy_chain.begin(); it != greedy_chain.end(); it++)
    {
        printf("%d: weight = %d, list = ",*it, H_weight[*it]);
        print_phasing_path(H[*it]);
        printf("\n");
    }

    assert(greedy_chain.size() == N1);
    AntiChain = greedy_chain;
    return AntiChain.size();
}

int ilp::max_antichain()
{
    bgr.clear();
    //cout << M << endl;
    for(int i = 0; i < 2*M; i++)
        bgr.add_vertex();
    //cout << M << endl;
    
    int s = bgr.num_vertices();
    bgr.add_vertex();
    int t = bgr.num_vertices();
    bgr.add_vertex();
   
    set< pair<int, int> > cnt_vid;
    set< pair<int, int> >::reverse_iterator itp;
    //set<int> minor_phasing_path;

    set<int>::const_iterator it_intron;
    /*cout << "Possible intron: ";
    for(it_intron = POSSIBLE_INTRON.begin(); it_intron != POSSIBLE_INTRON.end(); it_intron++)
    {
        cout << *it_intron << endl;
    }

    cout << endl;*/
    //cout << "Minor phasing path: ";
    //int cnt_all = 0;
    INTRON_PATH = 0;
    INTRON_PHASING_PATH.clear();
    double ignore_weight = 0.0;
    for(int i = 0; i < M; i++)
    {
        cnt_vid.insert(make_pair(H_weight[i], i));
        //cnt_all += H_weight[i];
        //print_phasing_path(H[i]);
        //if(H_weight[i]>=2*PHASING_PATH_THRESHOLD)
            //continue;

        /*int p_end = gr.get_vertex_info(Hv[i][0]).lpos;
        int junction = 0;
        for(vector<int>::const_iterator it = Hv[i].begin(); it != Hv[i].end(); it++)
        {
            int v_left = gr.get_vertex_info(*it).lpos;
            int v_right = gr.get_vertex_info(*it).rpos;
            //printf("p_end: %d\tv_left: %d", p_end, v_left);
            if(p_end != v_left )
            {
                junction = 1;
                break;
            }
            p_end = v_right;
        }
        if(!junction)
        {   
            ignore_weight += H_weight[i];
            minor_phasing_path.insert(i);
        }*/

        /*for(set<int>::const_iterator itv=POSSIBLE_INTRON.begin(); itv != POSSIBLE_INTRON.end(); itv++)
        {
            //cout << "lala" << *itv << endl;
            if(find(Hv[i].begin(), Hv[i].end(), *itv) != Hv[i].end())
            {
                minor_phasing_path.insert(i);
                cout << i << '\t';
                INTRON_PATH += 1;
                INTRON_PHASING_PATH.insert(i);
                break;
            }
        }*/
        
        /*for(set<int>::const_iterator ite = NON_OVLP_EDGE.begin(); ite != NON_OVLP_EDGE.end(); ite++)
        {
            //cout << "lala" << *itv << endl;
            if(find(H[i].begin(), H[i].end(), *ite) != H[i].end())
            {
                minor_phasing_path.insert(i);
                cout << "lala" << *ite << '\t';
                cout << i << '\t';
                INTRON_PATH += 1;
                INTRON_PHASING_PATH.insert(i);
                break;
            }
        }*/


    }
    //cout << endl;
    //cout << "INTRON_PATH: " << INTRON_PATH << endl;

    int cnt_acc = 0;
    int cnt = 0;
    
    //secondary antichain
    //cout << "Minor phasing path: ";
    for(set<int>::const_iterator it = AntiChain.begin(); it != AntiChain.end(); it++)
    {
        minor_phasing_path.insert(*it);
        //cout << *it << '\t';
    }
    
    int total_phasing_count = accumulate(H_weight.begin(), H_weight.end(), 0);
    //cout << "Total phasing count: " << total_phasing_count << endl;
    //cout << "PHASING_WEIGHT_RATIO: " << PHASING_WEIGHT_RATIO << endl;
    //cout << "PHASING_COUNT_RATIO:" << PHASING_COUNT_RATIO << endl;
    
    if(AntiChain.empty())
    {
        for(itp = cnt_vid.rbegin(); itp != cnt_vid.rend(); itp++)
        {
            int vid = itp->second;
            if(minor_phasing_path.find(vid)!=minor_phasing_path.end())
                continue;

            cnt_acc += itp->first;
            cnt++;
            if(cnt_acc > PHASING_WEIGHT_RATIO*total_phasing_count && cnt > PHASING_COUNT_RATIO*M)
            {
                minor_phasing_path.insert(vid);
            }
        }
    }
    
    PHASING_COVER_RATIO = 1.0;
    if(AntiChain.empty())
        MINOR_PATH = minor_phasing_path.size();
    //cout << "MINOR_PATH: " << MINOR_PATH << "\n\n";

    /*printf("\nPhasing path for anti-chain:\n");
    for(int i = 0; i < M; i++)
    {
        if(minor_phasing_path.find(i)!=minor_phasing_path.end())
            continue;
        printf("%d: weight = %d, list = ",i, H_weight[i]);
        print_phasing_path(H[i]);
		printf("\n");

    }*/
    //printf("\n");

    set<int> left, right;
    for(int i = 0; i < M; i++)
    {
         bgr.add_edge(s,i);
         left.insert(i);
    }
    for(int j = M; j < 2*M; j++)
    {
        bgr.add_edge(j,t);
        right.insert(j);
    }

    for(int i = 0; i < M; i++)
    {
        vector<int> vec_conflict;
        if(minor_phasing_path.find(i)!=minor_phasing_path.end()) 
        {
            if(H_conflict.size()<M) 
                H_conflict.push_back(vec_conflict);
            continue;
        }
        for(int j = M; j < 2*M; j++)
        {
            //printf("phasing %d <-> %d\n", i, j-M);
            if(i==j-M) continue;
            if(minor_phasing_path.find(j-M)!=minor_phasing_path.end()) continue;

            int end1 = Hv[i].size()-1;
            int end2 = Hv[j-M].size()-1;
            if(Hv[i][0]>=Hv[j-M][0] && Hv[i][end1]<=Hv[j-M][end2])
            {
                if(whether_compatible(Hv[j-M], Hv[i]))
                {
                    H_nest.insert(make_pair(i,j-M));
                    //cout << "Clear: " << i << " by " << j-M << endl;
                    //break;
                    continue;
                }
                else
                {
                    vec_conflict.push_back(j-M);
                    continue;
                }
            }
            if(Hv[i][end1]>Hv[j-M][end2])
            {
                if(Hv[i][0]<=Hv[j-M][0])
                {
                    if(!whether_compatible(Hv[i], Hv[j-M]))
                        vec_conflict.push_back(j-M);
                }
                else
                {
                    if(!whether_compatible(Hv[j-M], Hv[i]))
                        vec_conflict.push_back(j-M);
                }
                continue;
            }

            if(whether_compatible(Hv[i], Hv[j-M]))
            {
                //if(Hv[i][end1] < Hv[j-M][0]) continue;//no check_path
                bgr.add_edge(i,j);
                //printf("Add edge %d -> %d\n", i, j-M);
            }
            else
            {
                vec_conflict.push_back(j-M);
            }
        }
        if(H_conflict.size()<M) 
            H_conflict.push_back(vec_conflict);
    }
    assert(H_conflict.size() == M);
    
    for(map<int, int>::const_iterator ithn = H_nest.begin(); ithn != H_nest.end(); ithn++)
    {
        left.erase(ithn->first);
        right.erase(ithn->first+M);
        bgr.clear_vertex(ithn->first);
        bgr.clear_vertex(ithn->first+M);
    }
    
    for(set<int>::const_iterator it = minor_phasing_path.begin(); it != minor_phasing_path.end(); it++)
    {
        left.erase(*it);
        right.erase(*it);
        bgr.clear_vertex(*it);
        bgr.clear_vertex((*it)+M);
    }

    //draw_bipartite_graph("bgr-before.tex");
    
    vector<int> p;
    int loop_num = 0;
    while(bgr.compute_shortest_path(s,t,p))
    {
        assert(p.size()>=2);
        for(int i=0; i<p.size()-1; i++)
        {
            //printf("Rotate %d -> %d\n", p[i],p[i+1]);
            assert(bgr.check_path(p[i],p[i+1]));
            bgr.remove_edge(p[i],p[i+1]);
            bgr.add_edge(p[i+1],p[i]);
        }
        //cout << endl;
        loop_num++;
    }
    assert(loop_num==bgr.in_degree(s));
    assert(loop_num==bgr.out_degree(t));
    //draw_bipartite_graph("bgr-after.tex");
    
    //bgr.print();
    /*printf("Matching:\n");
    map<int, int> matching;
    for(int i = 0; i < M; i++)
    {
        if(bgr.in_degree(i)>0)
        {
            assert(bgr.in_degree(i)==1);
            PEEI pei = bgr.in_edges(i);
            edge_descriptor e = *(pei.first);
            if(e->source() == s) continue;
            printf("%d -- %d(%d)\n", e->target(), e->source()-M, e->source());
            matching.insert(make_pair(e->target(), e->source()-M));
        }
    }
    
    int key, new_key;
    vector<int> chain;
    if(!matching.empty())
    {
        key = matching.begin()->first;
    }
    while(!matching.empty())
    {
        chain.push_back(key);
        new_key = matching[key];
        matching.erase(key);
        key = new_key;
        if(matching.find(key) == matching.end())
        {
            chain.push_back(key);
            cout << endl;
            for(vector<int>::const_iterator itm = chain.begin(); itm != chain.end(); itm++)
            {
                print_phasing_path(H[*itm]);
            }
            chain.clear();
            if(!matching.empty())
            {
                key = matching.begin()->first;
            }
        }
    }*/

    vector<int> cc;
    bgr.bfs(s,cc);
    set<int> cc_set(cc.begin(), cc.end());
    
    set<int>::const_iterator itcc;
    /*cout << "Connected Components: " << endl;
    for(itcc = cc_set.begin(); itcc != cc_set.end(); itcc++)
    {
        cout << *itcc << '\t';
    }
    cout << endl;*/
    
    set<int> vertex_cover;
    set_difference(left.begin(), left.end(), cc_set.begin(), cc_set.end(), inserter(vertex_cover,vertex_cover.end()));
    set_intersection(right.begin(), right.end(), cc_set.begin(), cc_set.end(), inserter(vertex_cover,vertex_cover.end()));
    
    set<int>::const_iterator itvc;
    set<int> left_vc;
    //cout << "Vertex Cover: " << endl;
    for(itvc = vertex_cover.begin(); itvc != vertex_cover.end(); itvc++)
    {
        //cout << *itvc << '\t';
        if(*itvc >= M)
        {
            left_vc.insert(*itvc-M);
        }
        else
        {
            left_vc.insert(*itvc);
        }
    }
    //cout << endl;
    
    /*set<int> right_to_left;
    for(set<int>::const_iterator itrl=right.begin(); itrl != right.end(); itrl++)
    {
        right_to_left.insert(*itrl-M);
    }*/
    
    //cout << "Anti-Chain:\n";
    //AntiChain.clear();
    set_difference(left.begin(), left.end(), left_vc.begin(), left_vc.end(), inserter(AntiChain,AntiChain.end()));

    /*set<int>::const_iterator itac;
    for(itac = AntiChain.begin(); itac != AntiChain.end(); itac++)
    {
        printf("%d: weight = %d, list = ",*itac, H_weight[*itac]);
        print_phasing_path(H[*itac]);
        printf("\n");
    }*/

    return AntiChain.size();
}

int ilp::draw_bipartite_graph(const string &file)
{
    MIS mis;
    char buf[10240];
    for(int i = 0; i < bgr.num_vertices(); i++)
    {
        sprintf(buf, "%d", i);
        mis.insert(PIS(i, buf));
    }
    
    MES mes;
    edge_iterator it;
    for(it=(bgr.edges()).first; it != (bgr.edges()).second; it++)
    {
        sprintf(buf, " ");
        mes.insert(PES(*it, buf));
    }
    
    bgr.draw(file, mis, mes, 4.5);
    return 0;
}

bool ilp::whether_compatible(vector<int> &h1, vector<int> &h2)//h1<h2
{
    int end1 = h1.size()-1;
    int end2 = h2.size()-1;
    if(h1[end1] <= h2[0])
    {
        if(gr.check_path(h1[end1], h2[0]))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    int cmp = -1;
    for(int i = 0; i <= end1; i++)
    {
        if(h1[i]<h2[0])
        {
            continue;
        }
        else if(h1[i]>h2[0])
        {
            return false;
        }
        else
        {
            assert(h1[i]==h2[0]);
            int j = 0;
            while(i<=end1 && j<=end2)
            {
                if(h1[i]==h2[j])
                {
                    i++;
                    j++;
                }
                else
                    return false;
            }
            break;
        }
    }
    return true;
}

int ilp::print_phasing_path(const vector<int> &phasing_edge)
{
    set<int> vpath;
    vector<int>::const_iterator it;
    //cout << phasing_edge.size() << endl;
    for(it = phasing_edge.begin(); it != phasing_edge.end(); it++)
    {
        assert(*it < i2e.size());
        edge_descriptor e = i2e[*it];
        int s = e->source();
        int t = e->target();
        //cout << s << '\t' << t << endl;
        vpath.insert(s);
        vpath.insert(t);
    }
    //cout << "set: " << vpath.size() << endl;
    set<int>::const_iterator vit;

	assert(phasing_edge.size() + 1 == vpath.size());

    for(vit=vpath.begin(); vit != vpath.end(); vit++)
    {
        cout << *vit << ' ';
    }
    //printf("\n");
    return 0;
}


int ilp::add_edge_support_variables()
{
    X.clear();
    char name[1024];
    for(int i = 0; i < N; i++)
    {
        vector<GRBVar> Xi;
        for(int j = 0; j < gr.num_edges(); j++)
        {
            sprintf(name, "x-%d-%d", i,j);
            GRBVar var = model->addVar(0, 1, 0, GRB_BINARY, name);
            Xi.push_back(var);
        }
        assert(Xi.size()==gr.num_edges());
        X.push_back(Xi);
    }
    assert(X.size()==N);
    model->update();

    printf("Total %d*%d variables for supported edge\n", N, (int)gr.num_edges());
    return 0;
}

int ilp::add_edge_abundance_variables()
{
    T.clear();
    char name[1024];
    for(int i = 0; i < N; i++)
    {
        vector<GRBVar> Ti;
        for(int j = 0; j < gr.num_edges(); j++)
        {
            sprintf(name, "t-%d-%d", i,j);
            GRBVar var = model->addVar(0, MAX_EDGE_WEIGHT, 0, GRB_CONTINUOUS, name);
            //GRBVar var = model->addVar(0, gr.get_edge_weight(i2e[j]), 0, GRB_CONTINUOUS, name);
            Ti.push_back(var);
        }
        assert(Ti.size()==gr.num_edges());
        T.push_back(Ti);
    }
    assert(T.size()==N);
    model->update();
    
    printf("Total %d*%d edge abundance variables\n", N, (int)gr.num_edges());
    return 0;
}


int ilp::add_predicted_edge_coverage_variables()
{
    EDGE_SUPPORT.clear();
    char name[1024];
    for(int j = 0; j < gr.num_edges(); j++)
    {
        sprintf(name, "predicted-%d", j);
        double ew = gr.get_edge_weight(i2e[j]);

        double min_weight = 0;
        int s = i2e[j]->source();
        int t = i2e[j]->target();
        if(POSSIBLE_INTRON.find(s) == POSSIBLE_INTRON.end() && POSSIBLE_INTRON.find(t) == POSSIBLE_INTRON.end() && NON_OVLP_EDGE.find(j) == NON_OVLP_EDGE.end())
        {
            //min_weight = max(0.0, 0.5*(ew-PHASING_PATH_THRESHOLD));
            //cout << s << '\t' << t << "\tWeight: " << ew << endl;
        }

        //double min_weight = max(0.0, 0.5*(ew-PHASING_PATH_THRESHOLD));
        //double min_weight = 0;
        GRBVar var = model->addVar(min_weight, GRB_INFINITY , 0, GRB_CONTINUOUS, name);
        //model->addConstr(var, GRB_GREATER_EQUAL, 0.5*ew-1);
        EDGE_SUPPORT.push_back(var);
    }
    model->update();
    printf("Total %d predicted coverage variables.\n", (int)gr.num_edges());
    
    return 0;
}

int ilp::add_predicted_vertex_coverage_variables()
{
    VERTEX_SUPPORT.clear();
    char name[1024];
    for(int j = 0; j < gr.num_vertices(); j++)
    {
        sprintf(name, "v-abd-%d", j);
        //double min_weight = max(0.0, 0.5*(gr.get_vertex_weight(j)-PHASING_PATH_THRESHOLD));
        double min_weight = 0;
        GRBVar var = model->addVar(min_weight,GRB_INFINITY , 0, GRB_CONTINUOUS, name);
        VERTEX_SUPPORT.push_back(var);
    }
    model->update();
    printf("Total %d predicted vertex coverage variables.\n", (int)gr.num_vertices());

    return 0;
}

int ilp::add_phasing_support_variables()
{
    Y.clear();
    char name[1024];
    for(int i = 0; i < N; i++)
    {
        vector<GRBVar> Yi;
        for(int j = 0; j < M; j++)
        {
            sprintf(name, "y-%d-%d", i,j);
            GRBVar var = model->addVar(0, 1, 0, GRB_BINARY, name);
            Yi.push_back(var);
        }
        assert(Yi.size()==M);
        Y.push_back(Yi);
    }
    assert(Y.size()==N);
    model->update();
    
    if(verbose >= 2) printf("Total %d*%d variables for phasing path\n", N, M);
    return 0;
}


int ilp::add_path_abundance_variables()
{
    S.clear();
    char name[1024];
    for(int i = 0; i < N; i++)
    {
        vector<GRBVar> Si;
        for(int j = 0; j < M; j++)
        {
            sprintf(name, "s-%d-%d", i,j);
            GRBVar var = model->addVar(0, MAX_EDGE_WEIGHT, 0, GRB_CONTINUOUS, name);
            Si.push_back(var);
        }
        assert(Si.size()==M);
        S.push_back(Si);
    }
    assert(S.size()==N);
    model->update();

    printf("Total %d*%d path abundance variables\n", N, M);
    return 0;
}

int ilp::add_predicted_path_coverage_variables()
{
    Z.clear();
    char name[1024];
    for(int j = 0; j < M; j++)
    {
        sprintf(name, "z-%d", j);
        //cout << H_weight[j] << endl;
        //double min_count = max(0.5*(H_weight[j]-PHASING_PATH_THRESHOLD), H_weight[j]-2.0*PHASING_PATH_THRESHOLD);
        double min_count = 0;
        if(INTRON_PHASING_PATH.find(j) == INTRON_PHASING_PATH.end())
        {
            //min_count = max(0.0, (H_weight[j]-PHASING_PATH_THRESHOLD)*0.5);
            //min_count = max(0.0, sqrt(H_weight[j])-PHASING_PATH_THRESHOLD);

        }
        GRBVar var = model->addVar(min_count, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
        //GRBVar var = model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
        Z.push_back(var);
    }
    model->update();
    printf("Total %d predicted coverage variables.\n", M);

    return 0;
}

int ilp::add_phasing_cover_variables()
{
    PHASING_SUPPORT.clear();
    char name[1024];
    for(int j = 0; j < M; j++)
    {
        sprintf(name, "phasing-support-%d", j);
        GRBVar var = model->addVar(0, 1, 0, GRB_BINARY, name);
        PHASING_SUPPORT.push_back(var);
    }
    model->update();
    printf("Total %d phasing cover variables.\n", M);

    return 0;

}

int ilp::add_unbalanced_junction_variables()
{
    UBJ.clear();
    for(int i = 0; i < N; i++)
    {
        vector< vector<GRBVar> > uj;
        for(int j = 0; j < gr.num_edges(); j++)
        {
            vector<GRBVar> uk;
            for(int k = 0; k < gr.num_edges(); k++)
            {
                GRBVar var = model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
                uk.push_back(var);
            }
            uj.push_back(uk);
        }
        UBJ.push_back(uj);
    }
    return 0;
}

int ilp::add_duplicate_path_variables()
{
    D.clear();
    for(int i = 0; i < N; i++)
    {
        vector<GRBVar> Di;
        for(int j = i+1; j < N; j++)
        {
            GRBVar var = model->addVar(0, 1, 0, GRB_BINARY);
            Di.push_back(var);
        }
        D.push_back(Di);
    }
    return 0;
}

int ilp::add_non_phasing_edge_variables()
{
    NPE.clear();
    char name[1024];
    for(int i = 0; i < N; i++)
    {
        vector<GRBVar> NPEi;
        for(int j = 0; j < gr.num_edges(); j++)
        {
            sprintf(name, "npe-%d-%d", i,j);
            GRBVar var = model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, name);
            NPEi.push_back(var);
        }
        assert(NPEi.size()==gr.num_edges());
        NPE.push_back(NPEi);
    }
    assert(NPE.size()==N);
    model->update();
    
    printf("Total %d*%d non-phasing edge variables\n", N, (int)gr.num_edges());
    return 0;
}


int ilp::add_valid_path_constraint()
{
    char name[1024];
    edge_iterator ei1, ei2;
    int cnt = 0;
    for(int i = 0; i < N; i++)
    {
        GRBLinExpr expr_end;
        for(int v = 0; v < gr.num_vertices(); v++)
        {
            //cout << "v: " << v << "\tWeight: " << gr.get_vertex_weight(v) << endl;
            
            GRBLinExpr expr_in;
            PEEI pei = gr.in_edges(v);
            for(ei1 = pei.first, ei2 = pei.second; ei1 != ei2; ei1++)
            {
                int j = e2i[*ei1];
                expr_in += X[i][j];
                //printf("X[%d][%d]\t",i,j);
            }
            //cout << endl;
            
            sprintf(name, "sp-in-%d-%d", i, v);//simple path
            model->addConstr(expr_in, GRB_LESS_EQUAL, 1, name);
            cnt++;
            
            GRBLinExpr expr_out;
            pei = gr.out_edges(v);
            for(ei1 = pei.first, ei2 = pei.second; ei1 != ei2; ei1++)
            {
                int j = e2i[*ei1];
                expr_out += X[i][j];
                //printf("X[%d][%d]\t",i,j);
            }
            //cout << endl;
            
            sprintf(name, "sp-out-%d-%d", i, v);//simple path
            model->addConstr(expr_in, GRB_LESS_EQUAL, 1, name);
            cnt++;
            
            GRBVar in_out_diff = model->addVar(-1, 1 , 0, GRB_CONTINUOUS);
            model->addConstr(in_out_diff, GRB_EQUAL, expr_in-expr_out);
            
            GRBVar in_out_abs = model->addVar(0, 1, 0, GRB_CONTINUOUS);
            model->addGenConstrAbs(in_out_abs,in_out_diff);
            expr_end += in_out_abs;
        }
        
        sprintf(name, "ve-in-%d", i);//valid end
        model->addConstr(expr_end, GRB_EQUAL, 2, name);
        cnt++;
    }
    
    printf("Total %d valid path constraints.\n", cnt);
    return 0;
}

int ilp::add_valid_path_constraint_st()
{
    char name[1024];
    edge_iterator ei1, ei2;
    PEEI pei;
    int cnt = 0;
    for(int i = 0; i < N; i++)
    {
        for(int v = 0; v < gr.num_vertices(); v++)
        {
            //cout << "v: " << v << "\tWeight: " << gr.get_vertex_weight(v) << endl;

            if(gr.degree(v)==0)
                continue;

            GRBLinExpr expr_in;
            pei = gr.in_edges(v);
            for(ei1 = pei.first, ei2 = pei.second; ei1 != ei2; ei1++)
            {
                int j = e2i[*ei1];
                expr_in += X[i][j];
                //printf("X[%d][%d]\t",i,j);
            }
            //cout << endl;
            if(gr.out_degree(v)==0)//t
            {
                sprintf(name, "sp-in-%d-%d", i, v);//simple path
                model->addConstr(expr_in, GRB_EQUAL, 1, name);
                cnt++;
                continue;
            }
            else
            {
                sprintf(name, "sp-in-%d-%d", i, v);//simple path
                model->addConstr(expr_in, GRB_LESS_EQUAL, 1, name);
                cnt++;
            }

            GRBLinExpr expr_out;//valid path support
            pei = gr.out_edges(v);
            for(ei1 = pei.first, ei2 = pei.second; ei1 != ei2; ei1++)
            {
                int j = e2i[*ei1];
                expr_out += X[i][j];
                //printf("X[%d][%d]\t",i,j);
            }
            //cout << endl;

            if(gr.in_degree(v)==0)
            {
                sprintf(name, "sp-out-%d-%d", i, v);//simple path
                model->addConstr(expr_out, GRB_EQUAL, 1, name);
                cnt++;
                continue;
            }
            else
            {
                sprintf(name, "sp-out-%d-%d", i, v);//simple path
                model->addConstr(expr_out, GRB_LESS_EQUAL, 1, name);
                cnt++;
            }

            sprintf(name, "balance-%d-%d", i, v);//valid path support
            model->addConstr(expr_in, GRB_EQUAL, expr_out, name);
            cnt++;
        }
    }

    printf("Total %d valid path constraints.\n", cnt);
    return 0;
}



int ilp::add_x_to_y_constraint()
{
    char name[1024];
    int cnt = 0;
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < M; j++)
        {
            GRBLinExpr expr;
            vector<int> phasing_edge = H[j];
            vector<int>::const_iterator it;
            //printf("Phasing path %d contains:\n", j);
            for(it = phasing_edge.begin(); it != phasing_edge.end(); it++)
            {
                expr += X[i][*it];
                //printf("X[%d][%d]\t", i, *it);
                
                sprintf(name, "xy1-%d-%d-%d", i,j,*it);//simple path
                model->addConstr(Y[i][j], GRB_LESS_EQUAL, X[i][*it], name);
                cnt++;
            }
            //cout << endl;
            
            expr = expr - phasing_edge.size() + 1;
            sprintf(name, "xy2-%d-%d", i,j);//simple path
            model->addConstr(Y[i][j], GRB_GREATER_EQUAL, expr, name);
            cnt++;
            
        }
    }
    
    printf("Total %d x to y constraints.\n", cnt);
    return 0;
}

int ilp::add_phasing_cover_constraint()
{
    char name[1024];
    int cnt = 0;
    
    GRBLinExpr expr_total;
    for(int j = 0; j < M; j++)
    {
        GRBLinExpr expr;
        for(int i = 0; i < N; i++)
        {
            expr += Y[i][j];
            model->addConstr(PHASING_SUPPORT[j], GRB_GREATER_EQUAL, Y[i][j]);
        }
        sprintf(name, "pc-%d",j);//phasing cover
        //model->addConstr(expr, GRB_GREATER_EQUAL, 1, name);
        model->addConstr(PHASING_SUPPORT[j], GRB_LESS_EQUAL, expr, name);

        //int pp_cnt = H_weight[j];
        //model->addConstr(expr, GRB_LESS_EQUAL, pp_cnt/(2*PHASING_PATH_THRESHOLD)+1);
        cnt++;
        expr_total += PHASING_SUPPORT[j];
        
        //if(minor_phasing_path.find(j) == minor_phasing_path.end())
            //model->addConstr(PHASING_SUPPORT[j], GRB_GREATER_EQUAL, 0.9);

    }

    //model->addConstr(expr_total, GRB_GREATER_EQUAL, PHASING_COVER_RATIO*(M-MINOR_PATH));
    //model->addConstr(expr_total, GRB_GREATER_EQUAL, 0.95*M);

    if(verbose >= 2) printf("Total %d phasing path cover constraints.\n", cnt);
    return 0;
}

int ilp::set_one_spot_phasing_path(int pred, int pid)
{
    for(int i = 0; i < H_conflict[pid].size(); i++)
    {
        int conf = H_conflict[pid][i];
        //cout << "conf" << conf << endl;
        model->addConstr(Y[pred][conf], GRB_EQUAL, 0);
    }

    for(int i = 0; i < H[pid].size(); i++)
    {
        int eid = H[pid][i];
        //cout << "eid: " << eid << endl;
        model->addConstr(X[pred][eid], GRB_EQUAL, 1);
        for(int j = 0; j < H_X_conflict[eid].size(); j++)
        {
            int econf = H_X_conflict[eid][j];
            //cout << "econf: " << econf << endl
            model->addConstr(X[pred][econf], GRB_EQUAL, 0);
        }
    }
    return 0;
}

int ilp::add_one_spot_constraint()
{
    char name[1024];
    int cnt = 0;

    /*set< pair<int,int> > cnt_pid;
    for(int i = 0; i < H_weight.size(); i++)
    {
        cnt_pid.insert(make_pair(H_weight[i], i));
    }
    //sort(cnt_pid.begin(), cnt_pid.end());

    set<int> one_spot = AntiChain;
    int dtm = 0;
    for(set<int>::const_iterator it = AntiChain.begin(); it != AntiChain.end(); it++, dtm++)
    {
        if(H_weight[*it] < 2)
        {
            set< pair<int, int> >::reverse_iterator itp;
            for(itp = cnt_pid.rbegin(); itp != cnt_pid.rend(); itp++)
            {
                int partner = itp->second;
                int new_cnt = H_weight[partner];
                if(Hv[*it][0] >= Hv[partner][0] && whether_compatible(Hv[partner], Hv[*it]) && new_cnt >= H_weight[*it])
                {
                    sprintf(name, "one-spot-%d-%d", *it, partner);
                    model->addConstr(Y[dtm][partner], GRB_EQUAL, 1, name);
                    cout << *it << " changed to " << partner << ":" << endl;
                    print_phasing_path(H[partner]);
                    cnt_pid.erase(*itp);
                    cnt_pid.insert(make_pair(new_cnt/2, partner));
                    set_one_spot_phasing_path(dtm, partner);
                    break;
                }
                if(Hv[*it][0] < Hv[partner][0] && whether_compatible(Hv[*it], Hv[partner]) && H_weight[partner] >= H_weight[*it])
                {
                    sprintf(name, "one-spot-%d-%d", *it, partner);
                    model->addConstr(Y[dtm][partner], GRB_EQUAL, 1, name);
                    cout << *it << " changed to " << partner << ":" << endl;
                    print_phasing_path(H[partner]);
                    cnt_pid.erase(*itp);
                    cnt_pid.insert(make_pair(new_cnt/2, partner));

                    set_one_spot_phasing_path(dtm, partner);
                   break;
                }

            }
        }
        else
        {
            sprintf(name, "one-spot-%d", *it);
            model->addConstr(Y[dtm][*it], GRB_EQUAL, 1, name);
            set_one_spot_phasing_path(dtm, *it);

        }
    }*/
    
    set< pair<int, int> > ac_cnt;
    for(set<int>::const_iterator it = AntiChain.begin(); it != AntiChain.end(); it++)
    {
        ac_cnt.insert(make_pair(H_weight[*it], *it));
    }
    
    vector<int> antichain_vec;
    set< pair<int, int> >::reverse_iterator itp;
    for(itp = ac_cnt.rbegin(); itp != ac_cnt.rend(); itp++)
    {
        //cout << itp->second << ": " << itp->first << endl;
        antichain_vec.push_back(itp->second);
    }

    int ac_num = antichain_vec.size();
    int ac = 0;
    for(int i = 0; i < N; i++, ac++)
    {
        if(ac >= ac_num)
        {
            ac -= ac_num;
            antichain_vec.pop_back();
            ac_num --;
            assert(ac_num);
        }
        int ac_pp = antichain_vec[ac];
        model->addConstr(Y[i][ac_pp], GRB_EQUAL, 1);
        set_one_spot_phasing_path(i, ac_pp);
        //cout << ac_pp << endl;
    }


    /*for(int k = 0; k < AntiChain.size(); k++)
    {
        int dtm = 0;

        if(H_weight[antichain_vec[k]] < 1) 
        {
            cout << "Low coverage phasing path: ";
            print_phasing_path(H[antichain_vec[k]]);
            continue;
        }

        for(set<int>::const_iterator it = AntiChain.begin(); it != AntiChain.end(); it++, dtm++)
        {
            if(dtm==k)
            {
                sprintf(name, "ac-%d-%d",k,*it);
                model->addConstr(Y[k][*it], GRB_EQUAL, 1, name);
                //cout << name << ": " << 1 << endl;
            }
            else
            {
                sprintf(name, "ac-%d-%d",k,*it);
                model->addConstr(Y[k][*it], GRB_EQUAL, 0, name);
                //cout << name << ": " << 0 << endl;
            }
            cnt++;
        }
    }*/
    printf("Total %d anti-chain constraints.\n", cnt);
    return 0;
}

int ilp::add_edge_conflict_constraint()
{
    for(int i = 0; i < gr.num_edges(); i++)
    {
        vector<int> hx;
        edge_descriptor e1 = i2e[i];
        int s1 = e1->source(), t1 = e1->target();
        for(int j = 0; j < gr.num_edges(); j++)
        {
            if(i == j) continue;
            edge_descriptor e2 = i2e[j];
            int s2 = e2->source(), t2 = e2->target();
            if(t2 <= s1 && gr.check_path(t2,s1)) continue;
            if(t1 <= s2 && gr.check_path(t1,s2)) continue;
            hx.push_back(j);
        }
        H_X_conflict.push_back(hx);
    }
    
    char name[1024];
    int cnt = 0;
    for(int i = 0; i < N; i++)
    {
        for(int j=0; j<H_X_conflict.size();j++)
        {
            if(H_X_conflict[j].size()==0)
            {
                //X[i][j] = 1;
                model->addConstr(X[i][j], GRB_EQUAL, 1);
                //printf("X[%d][%d]\n", i, j);
                continue;
            }
            
            sprintf(name, "edge-conflict-%d-%d",i,j);
            //cout << name << endl;
            
            GRBLinExpr expr;
            for(int k = 0; k<H_X_conflict[j].size(); k++)
            {
                expr += X[i][H_X_conflict[j][k]];
                //cout << H_conflict[j][k] << endl;
            }
            model->addConstr(expr, GRB_LESS_EQUAL, H_X_conflict[j].size()*(1-X[i][j]), name);
            cnt++;
        }
    }
    printf("Total %d edge conflict constraints.\n", cnt);
    
    return 0;
}

int ilp::add_conflict_nest_constraint()
{
    char name[1024];
    int cnt = 0;
    for(int i = 0; i < N; i++)
    {
        for(int j=0; j<H_conflict.size();j++)
        {
            if(H_conflict[j].size()==0)
            {
                continue;
            }
            
            sprintf(name, "conflict-%d-%d",i,j);
            //cout << name << endl;
            
            GRBLinExpr expr;
            for(int k = 0; k<H_conflict[j].size(); k++)
            {
                expr += Y[i][H_conflict[j][k]];
                //cout << H_conflict[j][k] << endl;
            }
            model->addConstr(expr, GRB_LESS_EQUAL, H_conflict[j].size()*(1-Y[i][j]), name);
            cnt++;
        }
        
        for(map<int, int>::const_iterator ithn = H_nest.begin(); ithn != H_nest.end(); ithn++)
        {
            sprintf(name, "nest-%d-%d-%d",i,ithn->first,ithn->second);
            model->addConstr(Y[i][ithn->first], GRB_GREATER_EQUAL, Y[i][ithn->second], name);
            //cout << name << endl;
            cnt++;
        }
    }
    printf("Total %d conflict&nest constraints.\n", cnt);
    
    return 0;
}

int ilp::add_t0x0_constraint()
{
    char name[1024];
    int cnt = 0;
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < gr.num_edges(); j++)
        {
            sprintf(name, "t0x0-%d-%d", i,j);//simple path
            //cout << name << endl;
            model->addConstr(T[i][j], GRB_LESS_EQUAL, MAX_EDGE_WEIGHT*X[i][j], name);
            cnt++;
        }
    }
    
    printf("Total %d t0x0 constraints.\n", cnt);
    return 0;
}

int ilp::add_x_to_t_constraint()
{
    char name[1024];
    int cnt = 0;
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < gr.num_edges(); j++)
        {
            for(int k = j+1; k < gr.num_edges(); k++)
            {
                GRBLinExpr expr;
                
                GRBVar tdiff = model->addVar(-MAX_EDGE_WEIGHT, MAX_EDGE_WEIGHT, 0, GRB_CONTINUOUS);
                model->addConstr(tdiff, GRB_EQUAL, T[i][j]-T[i][k]);
                
                GRBVar tabs = model->addVar(0, MAX_EDGE_WEIGHT, 0, GRB_CONTINUOUS);
                model->addGenConstrAbs(tabs,tdiff);
                
                expr = MAX_EDGE_WEIGHT*(2-X[i][j]-X[i][k]);
                
                sprintf(name, "xt-%d-%d-%d", i,j,k);//simple path
                model->addConstr(tabs, GRB_LESS_EQUAL, expr, name);
                cnt++;
            }
        }
    }
    
    printf("Total %d x to t constraints.\n", cnt);
    return 0;
}


int ilp::add_s0y0_constraint()
{
    char name[1024];
    int cnt = 0;
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < M; j++)
        {
            sprintf(name, "s0y0-%d-%d", i,j);//simple path
            model->addConstr(S[i][j], GRB_LESS_EQUAL, MAX_EDGE_WEIGHT*Y[i][j], name);
            cnt++;
            model->addConstr(S[i][j], GRB_GREATER_EQUAL, 0.5*Y[i][j]);
            //cnt++;
        }
    }

    printf("Total %d s0y0 constraints.\n", cnt);
    return 0;
}

int ilp::add_y_to_s_constraint()
{
    char name[1024];
    int cnt = 0;
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < M; j++)
        {
            for(int k = j+1; k < M; k++)
            {
                GRBLinExpr expr;
                
                GRBVar sdiff = model->addVar(-MAX_EDGE_WEIGHT, MAX_EDGE_WEIGHT, 0, GRB_CONTINUOUS);
                model->addConstr(sdiff, GRB_EQUAL, S[i][j]-S[i][k]);
                
                GRBVar sabs = model->addVar(0, MAX_EDGE_WEIGHT, 0, GRB_CONTINUOUS);
                model->addGenConstrAbs(sabs,sdiff);
                
                expr = MAX_EDGE_WEIGHT*(2-Y[i][j]-Y[i][k]);
                sprintf(name, "ys-%d-%d-%d", i,j,k);//simple path
                model->addConstr(sabs, GRB_LESS_EQUAL, expr, name);
                cnt++;
            }
        }
    }
    
    printf("Total %d y to s constraints.\n", cnt);
    return 0;
}

/*int ilp::add_s_to_t_constraint()
{
    char name[1024];
    int cnt = 0;
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < M; j++)
        {
            vector<int> phasing_edge = H[j];
            vector<int>::const_iterator it;
            //printf("Phasing path %d contains:\n", j);
            for(it = phasing_edge.begin(); it != phasing_edge.end(); it++)
            {
                sprintf(name, "st-%d-%d-%d", i,j,*it);//simple path
                model->addConstr(S[i][j], GRB_LESS_EQUAL, T[i][*it], name);
                cnt++;
            }
        }
    }

    printf("Total %d s to t constraints.\n", cnt);
    return 0;
}*/

int ilp::add_s_to_t_constraint()
{
    char name[1024];
    int cnt = 0;
    for(int i = 0; i < N; i++)
    {
        for(int k = 0; k < gr.num_edges(); k++)
        {
            GRBLinExpr expr;
            set<int>::const_iterator it;
            set<int> phasing = map_e_to_phasing[k];
            for(it=phasing.begin(); it != phasing.end(); it++)
            {
                expr += S[i][*it];
                model->addConstr(S[i][*it], GRB_LESS_EQUAL, T[i][k]);
                //model->addConstr(S[i][*it], GRB_EQUAL, T[i][k]);
                
                edge_descriptor e = i2e[k];
                //printf("[%d -> %d] supported by: \n", e->source(), e->target());
                //print_phasing_path(H[*it]);
            }
            //model->addConstr(expr, GRB_LESS_EQUAL, T[i][k]*(1+PHASING_DIFF));
            //if(!phasing.empty())
                //model->addConstr(expr, GRB_GREATER_EQUAL, T[i][k]*(1-PHASING_DIFF)-0.5);
        }
    }
    return 0;
}

/*int ilp::add_s_to_t_constraint()
{
    char name[1024];
    int cnt = 0;
    for(int i = 0; i < N; i++)
    {
        GRBLinExpr expr1;
        for(int k = 0; k < gr.num_edges(); k++)
        {
            expr1 += T[i][k];
        }

        GRBLinExpr expr2;
        for(int j = 0; j < M; j++)
        {
            expr2 += S[i][j]*(H[j].size());
            //cout << H[j].size() << endl;
        }
        model->addConstr(expr1, GRB_EQUAL, expr2);
    }
    return 0;
}*/

int ilp::add_unbalanced_junction_constraint()
{
    char name[1024];
    int cnt = 0;
    for(int j = 0; j < gr.num_edges(); j++)
    {
        //vector<double> eubj;
        edge_descriptor e1 = i2e[j];
        int v = e1->target();
        double ew1 = gr.get_edge_weight(e1);
        for(int k = 0; k < gr.num_edges(); k++)
        {
            edge_descriptor e2 = i2e[k];
            if(v == e2->source())
            {
                double in_max = gr.get_max_in_weight(v);
                double out_max = gr.get_max_out_weight(v);
                double ew2 = gr.get_edge_weight(e2);
                //double ubj = (in_max-ew1) + (out_max-ew2)+fabs(ew1-ew2);
                double ubj = max(in_max-ew1, out_max-ew2);
                //eubj.push_back(ubj)
                //printf("[%d -> %d -> %d] U: %f\n", e1->source(), v, e2->target(), ubj);
                for(int i = 0 ; i < N; i++)
                {
                    model->addConstr(UBJ[i][j][k], GRB_GREATER_EQUAL, ubj*(X[i][j]+X[i][k]-1));
                }
            }
            else
            {
                for(int i = 0; i < N; i++)
                {
                    model->addConstr(UBJ[i][j][k], GRB_EQUAL, 0.0);
                }
                //eubj.push_back(0.0);
            }
        }
    }
    return 0;
}

int ilp::add_duplicate_path_constraint()
{
    char name[1024];
    int cnt = 0;
    for(int i = 0; i < N; i++)
    {
        for(int j = i+1; j < N; j++)
        {
            GRBLinExpr expr;
            for(int k = 0; k < gr.num_edges(); k++)
            {
                GRBVar diff = model->addVar(-1, 1, 0, GRB_CONTINUOUS);
                model->addConstr(diff, GRB_EQUAL, X[i][k]-X[j][k]);
                GRBVar dabs = model->addVar(0, 1, 0, GRB_CONTINUOUS);
                model->addGenConstrAbs(dabs, diff);
                expr += dabs;
            }
            model->addConstr(D[i][j-i-1], GRB_LESS_EQUAL, expr);
        }
    }
    return 0;
}

int ilp::add_adjacent_edge_constraint()
{
    for(int v = 0; v < gr.num_vertices(); v++)
    {
        if(gr.in_degree(v) == 1 || gr.out_degree(v) == 1)//trivial vertices
            continue;
        PEEI pei1 = gr.in_edges(v);
        PEEI pei2 = gr.out_edges(v);

        for(edge_iterator ei1 = pei1.first; ei1 != pei1.second; ei1++)
        {
            if(gr.get_max_in_weight(v) - gr.get_edge_weight(*ei1)<1)
                continue;

            for(edge_iterator ei2 = pei2.first; ei2 != pei2.second; ei2++)
            {
                if(gr.get_max_out_weight(v) - gr.get_edge_weight(*ei2)<1)
                    continue;
                int eid1 = e2i[*ei1], eid2 = e2i[*ei2];
                set<int> phasing1 = map_e_to_phasing[eid1];
                set<int> phasing2 = map_e_to_phasing[eid2];
                set<int> common;
                set_intersection(phasing1.begin(), phasing1.end(), phasing2.begin(), phasing2.end(), inserter(common,common.end()));

                if(common.size() == 0)
                {
                    printf("Edge without junction: [%d->%d] and [%d->%d]\n", (*ei1)->source(), (*ei1)->target(), (*ei2)->source(), (*ei2)->target());
                    for(int i = 0; i < N; i++)
                    {
                        model->addConstr(X[i][eid1]+X[i][eid2], GRB_LESS_EQUAL, 1);
                    }
                }


            }
        }
    }
    return 0;
}

int ilp::allele_phase_constraint()
{
    int NUM_PHASE_SET = 2; // for human, number of phase set is 2 (diploid)

	vector<PI32> as_exons;
	for(int i = 0; i < gr.vinf.size(); ++i)
	{
		for (auto v: gr.vinf)
		{
            as_pos32 j = v.lpos;
			if (j.ale !="$") 
			{
				as_exons.push_back(make_pair(v.lpos, v.rpos));
			}
		}
	}
	sort(as_exons.begin(), as_exons.end());
	map<PI32, int> as_exons_enum;
	for(int i = 0; i < as_exons.size(); ++i) 
	{
		as_exons_enum.insert(make_pair(as_exons[i], i));
	}

	// map<int, vector<int> > trs_as_exon_map; 
	// for(int i = 0; i < N; ++i)
	// {
	// 	for (int k = 0; k < as_exons.size(); ++k)
	// 	{
	// 		if (j.ale !="$") 
	// 		{
	// 			e.push_back(as_exons_enum.find(make_pair(v.lpos, v.rpos))->second);
	// 		}
	// 	}
	// 	trs_as_exon_map.insert(make_pair(i, e));
	// }

	// map<int, vector<int> > as_exon_trs_map;
	// for (auto it = trs_as_exon_map.begin(); it != trs_as_exon_map.end(); ++it)
	// {
	// 	vector<int> asexons = it->second;
	// 	for (auto jt = asexons.begin(); jt != asexons.end(); ++jt)
	// 	{
	// 		if (as_exon_trs_map.find(*jt) == as_exon_trs_map.end())
	// 		{
	// 			vector<int> v;
	// 			v.push_back(it->first);
	// 			as_exon_trs_map.insert(make_pair(*jt, v));
	// 		}
	// 		else 
	// 		{
	// 			as_exon_trs_map.find(*jt)->second.push_back(it->first);
	// 		}
	// 	}
	// }
	// try 
	{		
		// Create variables, all variables are binary
		// x_{i,k} 		- whether i-th AS-exon is in k-th phase set 
		// v_i	 		- whether i-th AS-exon is retained after filter
		// y_m 			- whether m-th transcript is retained after filter
		// z_{m,k} 		- whether k-th phase set satisfies m-th transcript
		GRBVar x[as_exons.size()][NUM_PHASE_SET];
		GRBVar vi[as_exons.size()];
		GRBVar y[N];
		GRBVar z[N][NUM_PHASE_SET];
		// GRBLinExpr obj = GRBLinExpr();
		double total_cov = 0;
		// for (int m = 0; m < N; ++m) 
		{
			// total_cov += trs[m].coverage;
		}
		for (int i = 0; i < as_exons.size(); ++i)
		{
			for (int k = 0; k < NUM_PHASE_SET; ++k)
			{
				x[i][k] = model->addVar(0, 1, 0, GRB_BINARY);
			}
		}
		for (int i = 0; i < as_exons.size(); ++i) 
		{
			vi[i] = model->addVar(0, 1, 0, GRB_BINARY);
			// obj += vi[i];
		}
		for (int m = 0; m < N; ++m) 
		{
			y[m] = model->addVar(0, 1, 0, GRB_BINARY);
			// obj += trs[m].coverage * y[m] / total_cov;
		}
		for (int m = 0; m < N; ++m) 
		{
			for (int k = 0; k < NUM_PHASE_SET; ++k)
			{
				z[m][k] = model->addVar(0, 1, 0, GRB_BINARY);
			}
		}

		// Set Objective y[m].coverage * y[m]
		// model.setObjective(obj, GRB_MAXIMIZE);

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
						model->addConstr(x[i][k] + x[i + j][k] <= 1);
					}
				}
			}
		}

		// Transcript phase set constraint
		// A retained transcript must be satisfied by at least one phase set
		int m = 0;
		for (int m = 0; m < N; ++m)
		{
            for (int k = 0; k < NUM_PHASE_SET; ++k)
            {
                GRBLinExpr expr_xz = GRBLinExpr();  
                int n = 0;
                for (int j = 0; j < e2i.size(); ++j)
                {
                    edge_descriptor d = i2e[j];
                    vertex_info v1 = gr.vinf[d->source()];
                    vertex_info v2 = gr.vinf[d->target()];
                    auto p1 = as_exons_enum.find(make_pair(v1.lpos, v1.rpos));
                    auto p2 = as_exons_enum.find(make_pair(v2.lpos, v2.rpos));
                    if (p1 != as_exons_enum.end())
                    {
                        expr_xz += x[p1->second][k]; 
                        n++;
                    }
                    if (p2 != as_exons_enum.end())
                    {
                        expr_xz += x[p2->second][k]; 
                        n++;
                    }
                }
                expr_xz = expr_xz /double(n);
                expr_xz -= z[m][k];
                model->addConstr(expr_xz, GRB_GREATER_EQUAL, 0);  // Sum_i(x[i][k]) / n - z[m][k] >= 0
            }

            GRBLinExpr expr_zy = GRBLinExpr();
            for (int k = 0; k < NUM_PHASE_SET; ++k)
            {
                expr_zy += z[m][k];
            }
            expr_zy -= y[m];
            model->addConstr(expr_zy, GRB_GREATER_EQUAL, 0);  // Sum_k(z[m][k]) - y[m] >= 0 		
		}
		// assert(m == trs_as_exon_map.size());



		// for (auto it = as_exon_trs_map.begin(); it != as_exon_trs_map.end(); ++it)
		// {
		// 	GRBLinExpr expr_v = GRBLinExpr(); 
		// 	for (int ii = 0; ii < it->second.size(); ii ++)
		// 	{
		// 		expr_v += y[it->second[ii]];
		// 	}
		// 	model.addConstr(expr_v, GRB_GREATER_EQUAL, vi[it->first]);
		// }
		

		// Optimize
		// model.optimize();
		// vector<transcript> v;
		// for (int m = 0; m < trs.size(); ++m)
		// {
		// 	if (NUM_PHASE_SET == 2) trs[m].allele = to_string(int(round(z[m][0].get(GRB_DoubleAttr_X))));

		// 	if (y[m].get(GRB_DoubleAttr_X) >= 0.999) v.push_back(trs[m]); // precision issue
		// 	if (verbose >= 2) cout << m << " y_m: " <<y[m].get(GRB_DoubleAttr_X) << endl;
		// 	assert((y[m].get(GRB_DoubleAttr_X) <= 0.001 && y[m].get(GRB_DoubleAttr_X) >= -0.001) 
		// 		|| (y[m].get(GRB_DoubleAttr_X) >= 0.999 && y[m].get(GRB_DoubleAttr_X) <= 1.001)); 
		// 	if (verbose >= 2)
		// 	{
		// 		for (int k = 0; k < NUM_PHASE_SET; ++k) 
		// 			cout << m << ' ' <<k << " z_mk: " <<z[m][k].get(GRB_DoubleAttr_X) << endl;
		// 	}
		// }
		/* // print x_ik values
		for (int k = 0; k < NUM_PHASE_SET; ++k) 
		{
			for (int i = 0; i < as_exons.size(); ++i)
			{
				cout << i << ' ' <<k << " x_ik: " <<x[i][k].get(GRB_DoubleAttr_X) << endl;
			}
		}
		*/
		// trs = v;
	}
	// catch (GRBException e)
	// {
	// 	cout << "Error code = " << e.getErrorCode() << endl;
    // 	cout << e.getMessage() << endl;
	// 	throw BundleError();
	// }

	return 0;
}

int ilp::set_objective()
{
    char name[1024];
    int cnt = 0;
    GRBLinExpr expr_obj;

    /*GRBLinExpr expr_duplicate_path;
    for(int i = 0; i < N; i++)
    {
        for(int j = i+1; j < N; j++)
        {
            expr_duplicate_path += D[i][j-i-1];
        }
    }
    expr_obj -= expr_duplicate_path * gr.get_edge_weight(gr.compute_maximum_edge_w());*/

    //penalty on edges that are not phased in the predicted path
    /*GRBLinExpr expr_non_phasing_edge;
    for(int i = 0; i < N; i++)
    {
        for(int k = 0; k < gr.num_edges(); k++)
        {
            GRBLinExpr expr;
            set<int>::const_iterator it;
            set<int> phasing = map_e_to_phasing[k];
            edge_descriptor e = i2e[k];

            if(phasing.size() == 0) continue;
            
            GRBVar edge_phasing= model->addVar(0, 1, 0, GRB_BINARY);
            for(it=phasing.begin(); it != phasing.end(); it++)
            {
                expr += Y[i][*it];
                model->addConstr(edge_phasing, GRB_GREATER_EQUAL, Y[i][*it]);

                //edge_descriptor e = i2e[k];
                //printf("[%d -> %d] supported by: \n", e->source(), e->target());
                //print_phasing_path(H[*it]);
            }
            model->addConstr(edge_phasing, GRB_LESS_EQUAL, expr);
            model->addConstr(NPE[i][k], GRB_EQUAL, (X[i][k]-edge_phasing)*gr.get_edge_weight(gr.compute_maximum_edge_w())/gr.get_edge_weight(e));
            expr_non_phasing_edge += NPE[i][k];
            if(gr.get_edge_weight(e)<pow(PHASING_PATH_THRESHOLD,2))
                model->addConstr(T[i][k], GRB_LESS_EQUAL,0.5+ (1-X[i][k]+edge_phasing)*MAX_EDGE_WEIGHT);

        }
    }
    expr_obj += expr_non_phasing_edge;*/

    GRBLinExpr expr_phasing_cover;
    for(int i = 0; i < M; i++)
    {
        expr_phasing_cover += H_weight[i]*(1-PHASING_SUPPORT[i]);
    }
    expr_obj += expr_phasing_cover;
    obj_phasing_cover = model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
    model->addConstr(obj_phasing_cover, GRB_EQUAL, expr_phasing_cover);

    GRBLinExpr expr_unbalanced_junction;
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < gr.num_edges(); j++)
        {
            for(int k = 0; k < gr.num_edges(); k++)
            {
                int v1 = i2e[j]->target();
                int v2 = i2e[k]->source();
                int d1 = gr.in_degree(v1);
                int d2 = gr.out_degree(v2);

                expr_unbalanced_junction += UBJ[i][j][k];
                //expr_unbalanced_junction += UBJ[i][j][k]/(d1*d2);
            }
        }
    }
    expr_obj += expr_unbalanced_junction/gr.compute_average_edge_weight();
    //expr_obj += expr_unbalanced_junction*0.5;
    obj_UBJ = model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
    model->addConstr(obj_UBJ, GRB_EQUAL, expr_unbalanced_junction);



    GRBLinExpr expr_phasing_count;
    for(int j = 0; j < M; j++)
    {
        GRBLinExpr expr;
        for(int i = 0; i < N; i++)
        {
            expr += S[i][j];
        }
        
        sprintf(name, "sz-%d", j);//s to z
        model->addConstr(Z[j], GRB_EQUAL, expr, name);
        cnt++;
        GRBVar wz_diff = model->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS);
        model->addConstr(wz_diff, GRB_EQUAL, H_weight[j]-Z[j]);
        //cout << H_weight[j] << endl;
        
        GRBVar wz_abs = model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
        model->addGenConstrAbs(wz_abs,wz_diff);
        
        expr_phasing_count += wz_abs/sqrt(H_weight[j]);
        //expr_phasing_count += wz_abs;
        //expr_phasing_count += H_weight[j]-Z[j];

    }
    //expr_obj += expr_phasing_count*1.0*(gr.num_edges()+gr.num_vertices())/M;
    expr_obj += expr_phasing_count*1.0*gr.num_edges()/M;
    obj_phasing_count = model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
    model->addConstr(obj_phasing_count, GRB_EQUAL, expr_phasing_count);


    /*GRBLinExpr expr_phasing_balance;
    for(int i = 0; i < N; i++)
    {
        for(int k = 0; k < gr.num_edges(); k++)
        {
            GRBLinExpr expr;
            set<int>::const_iterator it;
            set<int> phasing = map_e_to_phasing[k];
            if(phasing.size() == 0) continue;
            for(it=phasing.begin(); it != phasing.end(); it++)
            {
                expr += S[i][*it];
                
                edge_descriptor e = i2e[k];
                //printf("[%d -> %d] supported by: \n", e->source(), e->target());
                //print_phasing_path(H[*it]);
            }
             GRBVar st_diff = model->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS);
             model->addConstr(st_diff, GRB_EQUAL, expr-T[i][k]);
             GRBVar st_abs = model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
             model->addGenConstrAbs(st_abs,st_diff);
             expr_phasing_balance += st_abs;
            //expr_phasing_balance += T[i][k]-expr;
        }
    }
    expr_obj += expr_phasing_balance;
    //expr_obj += expr_phasing_balance*1.0*(gr.num_edges()+gr.num_vertices())/(N*gr.num_edges());
    //expr_obj += expr_phasing_balance*1.0*gr.num_edges()/(N*gr.num_edges());*/

    GRBLinExpr expr_edge_weight;
    double max_edge = gr.get_edge_weight(gr.compute_maximum_edge_w());
    for(int j = 0; j < gr.num_edges(); j++)
    {
        edge_descriptor e = i2e[j];
        int s = e->source(), t = e->target();
        double weight = gr.get_edge_weight(e);

        //if(e->source()==0 || e->target()==gr.num_vertices()-1) continue;

        GRBLinExpr expr;
        for(int i = 0; i < N; i++)
        {
            expr += T[i][j];
        }
        sprintf(name, "te-%d", j);
        
        double penalty = 2;
        double alpha1 = 1, alpha2 = 0;
        if(weight >= 2*pow(PHASING_PATH_THRESHOLD,2))
        {
            //penalty = 3;
            alpha1 = (0.5+penalty)*0.5;
            alpha2 = (penalty-0.5)*0.5;

        }
        else if(weight >= pow(PHASING_PATH_THRESHOLD,2))
        {
            //penalty = sqrt(weight)/PHASING_PATH_THRESHOLD;
            //penalty = 2;
            alpha1 = (1+penalty)*0.5;
            alpha2 = (penalty-1)*0.5;
            //model->addConstr(EDGE_SUPPORT[j], GRB_EQUAL, expr, name);
        }
        else if(weight <= PHASING_PATH_THRESHOLD)
        {
            //penalty = PHASING_PATH_THRESHOLD/weight;
            //penalty = 2;
            SQRTDIFF_E[j] *= penalty;
            double p = sqrt(weight)/PHASING_PATH_THRESHOLD;
            alpha1 = (p+SQRTDIFF_E[j])*0.5;
            alpha2 = (p-SQRTDIFF_E[j])*0.5;

            //model->addConstr(EDGE_SUPPORT[j], GRB_EQUAL,expr, name);
        }
        else
        {
            double p = sqrt(weight)/PHASING_PATH_THRESHOLD;
            alpha1 = (p+SQRTDIFF_E[j])*0.5;
            alpha2 = (p-SQRTDIFF_E[j])*0.5;
        }
        model->addConstr(EDGE_SUPPORT[j], GRB_EQUAL, expr*SQRTDIFF_E[j], name);
        //model->addConstr(EDGE_SUPPORT[j], GRB_LESS_EQUAL, weight*2);
        //cout << e->source() << "->" << e->target() << " penalty:" << penalty << endl;

        GRBVar we_diff = model->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS);
        model->addConstr(we_diff, GRB_EQUAL, weight-EDGE_SUPPORT[j]);
        
        GRBVar we_abs = model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
        model->addGenConstrAbs(we_abs,we_diff);

        //cout << alpha1 << '\t' << alpha2 << endl;
         expr_edge_weight += (we_abs*alpha1 + alpha2*(weight-EDGE_SUPPORT[j]))/sqrt(weight)*(t-s);
        //expr_edge_weight += (we_abs*alpha1 + alpha2*(-EDGE_SUPPORT[j]))*(t-s);
        
        //GRBVar ratio = model->addVar(0.5, 1, 0, GRB_CONTINUOUS);
        //Ratio.push_back(ratio);
        //double weight = gr.get_edge_weight(i2e[j]);
        //model->addConstr(expr, GRB_GREATER_EQUAL, ratio*weight);
        //expr_obj += expr-ratio*weight;

        //cout << "Edge weight=" << weight << "\talpha1=" << alpha1 << "\talpha2=" << alpha2 << endl;

    }
    expr_obj += expr_edge_weight;
    obj_edge_weight = model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
    model->addConstr(obj_edge_weight, GRB_EQUAL, expr_edge_weight);

    /*GRBLinExpr expr_vertex_weight;
    for(int v = 0; v < gr.num_vertices(); v++)
    {
        if(gr.in_degree(v)==0 || gr.out_degree(v)==0) continue;
        PEEI pei = gr.out_edges(v);
        GRBLinExpr expr_v_sup;
        for(edge_iterator ei1 = pei.first, ei2 = pei.second; ei1 != ei2; ei1++)
        {
            int eid = e2i[*ei1];
            expr_v_sup += EDGE_SUPPORT[eid];
        }
        
        double penalty = 2;
        double weight = gr.get_vertex_weight(v);
        double alpha1 = 1, alpha2 = 0;
        if(weight >= 2*pow(PHASING_PATH_THRESHOLD,2))
        {
            //penalty = 3;
            alpha1 = (0.5+penalty)*0.5;
            alpha2 = (penalty-0.5)*0.5;

        }
        else if(weight >= pow(PHASING_PATH_THRESHOLD,2) && SQRTDIFF[v]<1.01)
        {
            //penalty = sqrt(weight)/PHASING_PATH_THRESHOLD;
            //penalty = 2; 
            alpha1 = (1+penalty) *0.5;
            alpha2 = (penalty-1)*0.5;

        }
        else if(weight <= PHASING_PATH_THRESHOLD)
        {
            //penalty = sqrt(PHASING_PATH_THRESHOLD/weight);
            //penalty = 2;
            //if(SQRTDIFF[v]<penalty) SQRTDIFF[v] = penalty;
            SQRTDIFF[v] *= penalty;
            
            alpha1 = (1+SQRTDIFF[v])*0.5;
            alpha2 = (1-SQRTDIFF[v])*0.5; 
        }
        else
        {
            alpha1 = (1+SQRTDIFF[v])*0.5;
            alpha2 = (1-SQRTDIFF[v])*0.5;
        }

        sprintf(name, "v-sup-%d", v);//s to z
        model->addConstr(VERTEX_SUPPORT[v], GRB_EQUAL, expr_v_sup*SQRTDIFF[v], name);
        cnt++;

        GRBVar wv_diff = model->addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS);
        model->addConstr(wv_diff, GRB_EQUAL, weight-VERTEX_SUPPORT[v]);
        //cout << gr.get_vertex_weight(v) << endl;

        GRBVar wv_abs = model->addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
        model->addGenConstrAbs(wv_abs,wv_diff);

        //double alpha1 = 0.5*(SQRTDIFF[v]+1);
        //double alpha2 = 0.5*(1-SQRTDIFF[v]);
        expr_vertex_weight += (wv_abs*alpha1 + alpha2*(weight-VERTEX_SUPPORT[v]))/sqrt(weight);
        
        //cout << "Vertex weight=" << weight << "\talpha1=" << alpha1 << "\talpha2=" << alpha2 << endl;

    }
    expr_obj += expr_vertex_weight;*/


    model->setObjective(expr_obj, GRB_MINIMIZE);
    return 0;
}

int ilp::collect_results(int suggestN)
{
    map< vector<int>, double> map_path;
    for(int i = 0; i < N; i++)
    {
        printf("\n-------Predicted path %d-------\n", i);
        set<int> spath;
        double abd = 0;
        for(int j = 0; j < gr.num_edges(); j++)
        {
            GRBVar x = X[i][j];
            edge_descriptor e = i2e[j];
            int s = e->source();
            int t = e->target();
            abd = max(abd, T[i][j].get(GRB_DoubleAttr_X));
            if(x.get(GRB_DoubleAttr_X) > 0.9)
            {
                //printf("[%d->%d] Predict: %.2f\n", s,t,T[i][j].get(GRB_DoubleAttr_X));
                spath.insert(s);
                spath.insert(t);
            }
        }

        cout << "Abundance: " << abd << endl;
        //if(abd<0.5) continue;
        if(*spath.begin() != 0) spath.insert(0);
        if(*spath.rbegin() != gr.num_vertices()-1) spath.insert(gr.num_vertices()-1);

        set<int>::const_iterator vit;
        vector<int> vpath;
        //path p;
        for(vit=spath.begin(); vit != spath.end(); vit++)
        {
            if(verbose >= 1)
            {
                if(vit == spath.begin())
                {
                    printf("%d ", *vit);
                }
                else
                {
                    printf("-> %d ", *vit);
                }
            }
            vpath.push_back(*vit);
        }
        printf("\n");

        printf("Edge not phasing:\n");
        for(int j = 0; j < gr.num_edges(); j++)
        {
            if(NPE[i][j].get(GRB_DoubleAttr_X)>0.1 && X[i][j].get(GRB_DoubleAttr_X)>0.9)
            {
                edge_descriptor e = i2e[j];
                cout << e->source() << "->" << e->target() << endl;
            }
        }

        //if(abd<=2) continue;
        //if(abd<=PHASING_PATH_THRESHOLD) continue;

        if(map_path.find(vpath) == map_path.end())
        {
            map_path.insert(make_pair(vpath, abd));
        }
        else
        {
            map_path[vpath] += abd;
            cout << "Duplicate path!" << endl;
            continue;
        }
        
        //cout << map_path.size() << endl;

        if(verbose >= 1)
        {
            printf("\nPhasing paths:\n");
            for(int j = 0; j < M; j++)
            {
                GRBVar y = Y[i][j];
                GRBVar s = S[i][j];
                //cout << y.get(GRB_DoubleAttr_X) << endl;
                if(y.get(GRB_DoubleAttr_X) > 0.9)
                {
                    //cout << y.get(GRB_StringAttr_VarName) << endl;
                    cout << y.get(GRB_StringAttr_VarName) << ":\t";
                    print_phasing_path(H[j]);
                    cout << "\tWeight: " << s.get(GRB_DoubleAttr_X) << endl;

                }
            }
        }
    }
    
    map< vector<int>, double >::const_iterator it;
    set< pair< double, vector<int> > > set_abd_path;
    for(it = map_path.begin(); it != map_path.end(); it++)
    {
        /*path p;
        p.v = it->first;
        p.abd = it->second;
        paths.push_back(p);*/
        set_abd_path.insert(make_pair(it->second, it->first));
    }

    set< pair< double, vector<int> > >::reverse_iterator it_abd;
    int n = 0;
    paths.clear();
    for(it_abd = set_abd_path.rbegin(); it_abd != set_abd_path.rend() &&  n<N2; it_abd++, n++)
    {
        path p;
        p.v = it_abd->second;
        p.abd = it_abd->first;
        //cout << p.abd << endl;
        paths.push_back(p);

    }

    printf("\n");
    for(int j = 0; j < M; j++)
    {
        GRBVar z = Z[j];
        cout << z.get(GRB_StringAttr_VarName) << " "
        << z.get(GRB_DoubleAttr_X) << '/' << H_weight[j] << endl;
        print_phasing_path(H[j]);
        printf("\n");
    }

    if(verbose >= 1)
    {
        printf("\n");
        for(int j = 0; j < gr.num_edges(); j++)
        {
            GRBVar edge_sup = EDGE_SUPPORT[j];
            edge_descriptor e = i2e[j];
            cout << e->source() << "->" << e->target() << " " << edge_sup.get(GRB_DoubleAttr_X)/SQRTDIFF_E[j] << '/' << gr.get_edge_weight(e) << endl;
            //cout << e->source() << "->" << e->target() << " " << edge_sup.get(GRB_DoubleAttr_X) << '/' << gr.get_edge_weight(e) << "\tRatio: " << Ratio[j].get(GRB_DoubleAttr_X) << endl;
        }

        printf("\n");
        /*for(int v = 0; v < gr.num_vertices(); v++)
        {
            GRBVar vertex_sup = VERTEX_SUPPORT[v];
            cout << v << ": " << vertex_sup.get(GRB_DoubleAttr_X)/SQRTDIFF[v] << '/' << gr.get_vertex_weight(v) << endl;
        }*/
    }
    
    return 0;
}
