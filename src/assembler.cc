/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <algorithm>
#include <cstdio>
#include <cassert>
#include <sstream>
#include <iostream>
#include <map>
#include <cstring>

#include "config.h"
#include "genome.h"
#include "assembler.h"
#include "bundle.h"
#include "scallop.h"
#include "sgraph_compare.h"
#include "super_graph.h"
#include "filter.h"
#include "vcf_data.h"
#include "phaser.h"
#include "util.h"

assembler::assembler()
{
    sfn = sam_open(input_file.c_str(), "r");
    hdr = sam_hdr_read(sfn);
    b1t = bam_init1();
	hid = 0;
	index = 0;
	terminate = false;
	qlen = 0;
	qcnt = 0;
	vmap_chrm = "";  // reset vcf pointers from previewer
	trsts.resize(3);
	nonfull_trsts.resize(3);
	if(DEBUG_MODE_ON)
	{
		assert(trsts.size() == 3);
		assert(nonfull_trsts.size() == 3);
		for(const auto & i: trsts) assert(i.size() == 0);
		for(const auto & i: nonfull_trsts) assert(i.size() == 0);
	}
}

assembler::~assembler()
{
    bam_destroy1(b1t);
    bam_hdr_destroy(hdr);
    sam_close(sfn);
	fai_destroy(fai);
}

int assembler::assemble()
{
    while(sam_read1(sfn, hdr, b1t) >= 0)
	{
		if(terminate == true) return 0;

		bam1_core_t &p = b1t->core;

		if(p.tid < 0) continue;
		if((p.flag & 0x4) >= 1) continue;										// read is not mapped
		if((p.flag & 0x100) >= 1 && use_second_alignment == false) continue;	// secondary alignment
		if(p.n_cigar > max_num_cigar) continue;									// ignore hits with more than max-num-cigar types
		if(p.qual < min_mapping_quality) continue;								// ignore hits with small quality
		if(p.n_cigar < 1) continue;												// should never happen


		char buf[1024];
		strcpy(buf, hdr->target_name[p.tid]);

		hit ht(b1t, string(buf), hid++);
		ht.set_tags(b1t);
		ht.set_strand();
		
		qlen += ht.qlen;
		qcnt += 1;

		// truncate
		if(ht.tid != bb1.tid || ht.pos > bb1.rpos + min_bundle_gap)
		{
			if (bb1.hits.size() >= 1) pool.push_back(bb1);
			bb1.clear();
		}
		if(ht.tid != bb2.tid || ht.pos > bb2.rpos + min_bundle_gap)
		{
			if(bb2.hits.size() >= 1) pool.push_back(bb2);
			bb2.clear();
		}

		// process
		process(batch_bundle_size);

		//printf("read strand = %c, xs = %c, ts = %c\n", ht.strand, ht.xs, ht.ts);

		// add hit
		if(uniquely_mapped_only == true && ht.nh != 1) continue;
		if(library_type != UNSTRANDED && ht.strand == '+' && ht.xs == '-') continue;
		if(library_type != UNSTRANDED && ht.strand == '-' && ht.xs == '+') continue;
		if(library_type != UNSTRANDED && ht.strand == '.' && ht.xs != '.') ht.strand = ht.xs;
		if(library_type != UNSTRANDED && ht.strand == '+') bb1.add_hit(ht);
		if(library_type != UNSTRANDED && ht.strand == '-') bb2.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '.') bb1.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '.') bb2.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '+') bb1.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '-') bb2.add_hit(ht);
	}

	pool.push_back(bb1);
	pool.push_back(bb2);
	process(0);

	assign_RPKM();

	if(DEBUG_MODE_ON && trsts[0].size() < 1 && trsts[1].size() < 1 && trsts[2].size() < 1) throw runtime_error("No AS transcript found!");

	// filter each allele
	for(int i = 0; i < 3; i++)
	{
		if (! use_filter) break;

		filter ft(trsts[i]);
		ft.merge_single_exon_transcripts();
		trsts[i] = ft.trs;

		filter ft1(nonfull_trsts[i]);
		ft1.merge_single_exon_transcripts();
		nonfull_trsts[i] = ft1.trs;
	}	

	write();
	
	return 0;
}

int assembler::process(int n)
{
	if(pool.size() < n) return 0;
	for(int i = 0; i < pool.size(); i++)
	{
		bundle_base &bb = pool[i];
		bb.buildbase();

		if(verbose >= 3) printf("bundle %d has %lu reads\n", index, bb.hits.size());

		int cnt1 = 0;
		int cnt2 = 0;
		for(int k = 0; k < bb.hits.size(); k++)
		{
			//counts += (1 + bb.hits[k].spos.size());
			if(bb.hits[k].spos.size() >= 1) cnt1 ++;
			else cnt2++;
		}
		if(cnt1 + cnt2 < min_num_hits_in_bundle) continue;

		if(bb.tid < 0) continue;

		char buf[1024];
		strcpy(buf, hdr->target_name[bb.tid]);
		bb.chrm = string(buf);

		// transcript_set ts1(bb.chrm, 0.9);	
		// transcript_set ts2(bb.chrm, 0.9);		
		vector<transcript_set> ts_full;		 	// full-length set; [0]NONSPECIFIC/UNPHASED, [1]ALLELE1, [2]ALLELE2
		vector<transcript_set> ts_nonfull;		// non-full-length set; [0]NONSPECIFIC/UNPHASED, [1]ALLELE1, [2]ALLELE2
		for(int i = 0; i < 3; i++) ts_full.push_back(transcript_set(bb.chrm, 0.9));
		for(int i = 0; i < 3; i++) ts_nonfull.push_back(transcript_set(bb.chrm, 0.9));

		bundle bd(bb);
		bd.build(1, true);
		bd.print(index++);
		assemble(bd.gr, bd.hs, bb.is_allelic, ts_full, ts_nonfull);

		bd.build(2, true);
		bd.print(index++);				
		assemble(bd.gr, bd.hs, bb.is_allelic, ts_full, ts_nonfull);


		
		// retrieve and filter transcripts
		// i = {0, 1, 2}, corresponds to NONSPECIFIC/UNPHASED, ALLELE1, ALLELE2
		for(int i = 0; i < 3; i++)
		{
			if (DEBUG_MODE_ON)	assert(trsts.size() == 3 && nonfull_trsts.size() == 3);
			
			int sdup = assemble_duplicates / 1 + 1;
			int mdup = assemble_duplicates / 2 + 0;
			vector<transcript> gv1 = ts_full[i].get_transcripts(sdup, mdup); 
			vector<transcript> gv2 = ts_nonfull[i].get_transcripts(sdup, mdup);

			for(int k = 0; k < gv1.size(); k++)
			{
				if(gv1[k].exons.size() >= 2) gv1[k].coverage /= (1.0 * assemble_duplicates);
			}
			for(int k = 0; k < gv2.size(); k++) 
			{
				if(gv2[k].exons.size() >= 2) gv2[k].coverage /= (1.0 * assemble_duplicates);
			}

			if (use_filter)
			{
				filter ft1(gv1);
				ft1.filter_length_coverage();
				ft1.remove_nested_transcripts();
				if(ft1.trs.size() >= 1) trsts[i].insert(trsts[i].end(), ft1.trs.begin(), ft1.trs.end());

				filter ft2(gv2);
				ft2.filter_length_coverage();
				ft2.remove_nested_transcripts();
				if(ft2.trs.size() >= 1) nonfull_trsts[i].insert(nonfull_trsts[i].end(), ft2.trs.begin(), ft2.trs.end());
			}
			else
			{
				trsts[i].insert(trsts[i].end(), gv1.begin(), gv1.end());
				nonfull_trsts[i].insert(nonfull_trsts[i].end(), gv2.begin(), gv2.end());
			}
		}		
	}
	pool.clear();
	return 0;
}

int assembler::assemble(const splice_graph &gr0, const hyper_set &hs0, bool is_allelic, vector<transcript_set> &ts_full, vector<transcript_set> &ts_nonfull)
{
	if(DEBUG_MODE_ON)
	{
		for (int i = 0  ;  i < gr0.vinf.size(); i++)
		{
			cout << "gr0 bef scallop first round: " << i << " " << gt_str(gr0.vinf[i].gt) << endl;
		}
	}

	super_graph sg(gr0, hs0);
	sg.build();

	for(int k = 0; k < sg.subs.size(); k++)
	{
		splice_graph& gr_alias = sg.subs[k];
		hyper_set &hs = sg.hss[k];

		if(determine_regional_graph(gr_alias)) continue;

		if(gr_alias.num_edges() <= 0) continue;
		if(debug_bundle_only) continue; //debug parameter to build bundle only and skip assembly, default: false
		
		for(int r = 0; r < assemble_duplicates; r++)
		{
			splice_graph gr = gr_alias;  // graph copy for different duplicates, copy constructor used, not move
			string gid = "gene." + tostring(index) + "." + tostring(k) + "." + tostring(r);
			gr.gid = gid;

			// partial decomp of non-AS nodes
			scallop sc(gr, hs, r == 0 ? false : true, true);
			sc.assemble(is_allelic);  

			if(verbose >=3 && DEBUG_MODE_ON) for(auto& i: sc.paths) i.print(index);
			
			for(int i = 0; i < sc.trsts.size(); i++)
			{
				ts_full[0].add(sc.trsts[i], 1, 0, TRANSCRIPT_COUNT_ADD_COVERAGE_MIN, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
			}
			for(int i = 0; i < sc.non_full_trsts.size(); i++)
			{
				ts_nonfull[0].add(sc.non_full_trsts[i], 1, 0, TRANSCRIPT_COUNT_ADD_COVERAGE_MIN, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
			}
			if(!is_allelic || sc.asnonzeroset.size() <= 0) 
			{
				assert(sc.asnonzeroset.size() <= 0);
				assert(sc.nsnonzeroset.size() <= 0);
				continue;
			}
			
			// assemble alleles in seperate splice graphs/ scallops
			phaser ph(sc, is_allelic);				
			vector<transcript>& trsts1 = ph.trsts1;
			vector<transcript>& trsts2 = ph.trsts2;
			vector<transcript>& non_full_trsts1 = ph.non_full_trsts1;
			vector<transcript>& non_full_trsts2 = ph.non_full_trsts2;

			// collect transcripts 
			/*
			if(verbose >= 2)
			{
				printf("assembly with r = %d, total %lu transcripts in ALLELE1:\n", r, trsts1.size());
				for(int i = 0; i < trsts1.size(); i++) trsts1[i].write(cout);
				printf("assembly with r = %d, total %lu transcripts in ALLELE2:\n", r, trsts2.size());
				for(int i = 0; i < sc2.trsts.size(); i++) sc2.trsts[i].write(cout);
			}
			*/

			// add transcripts to corresponding transcript_set
			for(const transcript& _t: trsts1)
			{
				transcript t(_t);
				ts_full[1].add(t, 1, 0, TRANSCRIPT_COUNT_ADD_COVERAGE_MIN, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);				
			}
			for(const transcript& _t: non_full_trsts1)
			{
				transcript t(_t);
				ts_nonfull[1].add(t, 1, 0, TRANSCRIPT_COUNT_ADD_COVERAGE_MIN, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
			}
			for(const transcript& _t: trsts2)
			{
				transcript t(_t);
				ts_full[2].add(t, 1, 0, TRANSCRIPT_COUNT_ADD_COVERAGE_MIN, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
			}
			for(const transcript& _t: non_full_trsts2)
			{
				transcript t(_t);
				ts_nonfull[2].add(t, 1, 0, TRANSCRIPT_COUNT_ADD_COVERAGE_MIN, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
			}

			// also add those to NONSPECIFIC transcript_set; coverage should add for both alleles
			for(const transcript& _t: trsts1)
			{
				transcript t(_t);
				t.make_non_specific();
				ts_full[0].add(t, 1, 0, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);				
			}
			for(const transcript& _t: non_full_trsts1)
			{
				transcript t(_t);
				t.make_non_specific();
				ts_nonfull[0].add(t, 1, 0, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
			}
			for(const transcript& _t:trsts2)
			{
				transcript t(_t);
				t.make_non_specific();
				ts_full[0].add(t, 1, 0, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
			}
			for(const transcript& _t: non_full_trsts2)
			{
				transcript t(_t);
				t.make_non_specific();
				ts_nonfull[0].add(t, 1, 0, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
			}
		}
	}
	return 0;
}

bool assembler::determine_regional_graph(splice_graph &gr)
{
	bool all_regional = true;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.get_vertex_info(i).regional == false) all_regional = false;
		if(all_regional == false) break;
	}
	return all_regional;
}

int assembler::assign_RPKM()
{
	double factor = 1e9 / qlen;
	for(int i = 0; i < trsts.size(); i++)
	{
		vector<transcript>& trsts_of_allele = trsts[i];
		for(int j = 0; j < trsts_of_allele.size(); j++)
		{
			trsts_of_allele[i].assign_RPKM(factor);
		}
	}
	return 0;
}

int assembler::write()
{

	if(verbose >= 1) printf("\nWriting outputs for all transcripts.\n");

	for(int a = 0; a < trsts.size(); a++)
	{
		assert(a >= 0 && a <= 3);
		string allele_name_fix = "merged";
		if (a == 1) allele_name_fix = "allele1";
		if (a == 2) allele_name_fix = "allele2";
		string outname_prefix = output_file + "." + allele_name_fix;

		if(verbose >= 2) 
		{
			if (a == 0) printf("\tWriting outputs for all transcripts of merged alleles.\n");
			if (a == 1) printf("\tWriting outputs for all transcripts of allele 1.\n");
			if (a == 2) printf("\tWriting outputs for all transcripts of allele 2.\n");
			assert(a <= 2);
		}

		vector<transcript> & trsts_of_allele = trsts[a];
		
		// write gtf w/o variants
		if (output_file != "")
		{
			ofstream fout((outname_prefix + ".gtf").c_str());
			if(!fout.fail()) for(const transcript &t : trsts_of_allele) t.write(fout);
			fout.close();
		}
		
		// write gvf w/ variants
		if (output_file != "" || a == 0)
		{
			ofstream gvfout((outname_prefix + ".gvf").c_str());
			if(!gvfout.fail()) 
			{
				gvfout << "#allele \"ALLELE1/2\" correspondes to the first/second allele of \"GT\" field in the vcf file input." << endl;
				for(const transcript &t : trsts_of_allele) t.write_gvf(gvfout);
			}
			gvfout.close();
		}
		
		// write fasta w/ variants
		if(fasta_input != "") 
		{
			ofstream faout((outname_prefix + ".fa").c_str());	
			// if(!faout.fail()) for(const transcript &t : trsts_of_allele) t.write_fasta(faout, 60, fai);
			cerr << "fasta output is not implemented yet." << endl; //TODO:
			faout.close();
		}			

		// write non-full-length gvf w/ variants
		if(output_file1 != "")
		{
			ofstream fout1((output_file1 + "." + allele_name_fix + ".gvf").c_str());
			if(!fout1.fail()) for(const transcript &t : nonfull_trsts[a]) t.write_gvf(fout1);
			fout1.close();
		}
	}
	return 0;
}