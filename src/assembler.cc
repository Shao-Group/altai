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
#include "specific_trsts.hpp"

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
	specific_full_trsts.resize(3);
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
	/* for(int i = 0; i < 3; i++)
	{
		if (! use_filter) break;

		filter ft(trsts[i]);
		ft.merge_single_exon_transcripts();
		trsts[i] = ft.trs;

		filter ft1(nonfull_trsts[i]);
		ft1.merge_single_exon_transcripts();
		nonfull_trsts[i] = ft1.trs;
	}	 */

	// get specific trsts
	specific_full_trsts.clear();

	for(transcript t: trsts[1])
	{
		if(t.gt == ALLELE1) specific_full_trsts[1].push_back(t);
		else if(t.gt == NONSPECIFIC || t.gt == UNPHASED) specific_full_trsts[0].push_back(t);
	}
	for(transcript t: trsts[2])
	{
		if(t.gt == ALLELE2) specific_full_trsts[2].push_back(t); 			
		else if(t.gt == NONSPECIFIC || t.gt == UNPHASED) specific_full_trsts[0].push_back(t);
	}

	for(transcript& t: specific_full_trsts[0]) t.make_non_specific();
	filter ft0(specific_full_trsts[0]);
	/* ft0.merge_single_exon_transcripts();
	ft0.filter_length_coverage();
	ft0.remove_nested_transcripts(); */
	specific_full_trsts[0].clear();
	specific_full_trsts[0] = ft0.trs;

	//TODO: nf trsts
	// if(recover_partial_tx_min_overlap_with_full_tx > 0)
	// {
	// 	double f = recover_partial_tx_min_overlap_with_full_tx;
	// 	recovered_allele1 = specific_trsts::recover_full_from_partial_transcripts(trsts[0], nonfull_trsts[1], f, true);
	// 	recovered_allele2 = specific_trsts::recover_full_from_partial_transcripts(trsts[0], nonfull_trsts[2], f, true);
	// }

	trsts[0] = trsts_collective;

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
		vector<transcript_set> ts_full;		 	// full-length set; [0]merged, [1]ALLELE1, [2]ALLELE2
		vector<transcript_set> ts_nonfull;		// non-full-length set; [0]merged, [1]ALLELE1, [2]ALLELE2
		for(int i = 0; i < 3; i++) ts_full.push_back(transcript_set(bb.chrm, 0.9));
		for(int i = 0; i < 3; i++) ts_nonfull.push_back(transcript_set(bb.chrm, 0.9));

		bundle bd(bb);
		if(bundle_mode == 1 || bundle_mode == 3)
		{
			bd.build(1, true);
			bd.print(index++);
			assemble(bd.gr, bd.hs, bb.is_allelic, ts_full, ts_nonfull);
		}
		if(bundle_mode == 2 || bundle_mode == 3)
		{
			bd.build(2, true);
			bd.print(index++);				
			assemble(bd.gr, bd.hs, bb.is_allelic, ts_full, ts_nonfull);
		}
		
		// get allele spec transcripts again
		//TODO: deal with nonfull
		// int sdup = assemble_duplicates / 1 + 1;
		// int mdup = assemble_duplicates / 2 + 0;
		int sdup = 0;
		int mdup = 0;
		vector<transcript> tx0 = ts_full[0].get_transcripts(sdup, mdup); 
		vector<transcript> tx1 = ts_full[1].get_transcripts(sdup, mdup); 
		vector<transcript> tx2 = ts_full[2].get_transcripts(sdup, mdup); 
		specific_trsts::get_allele_spec_trsts(tx1, tx2, min_allele_transcript_cov);


		// retrieve and filter transcripts
		// i = {0, 1, 2}, corresponds to merged, ALLELE1, ALLELE2
		for(int i = 0; i <= 2; i++)
		{
			genotype gg = UNPHASED;
			if (DEBUG_MODE_ON)	assert(trsts.size() == 3 && nonfull_trsts.size() == 3 && specific_full_trsts.size() == 3);
			if (i == 0) gg = NONSPECIFIC;
			if (i == 1) gg = ALLELE1;
			if (i == 2) gg = ALLELE2;
			
			vector<transcript>* gv1;
			if(i == 0) gv1 = &tx0;
			else if(i == 1) gv1 = &tx1;
			else if(i == 2) gv1 = &tx2;
			else assert(0);

			for(int k = 0; k < gv1->size(); k++)
			{
				if((gv1->at(k)).exons.size() >= 2) gv1->at(k).coverage /= (1.0 * assemble_duplicates);
				if((i == 1 || i == 2) && DEBUG_MODE_ON) assert(!(gt_conflict(gv1->at(k).gt, gg)));
			}

			/*
			if (use_filter)
			{
				filter ft1(*gv1);
				ft1.filter_length_coverage();
				ft1.remove_nested_transcripts();
				if(ft1.trs.size() >= 1) trsts[i].insert(trsts[i].end(), ft1.trs.begin(), ft1.trs.end());

				// filter ft2(gv2);
				// ft2.filter_length_coverage();
				// ft2.remove_nested_transcripts();
				// if(ft2.trs.size() >= 1) nonfull_trsts[i].insert(nonfull_trsts[i].end(), ft2.trs.begin(), ft2.trs.end());
			}
			else
			*/
			{
				trsts[i].insert(trsts[i].end(), gv1->begin(), gv1->end());
				// nonfull_trsts[i].insert(nonfull_trsts[i].end(), gv2.begin(), gv2.end());
			}
		}		
	}
	pool.clear();
	return 0;
}

int assembler::assemble(const splice_graph &gr0, const hyper_set &hs0, bool is_allelic, vector<transcript_set> &ts_full, vector<transcript_set> &ts_nonfull)
{
	string chrm = gr0.chrm;
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
			// 0: merged; 1: ALLELE1; 2: ALLELE2 
			transcript_set fl_add_0(chrm, 0.9);		// full length allele 0, mode1 = mode2 = TRANSCRIPT_COUNT_ONE_COVERAGE_ADD := count is always 1, cov +=
			transcript_set fl_add_1(chrm, 0.9);		// full length allele 1, mode1 = mode2 = TRANSCRIPT_COUNT_ONE_COVERAGE_ADD
			transcript_set fl_add_2(chrm, 0.9);		// full length allele 2, mode1 = mode2 = TRANSCRIPT_COUNT_ONE_COVERAGE_ADD

			transcript_set nf_add_0(chrm, 0.9);		// nonfull len allele 0, mode1 = mode2 = TRANSCRIPT_COUNT_ONE_COVERAGE_ADD
			transcript_set nf_add_1(chrm, 0.9);		// nonfull len allele 1, mode1 = mode2 = TRANSCRIPT_COUNT_ONE_COVERAGE_ADD
			transcript_set nf_add_2(chrm, 0.9);		// nonfull len allele 2, mode1 = mode2 = TRANSCRIPT_COUNT_ONE_COVERAGE_ADD


			splice_graph gr = gr_alias;  // graph copy for different duplicates, copy constructor used, not move
			string gid = "gene." + tostring(index) + "." + tostring(k) + "." + tostring(r);
			gr.gid = gid;

			// decompose a merged graph
			splice_graph gr_copy(gr);
			hyper_set hs_copy(hs);
			scallop sc0(gr_copy, hs_copy, false, false);
			sc0.assemble(false);
			for(const transcript& _t: sc0.trsts)
			{
				transcript t0(_t);
				t0.make_non_specific();
				trsts_collective.push_back(t0);
			}

			// partial decomp of non-AS nodes
			scallop sc(gr, hs, r == 0 ? false : true, true);
			sc.assemble(is_allelic);
			if(verbose >= 2)
			{
				printf("assembly with r = %d; %lu transcripts in partial decomposition of merged splice graph\n", r, sc.trsts.size());
			}
		
			for(const transcript& _t: sc.trsts)
			{
				trsts_collective.push_back(_t);
				fl_add_0.add(transcript(_t), 1, 0, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD);
				if(ALLELE2 != _t.gt) fl_add_1.add(transcript(_t), 1, 0, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD);
				if(ALLELE1 != _t.gt) fl_add_2.add(transcript(_t), 1, 0, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD);
			}
			for(const transcript& _t: sc.non_full_trsts)
			{
				nf_add_0.add(transcript(_t), 1, 0, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD);
				if(ALLELE2 != _t.gt) nf_add_1.add(transcript(_t), 1, 0, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD);
				if(ALLELE1 != _t.gt) nf_add_2.add(transcript(_t), 1, 0, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD);
			}

			// assemble alleles in seperate splice graphs/ scallops
			phaser ph(sc, is_allelic);				
			vector<transcript>& trsts1 = ph.trsts1;
			vector<transcript>& trsts2 = ph.trsts2;
			vector<transcript>& non_full_trsts1 = ph.non_full_trsts1;
			vector<transcript>& non_full_trsts2 = ph.non_full_trsts2;

			// collect transcripts 
			if(verbose >= 2)
			{
				printf("assembly with r = %d; another %lu transcripts in splice graph of ALLELE1\n", r, trsts1.size());
				printf("assembly with r = %d; another %lu transcripts in splice graph of ALLELE2\n", r, trsts2.size());
			}

			// add transcripts to corresponding transcript_set
			for(const transcript& _t: trsts1)
			{
				transcript t0(_t);
				t0.make_non_specific();
				trsts_collective.push_back(t0);
				fl_add_0.add(t0, 1, 0, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD);
				fl_add_1.add(transcript(_t), 1, 0, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD);
			}
			for(const transcript& _t: non_full_trsts1)
			{
				transcript t0(_t);
				t0.make_non_specific();
				nf_add_0.add(t0, 1, 0, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD);
				nf_add_1.add(transcript(_t), 1, 0, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD);
			}
			for(const transcript& _t: trsts2)
			{
				transcript t0(_t);
				t0.make_non_specific();
				trsts_collective.push_back(t0);
				fl_add_0.add(t0, 1, 0, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD);
				fl_add_2.add(transcript(_t), 1, 0, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD);
			}
			for(const transcript& _t: non_full_trsts2)
			{
				transcript t0(_t);
				t0.make_non_specific();
				nf_add_0.add(t0, 1, 0, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD);
				nf_add_2.add(transcript(_t), 1, 0, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD, TRANSCRIPT_COUNT_ONE_COVERAGE_ADD);
			}

			ts_full[0].add(fl_add_0, TRANSCRIPT_COUNT_ADD_COVERAGE_MIN, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
			ts_full[1].add(fl_add_1, TRANSCRIPT_COUNT_ADD_COVERAGE_MIN, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
			ts_full[2].add(fl_add_2, TRANSCRIPT_COUNT_ADD_COVERAGE_MIN, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
			ts_nonfull[0].add(nf_add_0, TRANSCRIPT_COUNT_ADD_COVERAGE_MIN, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
			ts_nonfull[1].add(nf_add_1, TRANSCRIPT_COUNT_ADD_COVERAGE_MIN, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
			ts_nonfull[2].add(nf_add_2, TRANSCRIPT_COUNT_ADD_COVERAGE_MIN, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
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
		string allele_name_fix;
		if (a == 0) allele_name_fix = "merged";
		if (a == 1) allele_name_fix = "allele1";
		if (a == 2) allele_name_fix = "allele2";
		string outname_prefix = output_file + "." + allele_name_fix;

		if(verbose >= 2) 
		{
			if (a == 0) printf("\tWriting outputs for all transcripts of merged/non-specific alleles.\n");
			if (a == 1) printf("\tWriting outputs for all transcripts of allele 1.\n");
			if (a == 2) printf("\tWriting outputs for all transcripts of allele 2.\n");
			assert(a <= 2);
		}
		
		// write gtf
		if (output_file != "")
		{
			ofstream fout((outname_prefix + ".gtf").c_str());
			if(!fout.fail()) for(const transcript &t : trsts[a]) t.write(fout);
			fout.close();
		}
		
		// write gvf w/ variants
		if (output_file != "")
		{
			ofstream gvfout((outname_prefix + ".gvf").c_str());
			if(!gvfout.fail()) 
			{
				gvfout << "#allele \"ALLELE1/2\" correspondes to the first/second allele of \"GT\" field in the vcf file input." << endl;
				for(const transcript &t : trsts[a]) t.write_gvf(gvfout);
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

		// write specific transcritps' gtf and gvf
		if (output_file != "")
		{
			if (a == 0) outname_prefix = output_file + "." + "nonspec.multi-exon";
			if (a == 1) outname_prefix = output_file + "." + "allele1spec.multi-exon";
			if (a == 2) outname_prefix = output_file + "." + "allele2spec.multi-exon";
			ofstream gvfout((outname_prefix + ".gvf").c_str());
			if(!gvfout.fail())	for(const transcript &t : specific_full_trsts[a]) t.write_gvf(gvfout);
			gvfout.close();
			ofstream gtfout((outname_prefix + ".gtf").c_str());
			if(!gtfout.fail())	for(const transcript &t : specific_full_trsts[a]) t.write(gtfout);
			gtfout.close();
		}

		// write non-full-length gvf w/ variants
		if(output_file1 != "")
		{
			ofstream fout1((output_file1 + "." + allele_name_fix + ".gvf").c_str());
			if(!fout1.fail()) for(const transcript &t : nonfull_trsts[a]) t.write_gvf(fout1);
			fout1.close();
		}

		// write recovered partial transcripts
		if(recover_partial_tx_min_overlap_with_full_tx > 0 && a != 0)
		{
			ofstream frecov((output_file1 + ".re." + allele_name_fix + ".gvf").c_str());
			const vector<transcript>& v = (a == 1)? recovered_allele1 : recovered_allele2;
			if(!frecov.fail()) for(const transcript &t : v) t.write_gvf(frecov);
			frecov.close();
		}
	}
	return 0;
}