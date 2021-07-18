/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Altai
(c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cassert>
#include <sstream>
#include <iostream>
#include <map>
#include <cstring>

#include "config.h"
#include "assembler.h"
#include "scallop.h"
#include "sgraph_compare.h"
#include "super_graph.h"
#include "filter.h"
#include "vcf_data.h"
#include "util.h"

assembler::assembler()
{
    sfn = sam_open(input_file.c_str(), "r");
    hdr = sam_hdr_read(sfn);
    b1t = bam_init1();
	index = 0;
	terminate = false;
	qlen = 0;
	qcnt = 0;
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

		hit ht(b1t, string(buf));
		ht.set_tags(b1t);
		ht.set_strand();
		
		if(verbose >= 4) ht.print();

		
		qlen += ht.qlen;
		qcnt += 1;

		// truncate
		if(ht.tid != bb1.tid || ht.pos > bb1.rpos + min_bundle_gap)
		{
			pool.push_back(bb1);
			bb1.clear();
		}
		if(ht.tid != bb2.tid || ht.pos > bb2.rpos + min_bundle_gap)
		{
			pool.push_back(bb2);
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

	filter ft(trsts);
	// ft.merge_single_exon_transcripts(); // FIXME: 
	ft.keep_as_transcripts_only();
	trsts = ft.trs;
	write();
	
	cout << "Altai finished running." << endl;

	return 0;
}

int assembler::process(int n)
{
	if(pool.size() < n) return 0;

	for(int i = 0; i < pool.size(); i++)
	{
		bundle_base &bb = pool[i];
		bb.buildbase();

		if(verbose >= 3) printf("bundle %d has %lu reads\n", i, bb.hits.size());

		if(bb.hits.size() < min_num_hits_in_bundle) continue;
		if(bb.tid < 0) continue;

		char buf[1024];
		strcpy(buf, hdr->target_name[bb.tid]);

		bundle bd(bb);

		bd.chrm = string(buf);
		bd.build();
		if(verbose >= 1) bd.print(index);
		
		assemble(bd.gr, bd.hs, bd.is_allelic);
		index++;
	}
	pool.clear();
	return 0;
}

int assembler::assemble(const splice_graph &gr0, const hyper_set &hs0, bool is_allelic)
{
	super_graph sg(gr0, hs0);
	sg.build();

	vector<transcript> gv;
	for(int k = 0; k < sg.subs.size(); k++)
	{
		string gid = "gene." + tostring(index) + "." + tostring(k);
		if(fixed_gene_name != "" && gid != fixed_gene_name) continue;

		if(verbose >= 2 && (k == 0 || fixed_gene_name != "")) sg.print();

		splice_graph &gr = sg.subs[k];
		hyper_set &hs = sg.hss[k];

		if(determine_regional_graph(gr) == true) continue;
		if(gr.num_edges() <= 0) continue;

		gr.gid = gid;
		scallop sc(gr, hs);
		sc.assemble(is_allelic);
		if(verbose >=3) for(auto i: sc.paths) i.print(index);

		if(verbose >= 2)
		{
			printf("transcripts:\n");
			for(int i = 0; i < sc.trsts.size(); i++) sc.trsts[i].write(cout);
		}

		filter ft(sc.trsts);
		if (FILTER_BY_COV) ft.simple_phase_set_by_coverage();
		else ft.simple_phase_set_by_variant_number();
		// ft.join_single_exon_transcripts(); //FIXME:
		ft.filter_length_coverage(); //FIXME:
		if(ft.trs.size() >= 1) gv.insert(gv.end(), ft.trs.begin(), ft.trs.end());

		if(verbose >= 2)
		{
			printf("transcripts after filtering:\n");
			for(int i = 0; i < ft.trs.size(); i++) ft.trs[i].write(cout);
		}

		if(fixed_gene_name != "" && gid == fixed_gene_name) terminate = true;
		if(terminate == true) return 0;
	}

	filter ft(gv);
	// ft.remove_nested_transcripts(); //FIXME:
	if(ft.trs.size() >= 1) trsts.insert(trsts.end(), ft.trs.begin(), ft.trs.end());

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
		trsts[i].assign_RPKM(factor);
	}
	return 0;
}

int assembler::write()
{
	ofstream fout((output_file+".gtf").c_str());
	ofstream gvfout((output_file+".gvf").c_str());
	ofstream faout((output_file+".fa").c_str());
	// ofstream asout((output_file+".ASOnly.fa").c_str());
	if(fout.fail()) return 0;
	if(faout.fail()) return 0;
	for(int i = 0; i < trsts.size(); i++)
	{
		transcript &t = trsts[i];
		t.write(fout);
		t.write_gvf(gvfout);
		if(fasta_input != "") t.write_fasta(faout, 60, fai);
		// if(fasta_input != "") t.write_fasta_AS_only(asout, 60, fai);
	}
	fout.close();
	faout.close();
	gvfout.close();
	// asout.close();
	return 0;
}