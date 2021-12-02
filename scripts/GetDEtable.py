# gtf 1
# gtf 2
# fa 1
# fa 2
import sys

def read_gff(file_name):
    with open(file_name) as f:
        lines = f.readlines()

    diction = {} # trst_name:[tuples]
    for i in lines:
        e = i.strip().split('\t')
        # get attributes
        attr = {}
#        print(e[8])
        for x in e[8].split(";"):
            if x != "":
#                print(x)
                x = x.strip()
                y = x.split("\"")
                k = y[0].strip()
                v = y[1].strip()
                assert (y[2].strip() == "")
                attr[k] = v

        trst_id = attr["transcript_id"]
        if trst_id not in diction:
            diction[trst_id] = []
        if e[2] == "exon" or e[2] == "variant":
            diction[trst_id].append((int(e[3]), int(e[4])))

    # merge overlapping exons
    for k in diction:
        es = diction[k]
        l = []
        last_r = -1
        for i in es:
            if i[0] > last_r + 1:
                l.append(i)
            else:
                l[-1] = (l[-1][0], i[1])
            last_r = i[1]
        for i in range(len(l) - 1):
            assert l[i][1] < l[i + 1][0]
        diction[k] = tuple(l)

    dd = {}
    for k,v in diction.items():
        # intron chain
        l = list(v)
        l[0] = tuple([l[0][1]])
        l[-1] = tuple([l[-1][0]])
        v = tuple(l)
        if v not in dd:
            dd[v] = [k]
        else:
            dd[v].append(k)
    return dd # tuple(intron_chain): [names]

def obtain_counts(list_of_d, list_of_counttables, ids):
    gene_set = set()
    for i in list_of_d:
        gene_set = gene_set.union(set(i.keys()))  #TODO: also by chr number

    dict_high = {}
    dict_low = {}
    gene_enu = 0
    for ii in gene_set:
        gene_enu += 1
        dict_high ['gene' + str(gene_enu)] = {}
        dict_low ['gene' + str(gene_enu)] = {}
        high_value = -1
        low_value = -1
        for j in range(len(list_of_d)):
            try:
                names = list_of_d[j][ii]
                if len(names) < 2:
                    continue
                counts = []
                for n in names:
                    try:
                        counts.append(float(list_of_counttables[j]['altai_'+ n]))
                    except KeyError:
                        pass
                #sort two lists together
                if len(names) == len(counts):
                    counts, names = zip(*sorted(zip(counts, names)))
                    dict_high['gene' + str(gene_enu)][j] = (counts[-1], names[-1])
                    dict_low['gene' + str(gene_enu)][j] = (counts[0], names[0])
                else:
                    continue
            except KeyError:
                pass

    # write out
    file_count = open('DE_table.counts.tsv', 'w')
    file_gname = open('DE_table.gname.tsv', 'w')
    h_header = [x+'_high' for x in ids]
    l_header = [x+'_low' for x in ids]
    header = 'gene_id\t' + '\t'.join(['\t'.join(x) for x in zip(h_header, l_header)])
    file_count.write(header + '\n')
    file_gname.write(header + '\n')

    for i in dict_high.keys():
        line_c = i + '\t' # gene id
        line_g = i + '\t'
        for j in range(len(list_of_d)):
            try:
                c, n = dict_high[i][j]
            except KeyError:
                c, n = ('na', 'na')
            line_c += str(c) + '\t'
            line_g += str(n) + '\t'
            try:
                c, n = dict_low[i][j]
            except KeyError:
                c, n = ('na', 'na')
            line_c += str(c) + '\t'
            line_g += str(n) + '\t'
        file_count.write(line_c + '\n')
        file_gname.write(line_g + '\n')
    file_count.close()
    file_gname.close()


if __name__ == "__main__":
    # read id
    ids_file = sys.argv[1]
    with open(ids_file, 'r') as f:
        ids = [ x.strip() for x in f.readlines()]

    d = []
    # only considered altai gvf
    c = []
    for id in ids:
        cc = {}
        d.append(read_gff('/gpfs/group/restricted/mxs2589/dbgap/altai-gtex/altai/' + id + '.altai.gvf'))
        with open('/gpfs/group/restricted/mxs2589/dbgap/altai-gtex/kallisto/' + id + '.combined.quant/abundance.tsv') as f:
            ll = f.readlines()[1:]
            for i in ll:
                j = i.split()
                cc[j[0]] = float(j[3])
        c.append(cc)
    assert len(d) == len(c)
    obtain_counts(d, c, ids)



