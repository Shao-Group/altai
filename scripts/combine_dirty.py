# gtf 1
# gtf 2
# fa 1
# fa 2
import sys

def read_gff(file_name):
    with open(file_name) as f:
        lines = f.readlines()

    dict = {} # trst_name:[tuples]
    for i in lines:
        e = i.strip().split('\t')
        # get attributes
        attr = {}
        print(e[8])
        for x in e[8].split(";"):
            if x != "":
                print(x)
                x = x.strip()
                y = x.split("\"")
                k = y[0].strip()
                v = y[1].strip()
                assert (y[2].strip() == "")
                attr[k] = v

        trst_id = attr["transcript_id"]
        if trst_id not in dict:
            dict[trst_id] = []
        if e[2] == "exon" or e[2] == "variant":
            dict[trst_id].append((int(e[3]), int(e[4])))

    # merge overlapping exons
    for k in dict:
        es = dict[k]
        l = []
        last_r = -1
        for i in es:
            if i[0] > last_r:
                l.append(i)
            else:
                l[-1] = (l[-1][0], i[1])
            last_r = i[1]
        for i in range(len(l) - 1):
            assert l[i][1] < l[i + 1][0]
        dict[k] = l

    return dict

# retain d2 - d1
def retain(d1, d2):
    s = set(d1.values())
    rt = []
    for k, v in d2.items():
        if v not in s:
            rt.append(k)
    return rt


if __name__ == "__main__":
    if len(sys.argv) != 3:
        exit()
    _dummy, g1, g2 = sys.argv
    d1 = read_gff(g1)
    d2 = read_gff(g2)
    rt = retain(d1, d2)
    for i in rt:
        print(i)