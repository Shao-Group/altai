"""
  Part of Altai
  (c) 2021 by Xiaofei Carl Zang, Mingfu Shao, and The Pennsylvania State University.
  See LICENSE for licensing.
"""

class Transcript:
    def __init__(self):
        self.__exons__ = []         # [ExonLike]
        self.__variants__ = []      # [ExonLike]

    def get_exons(self):
        return tuple(self.__exons__)

    def get_intron_chain(self):
        intron_chain_list = []
        for i in range(len(self.__exons__) - 1):
            intron = (self.__exons__[i].p2, self.__exons__][i+1].p1)
            intron_chain_list.append(intron)
        return tuple(intron_chain_list)

    def add_exon(self, left_in, right_ex):
        self.__exons__.append(ExonLike(left_in, right_ex, 0))

    def add_variant(self, left_in, right_ex, seq):
        self.__variants__.append(ExonLike(left_in, right_ex, 1, seq))

    def merge_subexon(self):
        l = []
        last_r = 0
        for i in self.__exons__:
            if i.p1 > last_r:
                l.append(i)
            else:
                l[-1] = ExonLike(l[-1].p1, i.p2, 0)
            last_r = i.p2
        for i in range(len(l) - 1):
            assert l[i].p2 < l[i + 1].p1
        self.__exons__ = l


class ExonLike:
    def __init__(self, left_in, right_ex, feature, seq=''):
        self.p1 = left_in
        self.p2 = right_ex
        self.type = int(feature)  # -1: null, 0: exon, 1: variant
        self.seq = seq

    def get_pos(self):
        return self.p1, self.p2


def read_gff(file_name):
    list_of_transcript = []
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
    return list_of_transcript


if __name__ == "__main__":
    exit()
