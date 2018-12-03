#!/usr/bin/env python3

import numpy as np
import io

with open("singleMismatchCount.txt") as file:
    mmCount = [x.split() for x in file.read().splitlines()]


with open("totalKmerCounts.txt") as file:
    qtCount = [x.split() for x in file.read().splitlines()]


tcount,qcount = {},{}
for kmer,query,target in qtCount[1:]:
    tcount[kmer] = float(target)
    qcount[kmer] = float(query)

withProps = []
for kmer1,kmer2,delta,position,count in mmCount[1:]:
    mmProp = float(count) / tcount[kmer1]
    withProps += [[kmer1,kmer2,delta,position,count,str(mmProp)]]

with open("mmProps.txt","w") as file:
    file.write("kmer1\tkmer2\tdelta\tposition\tcount\tproportion\n")
    for mm in withProps:
        file.write("\t".join(mm)+"\n")