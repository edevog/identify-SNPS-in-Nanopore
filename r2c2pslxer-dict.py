#!/usr/bin/env python3

import numpy as np
import io


bases = ["A","C","G","T"]
def checkbases(kmer1,kmer2):
    for b in kmer1+kmer2:
        if b not in bases:
            return False
    return True

def main():
    file = open("ONT_6mer.model")
    ont = [x.split() for x in file.read().splitlines()]
    mismatches = {}
    tcount,qcount = {},{}
    mmpos = {}
    leveldict = {}
    for i in range(1,len(ont)):
        tcount[ont[i][0]] = 0
        qcount[ont[i][0]] = 0
        leveldict[ont[i][0]] = float(ont[i][1])
        for j in range(1,len(ont)):
            kmer1 = ont[i][0]
            kmer2 = ont[j][0]
            same = [kmer1[x]==kmer2[x] for x in range(6)]
            if sum(same) == 5:
                mismatches[kmer1+kmer2]=0
                mmpos[kmer1+kmer2] = np.argmin(same)+1




    with io.open("R2C2.aligned.out.pslx.1",encoding="cp850") as file:
        counter=0
        counter2=0
        counter3=0
        print("pslx...")
        for line in file:
            line = line.split()
            if len(line) == 23:
                seq1,seq2=line[21:]
                query = seq1.split(sep=",")
                target = seq2.split(sep=",")
                if len(query) == len(target):
                    for i in range(len(query)):
                        qilen = len(query[i])
                        tilen = len(target[i])
                        if qilen >= 6 and qilen == tilen:
                            for j in range(qilen-5):
                                kmer1 = query[i][j:j+6].upper()
                                kmer2 = target[i][j:j+6].upper()

                                match = [kmer1[x] == kmer2[x] for x in range(6)]
                                if checkbases(kmer1,kmer2):
                                    qcount[kmer1] += 1
                                    tcount[kmer2] += 1
                                if sum(match) == 5 and checkbases(kmer1,kmer2):
                                    mismatches[kmer1+kmer2] += 1

            counter += 1
            if counter % 200 == 0:
                counter2 += 1
                if counter2 % 200 == 0:
                    print(counter)


    with open("singleMismatchCount.txt","w") as file:
        file.write("kmer1\tkmer2\tdelta\tposition\tcount\n")
        for mm,count in mismatches.items():
            pos = mmpos[mm]
            kmer1,kmer2 = mm[:6],mm[6:]
            diff = leveldict[kmer1]-leveldict[kmer2]
            file.write("\t".join([kmer1,kmer2,str(diff),str(pos),str(count)])+"\n")
    with open("totalKmerCounts.txt","w") as file:
        file.write("kmer\tquery\ttarget\n")
        for kmer,count in qcount.items():
            file.write(kmer+"\t"+str(count)+"\t"+str(tcount[kmer])+"\n")



main()