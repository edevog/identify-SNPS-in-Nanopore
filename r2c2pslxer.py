#!/usr/bin/env python3

import numpy as np
import io


bases = ["a","c","g","t"]
def checkbases(kmer1,kmer2):
    for b in kmer1+kmer2:
        if b not in bases:
            return False
    return True

def main():
    with io.open("R2C2.aligned.out.pslx.1",encoding="cp850") as file:
        mismatches = []
        nummismatches = [0,0,0,0,0,0]
        counter=0
        counter2=0
        counter3=0
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
                                kmer1 = query[i][j:j+6]
                                kmer2 = target[i][j:j+6]
                                match = [kmer1[x] == kmer2[x] for x in range(6)]
                                if sum(match) == 5 and checkbases(kmer1,kmer2):
                                    mismatches.append([kmer1,kmer2,np.argmin(match)+1])
            counter += 1
            if counter % 200 == 0:
                counter2 += 1
                if counter2 % 200 == 0:
                    print(counter,len(mismatches))
                    break
    with open("ONT_6mer.model") as file:
        mmlevels = [x.split() for x in file.read().splitlines()]
        leveldict = {}
        for line in mmlevels[1:]:
            leveldict[line[0]] = float(line[1])
    with open("diffoutnew.txt","w") as file:
        for mm in mismatches:
            kmer1,kmer2 = mm[0],mm[1]
            diff = leveldict[kmer1.upper()]-leveldict[kmer2.upper()]
            file.write(str(diff)+"\t"+str(mm[2])+"\n")


main()