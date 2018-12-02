#!/usr/bin/env python3
import sys
import random

'''
This program creates a subset of the dataset from the input data file by
randomly selecting half of the data points from the initial dataset

Input:
A text file of tab delimited data in the following example format:

kmer1	kmer2	delta	position	count	proportion
AAAAAA	AAAAAC	2.5374979999999994	6	16288	0.009487716878326788
AAAAAA	AAAAAG	1.0109679999999912	6	23350	0.013601313182031588
AAAAAA	AAAAAT	2.0624289999999945	6	9423	0.005488872553074247
AAAAAA	AAAACA	9.389054999999999	5	18501	0.010776783519518903

Output:
Text file with the number of data points and the list of data points with
the following data format:

Number of data points: 36864
kmer1	kmer2	delta	position	count	proportion
CTTATA	CTTCTA	-7.54505300000001	4	79	0.0006676695796217103
GATGGC	GATGGG	-1.5108359999999976	6	388	0.000697929595452665
CACACA	CACAGA	-5.464443000000003	5	1080	0.0016235084016559785
'''

class DataWrangling(object):
    """
    Attributes:
        kmer_dict = a dictionary with the keys as the categories kmer1, kmer2,
        delta, position, count, and proportion where kmer1 is the expected
        kmer, kmer2 is miscalled kmer, delta is the difference in signal
        between kmer1 and kmer2, position is position in the hexamer where the
        miscalled base occurs, count is the number of times
        the kmer pair occurs, proportion is the proportion of number of times
        the kmer pair occurs to the total number of miscalled bases
    """
    def __init__(self, kmer_dict):
        self.kmer_dict = kmer_dict

    def selectRandomData(self):
        '''
        randomly selects half of the dataset
        '''
        random_dict = {k : [] for k in self.kmer_dict}
        data_indices = random.sample(list(range(len(self.kmer_dict["kmer1"]))), int(len(self.kmer_dict["kmer1"])/2))#len(self.kmer_dict["kmer1"])/2)
        for i in data_indices:
            for k in random_dict:
                random_dict[k].append(self.kmer_dict[k][i])
        return random_dict



def main():
    keys = sys.stdin.readline().strip().split("\t")
    data = sys.stdin.read().splitlines()
    kmer_dict = {k : [] for k in keys}
    for d in data:
        d = d.split('\t')
        for i in range(len(d)):
            kmer_dict[keys[i]].append(d[i])

    dw = DataWrangling(kmer_dict)
    random_data  = dw.selectRandomData()
    #prints the dictionary in a readable format
    print("Number of data points: " + str(len(random_data["kmer1"])))
    print(*random_data, sep = "\t")
    for i in range(len(random_data["kmer1"])):
        t_values = []
        for k in kmer_dict:
            t_values.append(random_data[k][i])
        print(*t_values, sep="\t")



if __name__ == '__main__':
    main()
