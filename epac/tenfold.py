#!/usr/bin/env python
import os
import sys
import random
import math

def gentesting(ftaxa, fseq, fout, fold = 10):
    ftax = open(ftaxa)
    lines = ftax.readlines()
    ftax.close()
    
    fse = open(fseq)
    slines = fse.readlines()
    fse.close()
    
    idx = range(len(lines))
    random.seed(12345)
    random.shuffle(idx)
    
    numtaxa = len(lines)
    onefold = int(math.ceil(float(numtaxa) / fold))
    
    idx_list = []
    for i in range(fold):
        start = i * onefold
        end = (i + 1) * onefold
        if end > numtaxa:
            end = numtaxa
        if i == fold -1 :
            end = numtaxa
        idx_list.append(idx[start:end])
    
    for i in range(len(idx_list)):
        idxi = idx_list[i]
        f1 = open(fout + repr(i+1) + "testing.tax", "w")
        f2 = open(fout + repr(i+1) + "testing.fa", "w")
        for index in idxi:
             f1.write(lines[index])
             f2.write(slines[index*2])
             f2.write(slines[index*2 + 1])
        f1.close()
        f2.close()
        
        f1 = open(fout + repr(i+1) + "training.tax", "w")
        f2 = open(fout + repr(i+1) + "training.fa", "w")
        
        for j in range(len(idx_list)):
            if not i==j:
                idxj = idx_list[j]
                for index in idxj:
                    f1.write(lines[index])
                    f2.write(slines[index*2])
                    f2.write(slines[index*2 + 1])
        f1.close()
        f2.close()
        
if __name__ == "__main__":
    if len(sys.argv) < 3: 
        print(python tenfold.py taxonomy fasta fout)
        sys.exit()
    gentesting(ftaxa = sys.argv[1], fseq = sys.argv[2], fout = sys.argv[3], fold = 10)
