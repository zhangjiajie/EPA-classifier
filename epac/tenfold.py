#!/usr/bin/env python
import os
import sys
import random
import math
from ete2 import SeqGroup

def gentesting(ftaxa, fseq, fout, fold = 10):
    ftax = open(ftaxa)
    lines = ftax.readlines()
    ftax.close()
    
    seqs = SeqGroup(fseq)
    
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
             tax = lines[index]
             seqid = tax.split()[0]
             seq = seqs.get_seq(seqid)
             f1.write(tax)
             f2.write(">" + seqid + "\n")
             f2.write(seq + "\n")
        f1.close()
        f2.close()
        
        f1 = open(fout + repr(i+1) + "training.tax", "w")
        f2 = open(fout + repr(i+1) + "training.fa", "w")
        
        for j in range(len(idx_list)):
            if not i==j:
                idxj = idx_list[j]
                for index in idxj:
                    tax = lines[index]
                    seqid = tax.split()[0]
                    seq = seqs.get_seq(seqid)
                    f1.write(tax)
                    f2.write(">" + seqid + "\n")
                    f2.write(seq + "\n")
        f1.close()
        f2.close()
        
if __name__ == "__main__":
    if len(sys.argv) < 3: 
        print("python tenfold.py taxonomy fasta fout")
        sys.exit()
    gentesting(ftaxa = sys.argv[1], fseq = sys.argv[2], fout = sys.argv[3], fold = 10)
