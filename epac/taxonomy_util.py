#!/usr/bin/env python

def extrac_trainning_taxonomy(fin, fout_test, fout_train):
    with open(fin) as f:
        lines = f.readlines()
    species_list = []
    
    with open(fout_train, "w") as fo:
        with open(fout_test, "w") as fot:
            for line in lines:
                toks = line.strip().split("\t")
                seqid = toks[0]
                ranks = toks[1].split(";")
                spe = ranks[-1].strip()
                if not spe in species_list:
                    species_list.append(spe)
                    fo.write(line)
                else:
                    fot.write(line)
    

def extrac_sequences(fin_taxonomy, fin_db, fout):
    with open(fin_taxonomy) as ft:
        lines = ft.readlines()
    taxlist = []
    for line in lines:
        taxlist.append(line.split("\t")[0].strip())
    
    with open(fin_db) as fd:
        dbs = fd.readlines()
    with open(fout, "w") as fout:
        for tax in taxlist:
            seq = ""
            for i in range(len(dbs)):
                if i % 2 == 0:
                    sid = dbs[i].strip()[1:]
                    if sid == tax:
                        seq = dbs[i + 1]
                        break
            fout.write(">" + tax + "\n")
            fout.write(seq)
        
class TaxonomyUtils:
    EMPTY_RANK = "-"

if __name__ == "__main__":
    #extrac_trainning_taxonomy(fin = "testdata/gg_clostridia_tax.txt", fout_train = "training_tax.txt", fout_test = "testing_tax.txt")
    extrac_sequences(fin_taxonomy = "testing_tax.txt", fin_db = "test1_seq.fa", fout = "testing_seq.fa")
