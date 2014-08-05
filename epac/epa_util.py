#! /usr/bin/env python
import sys
import os
import json
import operator
import time
from epac.ete2 import SeqGroup, Tree
from subprocess import call
from json_util import jsonparser


class epa:
    def __init__(self, raxmlbin = "raxmlHPC-PTHREADS-SSE3"):
        self.basepath = os.path.dirname(os.path.abspath(__file__))
        self.raxmlpath = self.basepath + "/bin/" + raxmlbin
        self.tmppath = self.basepath + "/tmp"
        if not os.path.exists(self.raxmlpath):
            print("The pipeline uses RAxML to infer phylogenetic trees,")
            print("please download the latest source code from: ")
            print("https://github.com/stamatak/standard-RAxML")
            print("Please complie the SSE + PTHREAD version, ")
            print("rename the executable to raxmlHPC-PTHREADS-SSE3 and put it to bin/  \n")
            sys.exit() 
        if not os.path.exists(self.tmppath):
            print("The tmp folder for keeping intermediate files does not exit")
            print("Please create this folder:" + self.tmppath)
            sys.exit() 
        self.name = str(time.time())
    
    
    def run(self, reftree, alignment, num_thread = "2", model = None ):
        #./raxmlHPC-PTHREADS-SSE3 -f v -m GTRGAMMA -s ../example/t1.fa -t ../example/t1.test.tre -p 1234 -n jjj -T 2
        if os.path.exists(reftree):
            if model == None:
                call([self.raxmlpath,"-f", "v", "-m","GTRGAMMA","-s",alignment,"-t", reftree, "-n",self.name,"-p", "1234", "-T", str(num_thread), "-w", self.tmppath] , stdout=open(os.devnull, "w"), stderr=open(os.devnull, "w"))
            else:
                call([self.raxmlpath,"-f", "v", "-G", "0.01", "-R", model, "-m","GTRGAMMA","-s",alignment,"-t", reftree, "-n",self.name,"-p", "1234", "-T", str(num_thread), "-w", self.tmppath] , stdout=open(os.devnull, "w"), stderr=open(os.devnull, "w"))
        else:
            tmptree = self.tmppath+"/" + self.name + ".tre"
            with open(tmptree, "w") as fout:
                fout.write(reftree)
            if model == None:
                call([self.raxmlpath,"-f", "v", "-m","GTRGAMMA","-s",alignment,"-t", tmptree, "-n",self.name,"-p", "1234", "-T", str(num_thread), "-w", self.tmppath] , stdout=open(os.devnull, "w"), stderr=open(os.devnull, "w"))
            else:
                call([self.raxmlpath,"-f", "v", "-G", "0.01","-R", model, "-m","GTRGAMMA","-s",alignment,"-t", tmptree, "-n",self.name,"-p", "1234", "-T", str(num_thread), "-w", self.tmppath] , stdout=open(os.devnull, "w"), stderr=open(os.devnull, "w"))

        jp = jsonparser(self.tmppath + "/" + "RAxML_portableTree." + self.name + ".jplace")
        return jp
    
    
    def clean(self):
        os.remove(self.tmppath + "/" + "RAxML_classification." + self.name)
        os.remove(self.tmppath + "/" + "RAxML_classificationLikelihoodWeights." + self.name)
        os.remove(self.tmppath + "/" + "RAxML_entropy." + self.name)
        os.remove(self.tmppath + "/" + "RAxML_info." + self.name)
        os.remove(self.tmppath + "/" + "RAxML_labelledTree." + self.name)
        os.remove(self.tmppath + "/" + "RAxML_originalLabelledTree." + self.name)
        #os.remove(self.tmppath + "/" + "RAxML_portableTree." + self.name + ".jplace")
        os.remove(self.tmppath + "/" + self.name + ".tre")
        
     
    def dummy(self, reftree, alignment):
        seqs = SeqGroup(sequences=alignment, format='fasta')
        entries = seqs.get_entries()
        seq0 = entries[0][1]
        dummyseq = seq0[:-50] + "A"*50
        seqs.set_seq(name = "dummy", seq = dummyseq)
        fout = self.tmppath + "/dummy" + self.name + ".fa"
        seqs.write(format='fasta', outfile=fout) 
        self.run(reftree = reftree, alignment = fout)
        self.clean()
        os.remove(fout)
        return self.tmppath + "/" + "RAxML_portableTree." + self.name + ".jplace"


class raxml:
    def __init__(self, raxmlbin = "raxmlHPC-PTHREADS-SSE3"):
        self.basepath = os.path.dirname(os.path.abspath(__file__))
        self.raxmlpath = self.basepath + "/bin/" + raxmlbin
        self.tmppath = self.basepath + "/tmp"
        if not os.path.exists(self.raxmlpath):
            print("The pipeline uses RAxML to infer phylogenetic trees,")
            print("please download the latest source code from: ")
            print("https://github.com/stamatak/standard-RAxML")
            print("Please complie the SSE + PTHREAD version, ")
            print("rename the executable to raxmlHPC-PTHREADS-SSE3 and put it to bin/  \n")
            sys.exit() 
        if not os.path.exists(self.tmppath):
            print("The tmp folder for keeping intermediate files does not exit")
            print("Please create this folder:" + self.tmppath)
            sys.exit() 
        self.name = str(time.time())
    
    def resolve_mftree(self, mftree, alignment, num_thread = "2"):
        if os.path.exists(mftree):
            call([self.raxmlpath, "-m","GTRGAMMA","-s",alignment,"-g", mftree, "-n",self.name,"-p", "1234", "-T", num_thread, "-w", self.tmppath] ) #, stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
        else:
            tmptree = self.tmppath+"/" + self.name + ".tre"
            with open(tmptree, "w") as fout:
                fout.write(mftree)
            call([self.raxmlpath,"-m","GTRGAMMA","-s",alignment,"-g", tmptree, "-n",self.name,"-p", "1234", "-T", num_thread, "-w", self.tmppath] ) #, stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
        return self.tmppath + "/" + "RAxML_bestTree." + self.name
    
    def get_model_parameters(self, bftree, alignment, num_thread = "2"):
        #./raxmlHPC-AVX -f e -s 100.phy -m GTRGAMMA -t RAxML_bestTree.T1 -n Evaluate
        #RAxML_binaryModelParameters.Evaluate
        if os.path.exists(bftree):
            call([self.raxmlpath, "-f", "e", "-m","GTRGAMMA","-t", bftree, "-s",alignment, "-n",self.name,"-p", "1234", "-T", num_thread, "-w", self.tmppath] ) #, stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
        else:
            tmptree = self.tmppath+"/" + self.name + ".tre"
            with open(tmptree, "w") as fout:
                fout.write(bftree)
            call([self.raxmlpath, "-f", "e", "-m","GTRGAMMA","-t", bftree, "-s",alignment, "-n",self.name,"-p", "1234", "-T", num_thread, "-w", self.tmppath] ) #, stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
        return self.tmppath + "/" + "RAxML_binaryModelParameters." + self.name

        
    def clean(self):
        os.remove(self.tmppath + "/" + "RAxML_info." + self.name)
        os.remove(self.tmppath + "/" + "RAxML_log." + self.name)
        os.remove(self.tmppath + "/" + "RAxML_result." + self.name)


def epa_2_ptp(placements, min_lw = 0.5):
    placemap = {}
    """find how many edges are used for placement, and create a map to store """
    for placement in placements:
        edges = placement["p"]
        curredge = edges[0][0]
        lw = edges[0][2] 
        if lw >= min_lw:
            placemap[curredge] = placemap.get(curredge, [])

    """group taxa name by placement branch"""
    for placement in placements:
        edges = placement["p"]
        taxa_names = placement["n"]
        curredge = edges[0][0]
        lw = edges[0][2] 
        if lw >= min_lw:
            a = placemap[curredge] 
            a.extend(taxa_names)
            placemap[curredge]  = a

    

if __name__ == "__main__":
    if len(sys.argv) < 3: 
        print("usage: ./epa_util.py <multifurcating.tre> <alignment> <num_thread>")
        sys.exit()
    rx = raxml()
    tree = rx.resolve_mftree(sys.argv[1], sys.argv[2], sys.argv[3])
    rx.clean()
    print("The resolved tree has been written to: ")
    print(tree)
