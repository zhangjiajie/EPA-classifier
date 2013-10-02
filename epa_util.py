#! /usr/bin/env python
import sys
import os
import json
import operator
import time
from ete2 import Tree, TreeStyle, TextFace, SeqGroup
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
        self.name = str(time.time())
    
    
    def run(self, reftree, alignment, num_thread = "2"):
        #./raxmlHPC-PTHREADS-SSE3 -f v -m GTRGAMMA -s ../example/t1.fa -t ../example/t1.test.tre -p 1234 -n jjj -T 2
        if os.path.exists(reftree):
            call([self.raxmlpath,"-f", "v", "-m","GTRGAMMA","-s",alignment,"-t", reftree, "-n",self.name,"-p", "1234", "-T", num_thread, "-w", self.tmppath] ) #, stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
        else:
            tmptree = self.tmppath+"/" + self.name + ".tre"
            with open(tmptree, "w") as fout:
                fout.write(reftree)
            call([self.raxmlpath,"-f", "v", "-m","GTRGAMMA","-s",alignment,"-t", tmptree, "-n",self.name,"-p", "1234", "-T", num_thread, "-w", self.tmppath] ) #, stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
            
        jp = jsonparser(self.tmppath + "/" + "RAxML_portableTree." + self.name + ".jplace")
        return jp.get_placement()
    
    
    def clean(self):
        os.remove(self.tmppath + "/" + "RAxML_classification." + self.name)
        os.remove(self.tmppath + "/" + "RAxML_classificationLikelihoodWeights." + self.name)
        os.remove(self.tmppath + "/" + "RAxML_entropy." + self.name)
        os.remove(self.tmppath + "/" + "RAxML_info." + self.name)
        os.remove(self.tmppath + "/" + "RAxML_labelledTree." + self.name)
        os.remove(self.tmppath + "/" + "RAxML_originalLabelledTree." + self.name)
        #os.remove(self.tmppath + "/" + "RAxML_portableTree." + self.name + ".jplace")


if __name__ == "__main__":
    print("This is main")
    EPA = epa()
    place = EPA.run(reftree = "/home/zhangje/GIT/EPA-classifier/example/t1.test.tre", alignment = "/home/zhangje/GIT/EPA-classifier/example/t1.fa")
    EPA.clean()
    print(place)
