#! /usr/bin/env python
import sys
import os
import json
import operator
import time
import subprocess
from epac.ete2 import SeqGroup, Tree
from subprocess import call
#from json_util import jsonparser
from PTP_light.PTP import EPA_interface
from raxml_util import FileUtils


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
        self.alignment = ""
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
    
    
    def raxml(self, alignment, num_thread = "2"):
        call([self.raxmlpath, "-m","GTRGAMMA","-s",alignment, "-n",self.name,"-p", "1234", "-T", num_thread, "-w", self.tmppath], stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
        self.alignment = alignment
        return self.tmppath + "/" + "RAxML_bestTree." + self.name
    
    
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
        FileUtils.remove_if_exists(self.tmppath + "/" + "RAxML_info." + self.name)
        #print(self.tmppath + "/" + "RAxML_info." + self.name)
        FileUtils.remove_if_exists(self.tmppath + "/" + "RAxML_log." + self.name)
        FileUtils.remove_if_exists(self.tmppath + "/" + "RAxML_result." + self.name)
        FileUtils.remove_if_exists(self.tmppath + "/" + "RAxML_bestTree." + self.name)
        FileUtils.remove_if_exists(self.tmppath + "/" + "RAxML_parsimonyTree." + self.name)
        FileUtils.remove_if_exists(self.alignment)
        FileUtils.remove_if_exists(self.alignment+".reduced")
        


def epa_2_ptp(epa_jp, ref_jp, full_alignment, min_lw = 0.5, debug = False):
    placements = epa_jp.get_placement()
    reftree = Tree(epa_jp.get_std_newick_tree())
    allnodes = reftree.get_descendants()
    species_list = []
    
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

    groups = placemap.items()
    cnt_leaf = 0
    cnt_inode = 0
    
    """check each placement edge""" 
    for i,item in enumerate(groups):
        place_branch_name = item[0]
        seqset = item[1]
        if len(seqset) < 4:
            species_list.append(seqset)
        else:
            branch_alignment = SeqGroup()
            for taxa in seqset:
                branch_alignment.set_seq(taxa, full_alignment.get_seq(taxa))
            species = build_tree_run_ptp(branch_alignment, ref_jp.get_rate())
            species_list.extend(species)
    return species_list



def build_tree_run_ptp(alignment, sp_rate):
    """return a list of sets of taxa names, and remove the neighoring taxa name"""
    rml = raxml()
    name = str(time.time()) + ".afa"
    outname = rml.tmppath + "/" + name
    alignment.write(outfile = outname)
    all_taxa = []
    for seq in alignment:
        all_taxa.append(seq[0])
        
    outtree = rml.raxml(alignment = outname)
    if os.path.exists(outtree):
        slist = EPA_interface(tree = outtree, sp_rate = sp_rate, reroot = True, method = "H0", max_iters = 20000, min_brl = 0.0001, pvalue = 0.001)
        rml.clean()
        return slist
    else:
        rml.clean()
        return [all_taxa]
    
if __name__ == "__main__":
    if len(sys.argv) < 3: 
        print("usage: ./epa_util.py <multifurcating.tre> <alignment> <num_thread>")
        sys.exit()
    rx = raxml()
    tree = rx.resolve_mftree(sys.argv[1], sys.argv[2], sys.argv[3])
    rx.clean()
    print("The resolved tree has been written to: ")
    print(tree)
