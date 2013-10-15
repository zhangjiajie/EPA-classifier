#! /usr/bin/env python
import sys
import os
import json
import operator
import time
from ete2 import Tree, TreeStyle, TextFace, SeqGroup
from subprocess import call
from epac.json_util import jsonparser
from epac.taxonomy_parser import TreeBuilder
from epac.epa_util import raxml, epa
from epac.taxonomy_training import phylogeny_annotator, trainning


def auto_train(ref_taxonomy, ref_sequences, outfile):
    basepath = os.path.dirname(os.path.abspath(__file__))
    tmppath = basepath + "/epac/tmp"
    if not os.path.exists(tmppath):
        print("The tmp folder for keeping intermediate files does not exit")
        print("Please create this folder:" + tmppath)
        sys.exit()
    name = str(time.time())
    
    mftree_name = tmppath + "/" + name + ".tre"
    tb = TreeBuilder(flat_taxonomy_file = ref_taxonomy, fout = mftree_name)
    tb.build()
    
    rx = raxml()
    bf_tree_name = rx.resolve_mftree(mftree_name, ref_sequences)
    os.remove(mftree_name)
    rx.clean()
    
    #1. root this tree 
    pa = phylogeny_annotator(sphylogeny = bf_tree_name, s_seq_db = ref_taxonomy)
    pa.annotate()
    rooted_raxml_tree = pa.root.write(format=5)
    os.remove(bf_tree_name)
    
    #2. dummy epa to label branches
    EPA = epa()
    dummy_jplace = EPA.dummy(reftree = rooted_raxml_tree, alignment = ref_sequences)
    
    #3. annotate the epa labeld tree 
    t = trainning(temp_epa_json = dummy_jplace, ref_taxonomy = ref_taxonomy, ref_sequences = ref_sequences)
    t.trainning2json(fout = outfile)
    
    
    
if __name__ == "__main__":
    print("This is main")
    auto_train(ref_taxonomy = "/home/zhangje/GIT/EPA-classifier/example/training_tax.txt", ref_sequences = "/home/zhangje/GIT/EPA-classifier/example/training_seq.fa", outfile = "/home/zhangje/GIT/EPA-classifier/example/auto1.json")
