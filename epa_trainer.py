#! /usr/bin/env python
import sys
import os
import json
import operator
import time
from ete2 import Tree, TreeStyle, TextFace, SeqGroup
from subprocess import call
from json_util import jsonparser
from taxonomy_parser import TreeBuilder
from epa_util import raxml


def auto_train(ref_taxonomy, ref_sequences):
    basepath = os.path.dirname(os.path.abspath(__file__))
    tmppath = self.basepath + "/tmp"
    if not os.path.exists(self.tmppath):
        print("The tmp folder for keeping intermediate files does not exit")
        print("Please create this folder:" + self.tmppath)
        sys.exit()
    name = str(time.time())
    
    mftree_name = tmppath + "/" + name + ".tre"
    tb = TreeBuilder(flat_taxonomy_file = ref_taxonomy, fout = mftree_name)
    tb.build()
    
    rx = raxml()
    bf_tree_name = rx.resolve_mftree(mftree_name, ref_sequences)
    #os.remove(mftree_name)
    
    #1. root this tree 
    
    #2. dummy epa to label branches
    
    #3. annotate the epa labeld tree 
    
    #4. output json file  
    
    
if __name__ == "__main__":
    print("This is main")
