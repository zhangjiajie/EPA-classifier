#! /usr/bin/env python
import sys
import os
import json
import operator
from ete2 import Tree, TreeStyle, TextFace, SeqGroup
from subprocess import call


class jsonparser:
    def __init__(self, jsonfin):
        self.jdata = json.load(open(jsonfin))
    
    def get_raxml_readable_tree(self, fout = None):
        t = Tree(self.jdata["tree"], format=1)
        #t.unroot()
        if fout!=None:
            t.write(outfile=fout, format=5)
        else:
            return t.write(format=5)
    
    def get_reftree(self):
        t = Tree(self.jdata["tree"], format=1)
        return t
    
    def get_bid_tanomomy_map(self):
        return self.jdata["taxonomy"]
    
    def get_alignment(self, fout):
        soutput = ""
        entries = self.jdata["sequences"]
        for entr in entries:
            soutput = soutput + ">" + entr[0] + "\n" + entr[1] + "\n"
        with open(fout, "w") as fo:
            fo.write(soutput)
        return fout
    
    def get_alignment_list(self):
        return self.jdata["sequences"]
    
    def get_sequences_names(self):
        nameset = set()
        entries = self.jdata["sequences"]
        for entr in entries:
            nameset.add(entr[0])
        return nameset
    
    def get_alignment_length(self):
        entries = self.jdata["sequences"]
        return len(entries[0][1])
    
    def get_placement(self):
        return self.jdata["placements"]
        
    def get_std_newick_tree(self):
        tree = self.jdata["tree"]
        tree = tree.replace("{", "[&&NHX:B=")
        tree = tree.replace("}", "]")
        return tree
    
    def get_hmm_profile(self, fout):
        lines = self.jdata["hmmprofile"]
        with open(fout, "w") as fo:
            for line in lines:
                fo.write(line)

if __name__ == "__main__":
    jp = jsonparser("tt.json")
    jp.get_raxml_readable_tree("jsontt.tre")
    jp.get_alignment("jsontt.fa")
    print(jp.get_bid_tanomomy_map())
