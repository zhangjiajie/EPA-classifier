#! /usr/bin/env python
import sys
import os
import json
import operator
from epac.ete2 import Tree, SeqGroup
from subprocess import call


class RefJsonChecker:
    def __init__(self, jsonfin= None, jdata = None):
        if jsonfin!=None:
            self.jdata = json.load(open(jsonfin))
        else:
            self.jdata = jdata
    
    def valid(self, ver = "1.0"):
        #tree
        if "tree" in self.jdata:
            tree = self.jdata["tree"]
            if not isinstance(tree, unicode):
                print("Tree is")
                print(type(tree).__name__)
                return False
        else:
            return False
        
        #raxmltree
        if "raxmltree" in self.jdata:
            tree = self.jdata["raxmltree"]
            if not isinstance(tree, unicode):
                return False
        else:
            return False
            
        #rate:
        if "rate" in self.jdata:
            rate = self.jdata["rate"]
            if not isinstance(rate, float):
                return False
        else:
            return False
            
        #node_height
        if "node_height" in self.jdata:
            node_height = self.jdata["node_height"]
            if not isinstance(node_height, dict):
                return False
        else:
            return False
            
        #taxonomy
        if "taxonomy" in self.jdata:
            taxonomy = self.jdata["taxonomy"]
            if not isinstance(taxonomy, dict):
                return False
        else:
            return False

        #origin_taxonomy
        if "origin_taxonomy" in self.jdata:
            origin_taxonomy = self.jdata["origin_taxonomy"]
            if not isinstance(origin_taxonomy, dict):
                return False
        else:
            return False

        #sequences
        if "sequences" in self.jdata:
            sequences = self.jdata["sequences"]
            if not isinstance(sequences, list):
                return False
        else:
            return False
            
        #hmm_profile
        if "hmm_profile" in self.jdata:
            hmm_profile = self.jdata["hmm_profile"]
            if not isinstance(hmm_profile, list):
                return False
        else:
            return False
        
        return True

class EpaJsonParser:
    """This class parses the RAxML-EPA json output file"""
    def __init__(self, jsonfin):
        self.jdata = json.load(open(jsonfin))
    
    def get_placement(self):
        return self.jdata["placements"]
        
    def get_tree(self):
        return self.jdata["tree"]
        
    def get_std_newick_tree(self):
        tree = self.jdata["tree"]
        tree = tree.replace("{", "[&&NHX:B=")
        tree = tree.replace("}", "]")
        return tree

class RefJsonParser:
    """This class parses the EPA Classifier reference json file"""
    def __init__(self, jsonfin):
        self.jdata = json.load(open(jsonfin))
        
    def validate(self):
        jc = RefJsonChecker(jdata = self.jdata)
        if not jc.valid():
            print("Invalid reference database format")
            sys.exit()
    
    def get_rate(self):
        return self.jdata["rate"]
    
    def get_node_height(self):
        return self.jdata["node_height"]
    
    def get_raxml_readable_tree(self, fout_name = None):
        tree_str = self.jdata["raxmltree"]
        #t.unroot()
        if fout_name != None:
            with open(fout_name, "w") as fout:
                fout.write(tree_str)
        else:
            return tree_str
    
    def get_reftree(self):
        t = Tree(self.jdata["tree"], format=1)
        return t
    
    def get_bid_tanomomy_map(self):
        return self.jdata["taxonomy"]

    def get_origin_taxonomy(self):
        return self.jdata["origin_taxonomy"]
    
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
    
    def get_hmm_profile(self, fout):
        lines = self.jdata["hmm_profile"]
        with open(fout, "w") as fo:
            for line in lines:
                fo.write(line)
                
class RefJsonBuilder:
    """This class builds the EPA Classifier reference json file"""
    def __init__(self):
        self.jdata = {}
        self.jdata["version"] = "1.0"
        self.jdata["author"] = "Jiajie Zhang"
        
    def set_taxonomy(self, bid_ranks_map):
        self.jdata["taxonomy"] = bid_ranks_map

    def set_origin_taxonomy(self, orig_tax_map):
        self.jdata["origin_taxonomy"] = orig_tax_map

    def set_tree(self, tr):
        self.jdata["tree"] = tr
        self.jdata["raxmltree"] = Tree(tr, format=1).write(format=5)
        
    def set_sequences(self, seqs):    
        self.jdata["sequences"] = seqs
        
    def set_hmm_profile(self, fprofile):    
        with open(fprofile) as fp:
            lines = fp.readlines()
        self.jdata["hmm_profile"] = lines
       
    def set_rate(self, rate):    
        self.jdata["rate"] = rate
        
    def set_nodes_height(self, height):    
        self.jdata["node_height"] = height

    def dump(self, out_fname):
        with open(out_fname, "w") as fo:
            json.dump(self.jdata, fo, indent=4, sort_keys=True)                


if __name__ == "__main__":
    if len(sys.argv) < 2: 
        print("usage: ./json_util.py jsonfile")
        sys.exit()
    jc = json_checker(jsonfin = sys.argv[1])
    if jc.valid():
        print("The json file is OK for EPA-classifer")
    else:
        print("!!!Invalid json file!!!")
    
