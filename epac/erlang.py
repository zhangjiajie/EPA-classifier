#! /usr/bin/env python
import sys
import os
import json
import operator
import math
from ete2 import Tree

class erlang:
    def __init__(self):
        pass
    
    def one_tail_test(self, rate, k, x):
        """rate: estimated branching rate from reference tree
           k: node height
           x: placement branch length"""
        p = 0.0
        for n in range(k):
            p = p + (1.0/float(math.factorial(n))) * math.exp((-rate)*x) * math.pow(rate*x, n)
        return p

class tree_param:
    def __init__(self, tree, origin_taxonomy):
        """tree: rooted and branch labled tree in newick format
            origin_taxonomy: a dictionary of leaf name and taxonomy ranks"""
        self.tree = tree 
        self.taxonomy = origin_taxonomy 
    
    def get_speciation_rate(self):
        #pruning the input tree such that each speices only appear once
        species = set()
        keepseqs = []
        for name is self.taxonomy.keys():
            ranks = self.taxonomy[name]
            sp = ranks[-1]
            if sp == "-":
                keepseqs.append(name)
            else:
                if not sp in species:
                    keepseqs.append(name)
                    species.add(sp)
        #TODO
        




if __name__ == "__main__":
    print("This is erlang.py main")
    el = erlang()
    print el.one_tail_test(rate = 10, k = 2, x = 1)
