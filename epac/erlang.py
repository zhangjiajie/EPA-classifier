#! /usr/bin/env python
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
        for name in self.taxonomy.keys():
            ranks = self.taxonomy[name]
            sp = ranks[-1]
            if sp == "-":
                keepseqs.append(name)
            else:
                if not sp in species:
                    keepseqs.append(name)
                    species.add(sp)
        root = Tree(self.tree)
        root.prune(keepseqs, preserve_branch_length=True)
        sumbr = 0.0
        cnt = 0.0 
        for node in root.traverse(stratagy = "preorder"):
            sumbr = sumbr + node.dist
            cnt = cnt + 1.0
        return float(cnt) / float(sumbr)
       
    def get_nodesheight(self):
        nh_map = {}
        for node in root.traverse(stratagy = "preorder"):
            if hasattr(node, "B"):
                height = node.get_farthest_leaf(topology_only=True)
                nh_map[node.B] = height[1] + 1
        
        return nh_map

if __name__ == "__main__":
    print("This is erlang.py main")
    el = erlang()
    print el.one_tail_test(rate = 10, k = 2, x = 1)
