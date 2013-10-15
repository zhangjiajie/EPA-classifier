#!/usr/bin/env python
import sys
from ete2 import Tree, TreeNode


class Taxa:
    def __init__(self, ranks, node_list):
        """a list of ranks from high to low, and node_list of the entire tree"""
        self.ranks = ranks
        for i in range(len(self.ranks)):
            self.ranks[i] = self.ranks[i].strip()
        self.node_list = node_list
        
    def create_nodes(self, father):
        root = None
        lastNode = None
        i = 0
        has_found = False
        for rank in self.ranks:
            if rank == father:
                i = i + 1
                has_found = True
                break
            else:
                i = i + 1
        if not has_found:
            i = 0
        for rank in self.ranks[i:]:
            newnode = Tree(name = rank)
            self.node_list.append(newnode)
            if lastNode!=None:
                lastNode.add_child(child=newnode)
                lastNode = newnode
            else:
                root = newnode
                lastNode = newnode
        return root
        
    def __find_rank_node_list(self, rank):
        for node in self.node_list:
            if rank.strip() == node.name:
                return True, node
        return False, None
    
    def find_my_farther(self):
        for rank in reversed(self.ranks):
            hasfound, node = self.__find_rank_node_list(rank)
            if hasfound:
                return node
        return None


class TreeBuilder:
    """Build a multifurcating taxonomy tree from greengene flat taxonomy file"""
    def __init__(self, flat_taxonomy_file, fout):
        with open(flat_taxonomy_file) as f:
            self.lines = f.readlines()
        self.fout = fout
    
    
    def build(self):
        root = TreeNode(name = "root") 
        node_list = []
        node_list.append(root)
        for line in self.lines:
            toks = line.strip().split("\t")
            seqid = toks[0]
            ranks = toks[1].split(";")
            ranks.append(seqid)
            taxa = Taxa(ranks, node_list)
            farther_node = taxa.find_my_farther()
            if farther_node != None:
                new_subtree = taxa.create_nodes(farther_node.name)
                farther_node.add_child(child = new_subtree)
            else:
                new_subtree = taxa.create_nodes("root")
                root.add_child(child = new_subtree)
        root = root.children[0]
        root.up = None
        while len(root.children) == 1:
            root = root.children[0]
            root.up = None
        leaves = root.get_leaf_names()
        root.prune(leaves)
        root.unroot()
        root.write(outfile=self.fout, format=6)
        return root
                
            
if __name__ == "__main__":
    if len(sys.argv) != 3: 
        print("usage: ./taxonomy_parser.py <greengene_taxonomy> <outputfile>")
        sys.exit()
    
    tb = TreeBuilder(flat_taxonomy_file = sys.argv[1], fout = sys.argv[2])
    tree = tb.build()
    tree.show()
    
