#! /usr/bin/env python
import sys
import os
import json
import operator
from ete2 import Tree, TreeStyle, TextFace, SeqGroup
from subprocess import call
from epa_util import epa
from json_util import jsonparser
from pprint import pprint

class magic:
    def __init__(self, refjson, query):
        self.refjson = jsonparser(refjson)
        self.bid_taxonomy_map = self.refjson.get_bid_tanomomy_map()
        self.refree = self.refjson.get_reftree()
        self.query = query


    def alignment(self):
        pass


    def checkinput(self):
        return True


    def print_ranks(self, rks):
        ss = ""
        for rk in rks:
            ss = ss + rk + ";"
        return ss[:-1]


    def print_confi(self, confs):
        ss = ""
        for rk in confs:
            if rk == confs[0]:
                rk = 1.0
            ss = ss + repr(rk) + ";"
        return ss[:-1]


    def classify(self, fout = None):
        self.checkinput()
        EPA = epa()
        placements = EPA.run(reftree = self.refjson.get_raxml_readable_tree(), alignment = self.query)
        EPA.clean()
        if fout!=None:
            fo = open(fout, "w")
        
        for place in placements:
            taxa_name = place["n"][0]
            edges = place["p"]
            ranks, lws = self.assign_taxonomy(edges)
            output = taxa_name+ "\t" + self.print_ranks(ranks) + "\t" + self.print_confi(lws) + "\n"
            print(output) 
            if fout!=None:
                fo.write(output)
        if fout!=None:
            fo.close()


    def assign_taxonomy(self, edges):
        #Calculate the sum of likelihood weight for each rank
        taxonmy_sumlw_map = {}
        for edge in edges:
            edge_nr = str(edge[0])
            lw = edge[2]
            taxonomy = self.bid_taxonomy_map[edge_nr]
            for rank in taxonomy:
                if rank == "-":
                    taxonmy_sumlw_map[rank] = -1
                elif rank in taxonmy_sumlw_map:
                    oldlw = taxonmy_sumlw_map[rank]
                    taxonmy_sumlw_map[rank] = oldlw + lw
                else:
                    taxonmy_sumlw_map[rank] = lw
        
        #Assignment using the max likelihood placement
        ml_edge = edges[0]
        edge_nr = str(ml_edge[0])
        maxlw = ml_edge[2]
        ml_ranks = self.bid_taxonomy_map[edge_nr]
        ml_ranks_copy = []
        for rk in ml_ranks:
            ml_ranks_copy.append(rk)
        lws = []
        cnt = 0
        for rank in ml_ranks:
            lw = taxonmy_sumlw_map[rank]
            if lw > 1.0:
                lw = 1.0
            lws.append(lw)
            if rank == "-" and cnt > 0 :                
                for edge in edges[1:]:
                    edge_nr = str(edge[0])
                    taxonomy = self.bid_taxonomy_map[edge_nr]
                    newrank = taxonomy[cnt]
                    newlw = taxonmy_sumlw_map[newrank]
                    higherrank_old = ml_ranks[cnt -1]
                    higherrank_new = taxonomy[cnt -1]
                    if higherrank_old == higherrank_new and newrank!="-":
                        ml_ranks_copy[cnt] = newrank
                        lws[cnt] = newlw
            cnt = cnt + 1
            
        return ml_ranks_copy, lws


if __name__ == "__main__":
    print("This is main")
    m = magic("/home/zhangje/GIT/EPA-classifier/example/tt.json", "/home/zhangje/GIT/EPA-classifier/example/t1.fa")
    m.classify(fout = "/home/zhangje/GIT/EPA-classifier/example/taxout.txt")
    
