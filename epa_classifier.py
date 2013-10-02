#! /usr/bin/env python
import sys
import os
import json
import operator
from ete2 import Tree, TreeStyle, TextFace, SeqGroup
from subprocess import call
from epa_util import epa
from json_util import jsonparser

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
    
    def classify(self, fout = None):
        self.checkinput()
        EPA = epa()
        placements = EPA.run(reftree = self.refjson.get_raxml_readable_tree(), alignment = self.query)
        EPA.clean()
        














if __name__ == "__main__":
    print("This is main")
