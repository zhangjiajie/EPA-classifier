#!/usr/bin/env python

import os
import sys
import glob
import shutil
import datetime
import time
import ConfigParser

class DefaultedConfigParser(ConfigParser.SafeConfigParser):
    def get_param(self, section, option, ctype=str, default=None):
        if default is None:
            ret = self.get(section, option)
        else:
            confdict = self.__dict__.get('_sections')
            sec = confdict.get(section)
            if sec:
                if ctype == bool:
                    ret = self.getboolean(section, option)
                else:
                    ret = sec.get(option, default)
            else:
                ret = default
        return ctype(ret)

class EpacConfig:
    # this prefix will be added to every sequence name in reference to prevent 
    # name clashes with query sequences, which are coded with numbers
    REF_SEQ_PREFIX = "r_";
    QUERY_SEQ_PREFIX = "q_";

    def __init__(self):
        self.set_defaults()
        
    def __init__(self, args): 
        self.verbose = args.verbose
        self.debug = args.debug
        self.refjson_fname = args.ref_fname        
        self.basepath = os.path.dirname(os.path.abspath(__file__))
        self.epac_home = os.path.abspath(os.path.join(self.basepath, os.pardir)) + "/"
        self.reftree_home = os.path.abspath("reftree/") + "/"
        self.temp_dir = self.basepath + "/tmp/"
        self.raxml_outdir = self.temp_dir #"raxml_output/"
        self.raxml_outdir_abs = os.path.abspath(self.raxml_outdir)
        self.results_home = os.path.abspath("results/") + "/"
        self.set_defaults()
        if args.config_fname:
            self.read_from_file(args.config_fname)
        # command line setting has preference over config file and default
        if args.num_threads:
            self.num_threads = args.num_threads        
        self.check_raxml()    
        if args.output_name:
            self.name = args.output_name
        else:
            self.name = "%d" % (time.time()*1000) #datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        results_name = self.name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        self.results_dir = self.results_home + results_name + "/"

    def set_defaults(self):
        self.muscle_home = self.epac_home + "/epac/bin" + "/"
        self.hmmer_home = self.epac_home + "/epac/bin" + "/"
        self.raxml_home = self.epac_home + "/epac/bin" + "/"
        self.raxml_exec = "raxmlHPC-PTHREADS-SSE3"
        self.raxml_model = "GTRCAT"
        self.raxml_remote_host = ""
        self.raxml_remote_call = False        
        self.run_on_cluster = False
        self.cluster_epac_home = self.epac_home
        self.cluster_qsub_script = ""
        self.epa_load_optmod = True
        self.epa_use_heuristic = "auto"
        self.epa_heur_rate = 0.01
        self.min_confidence = 0.2
        self.num_threads = 2
        
    def check_raxml(self):
        self.raxml_exec_full = self.raxml_home + self.raxml_exec
        if self.raxml_remote_host in ["", "localhost"]:
            self.raxml_remote_call = False
            if not os.path.isdir(self.raxml_home):
                print "RAxML home directory not found: %s" % self.raxml_home
                sys.exit()
            if not os.path.isfile(self.raxml_exec_full):
                print "RAxML executable not found: %s" % self.raxml_exec_full
                sys.exit()
        else:
            self.raxml_remote_call = True
        self.raxml_cmd = [self.raxml_exec_full, "-p", "12345", "-T", str(self.num_threads), "-w", self.raxml_outdir_abs]
        
    def read_from_file(self, config_fname):
        if not os.path.exists(config_fname):
            print "Config file not found: " + config_fname
            sys.exit()

        parser = DefaultedConfigParser() #ConfigParser.SafeConfigParser()
        parser.read(config_fname)
        
        self.raxml_home = parser.get_param("raxml", "raxml_home", str, self.raxml_home) + "/"
        self.raxml_exec = parser.get_param("raxml", "raxml_exec", str, self.raxml_exec)
        self.raxml_remote_host = parser.get_param("raxml", "raxml_remote_host", str, self.raxml_remote_host)

        self.raxml_model = parser.get_param("raxml", "raxml_model", str, self.raxml_model)
        self.num_threads = parser.get_param("raxml", "raxml_threads", int, self.num_threads)

        self.epa_use_heuristic = parser.get_param("raxml", "epa_use_heuristic", str, self.epa_use_heuristic)
        self.epa_heur_rate = parser.get_param("raxml", "epa_heur_rate", float, self.epa_heur_rate)
        self.epa_load_optmod = parser.get_param("raxml", "epa_load_optmod", bool, self.epa_load_optmod)

        self.hmmer_home = parser.get_param("hmmer", "hmmer_home", str, self.hmmer_home)
        self.muscle_home = parser.get_param("muscle", "muscle_home", str, self.muscle_home)
        
        self.run_on_cluster = parser.get_param("cluster", "run_on_cluster", bool, self.run_on_cluster)
        self.cluster_epac_home = parser.get_param("cluster", "cluster_epac_home", str, self.cluster_epac_home) + "/"
        self.cluster_qsub_script = parser.get_param("cluster", "cluster_qsub_script", str, self.cluster_qsub_script)

        self.min_confidence = parser.get_param("assignment", "min_confidence", float, self.min_confidence)

        return parser

    def subst_name(self, in_str):
        """Replace %NAME% macros with an actual EPAC run name. Used to 
        generate unique run-specific identifiers (filenames, RAxML job names etc)"""
        return in_str.replace("%NAME%", self.name)
    
    def tmp_fname(self, fname):
        return self.temp_dir + self.subst_name(fname)
    
class EpacTrainerConfig(EpacConfig):
    
    def __init__(self, args):
        self.taxonomy_fname = args.taxonomy_fname
        self.align_fname = args.align_fname
        self.no_hmmer = args.no_hmmer
        EpacConfig.__init__(self, args)
        
    def set_defaults(self):
        EpacConfig.set_defaults(self)
        # default settings below imply no taxonomy filtering, 
        # i.e. all sequences from taxonomy file will be included into reference tree
        self.reftree_min_rank = 0
        self.reftree_max_seqs_per_leaf = 1e6
        self.reftree_clades_to_include=[]
        self.reftree_clades_to_ignore=[]

    def read_from_file(self, config_fname):
        parser = EpacConfig.read_from_file(self, config_fname)
        
        self.reftree_min_rank = parser.get_param("reftree", "min_rank", int, self.reftree_min_rank)
        self.reftree_max_seqs_per_leaf = parser.get_param("reftree", "max_seqs_per_leaf", int, self.reftree_max_seqs_per_leaf)
        clades_str = parser.get_param("reftree", "clades_to_include", str, "")
        self.reftree_clades_to_include = self.parse_clades(clades_str)
        clades_str = parser.get_param("reftree", "clades_to_ignore", str, "")
        self.reftree_clades_to_ignore = self.parse_clades(clades_str)
        
    def parse_clades(self, clades_str):
        clade_list = []
        try:        
            if clades_str:
                clades = clades_str.split(",")
                for clade in clades:
                    toks = clade.split("|")
                    clade_list += [(int(toks[0]), toks[1])]
        except:
            print "Invalid format in config parameter: clades_to_include"
            sys.exit()

        return clade_list
        
