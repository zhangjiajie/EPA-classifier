#!/usr/bin/env python

import os
import sys
import glob
import shutil
import datetime
import time
import ConfigParser

class EpacConfig:
    F_ALIGN_REF=1
    F_ALIGN_LBLQ=2
    F_ALIGN_QUERY=3

    F_TREE_MULTIF=11
    F_TREE_MULTIF_TAX=12
    F_TREE_OUTGR=13
    F_TREE_BIF_UNROOTED=14
    F_TREE_BIF_UNROOTED_LBL=15
    F_TREE_BIF_ROOTED_LBL=16
    F_TREE_BIF_ROOTED=17
    F_TREE_BIF_ROOTED_TAX=18
    F_TREE_EPA_RESULT=19
    F_TREE_EPA_RESULT_TAX=20

    F_SEQLIST_REF=31
    F_SEQLIST_LBLQ=32

    F_TAXONOMY_REF=41

    F_MAP_BRANCH_RANK=51
    F_OPT_MODEL=52

    F_QSUB_SCRIPT = 61,

    F_RESULT_TAX_ASSIGN=71
    F_RESULT_STATS=72
    F_RESULT_MIS=73
    F_RESULT_BRANCH_ASSIGN=74
    F_RESULT_LH_WEIGHTS=75

    def __init__(self):
        self.set_defaults()
        
    def __init__(self, args): 
        self.verbose = args.verbose
        self.debug = args.debug
        self.refjson_fname = args.ref_fname        
        self.num_threads = args.num_threads        
        self.epac_home = os.path.abspath("") + "/"
        self.basepath = os.path.dirname(os.path.abspath(__file__))
        self.reftree_home = os.path.abspath("reftree/") + "/"
        self.temp_dir = self.basepath + "/tmp/"
        self.raxml_outdir = self.temp_dir #"raxml_output/"
        self.raxml_outdir_abs = os.path.abspath(self.raxml_outdir)
        self.results_home = os.path.abspath("results/") + "/"
        self.set_defaults()
        if args.config_fname:
            self.read_from_file(args.config_fname)
        self.name = datetime.datetime.now().strftime("%Y%m%d_%H%M%S") #str(time.time())
        results_name = self.name + "_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        self.results_dir = self.results_home + results_name + "/"

    def set_defaults(self):
        self.muscle_home = self.epac_home + "/epac/bin" + "/"
        self.hmmer_home = self.epac_home + "/epac/bin" + "/"
        self.raxml_home = self.epac_home + "/epac/bin" + "/"
        self.raxml_exec = "raxmlHPC-PTHREADS-SSE3"
        self.raxml_remote_host = ""
    
    def read_from_file(self, config_fname):
        if not os.path.exists(config_fname):
            print "Config file not found: " + config_fname
            sys.exit()

        parser = ConfigParser.SafeConfigParser()
        parser.read(config_fname)
        
        try:
            self.raxml_home = parser.get("raxml", "raxml_home") + "/"
            self.raxml_exec = parser.get("raxml", "raxml_exec")
            self.raxml_remote_host = parser.get("raxml", "raxml_remote_host")
        except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
            pass

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

        self.raxml_model = parser.get("raxml", "raxml_model")
        if not self.num_threads:
            self.num_threads = parser.getint("raxml", "raxml_threads")
        self.raxml_cmd = [self.raxml_exec_full, "-p", "12345", "-T", str(self.num_threads), "-w", self.raxml_outdir_abs]

        self.epa_use_heuristic = parser.getboolean("raxml", "epa_use_heuristic")
        self.epa_heur_rate = parser.getfloat("raxml", "epa_heur_rate")
        self.epa_load_optmod = parser.getboolean("raxml", "epa_load_optmod")

        try:
            self.hmmer_home = parser.get("hmmer", "hmmer_home")
            self.muscle_home = parser.get("muscle", "muscle_home")
        except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
            pass
        
        try:        
            self.run_on_cluster = parser.getboolean("cluster", "run_on_cluster")
            self.cluster_epac_home = parser.get("cluster", "cluster_epac_home") + "/"
            self.cluster_qsub_script = parser.get("cluster", "cluster_qsub_script")
        except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
            self.run_on_cluster = False

        self.min_confidence = parser.getfloat("assignment", "min_confidence")
        return parser

    def tmp_fname(self, fname):
        return self.temp_dir + fname.replace("%NAME%", self.name)
    
    def get_fname(self, file_type):
        return {
            EpataxConfig.F_ALIGN_REF: self.reftree_dir + self.reftree_name + ".fasta",
            EpataxConfig.F_ALIGN_LBLQ: self.temp_dir + self.reftree_name + "_lblq.fasta",
            EpataxConfig.F_ALIGN_QUERY: self.temp_dir + self.reftree_name + "_query.fasta",

            EpataxConfig.F_TREE_MULTIF: self.temp_dir + self.reftree_name + "_mfu.tre",
            EpataxConfig.F_TREE_MULTIF_TAX: self.temp_dir + self.reftree_name + "_mfu_tax.tre",
            EpataxConfig.F_TREE_OUTGR: self.temp_dir + self.reftree_name + "_outgr.tre",
            EpataxConfig.F_TREE_BIF_UNROOTED: self.reftree_dir + self.reftree_name + ".tre",
            EpataxConfig.F_TREE_BIF_UNROOTED_LBL: self.temp_dir + self.reftree_name + "_bfu_lbl.tre",
            EpataxConfig.F_TREE_BIF_ROOTED_LBL: self.temp_dir + self.reftree_name + "_bfr_lbl.tre",
            EpataxConfig.F_TREE_BIF_ROOTED: self.temp_dir + self.reftree_name + "_bfr.tre",
            EpataxConfig.F_TREE_BIF_ROOTED_TAX: self.temp_dir + self.reftree_name + "_bfr_tax.tre",
            EpataxConfig.F_TREE_EPA_RESULT: self.results_dir + self.reftree_name + "_epa.tre",
            EpataxConfig.F_TREE_EPA_RESULT_TAX: self.results_dir + self.reftree_name + "_epa_tax.tre",

            EpataxConfig.F_SEQLIST_REF: self.reftree_dir + self.reftree_name + "_seq.txt",
            EpataxConfig.F_SEQLIST_LBLQ: self.data_dir + "greengenes_lblq_seq.txt",

            EpataxConfig.F_TAXONOMY_REF: self.reftree_dir + self.reftree_name + "_tax.txt",

            EpataxConfig.F_MAP_BRANCH_RANK: self.reftree_dir + self.reftree_name + ".map",
            EpataxConfig.F_OPT_MODEL: self.reftree_dir + self.reftree_name + ".opt",

            EpataxConfig.F_QSUB_SCRIPT: self.temp_dir + self.reftree_name + "_sub.sh",

            EpataxConfig.F_RESULT_TAX_ASSIGN: self.results_dir + self.reftree_name + ".assigned.txt",
            EpataxConfig.F_RESULT_STATS: self.results_dir + self.reftree_name + ".stats.txt",
            EpataxConfig.F_RESULT_MIS: self.results_dir + self.reftree_name + ".mis.txt"
            }[file_type]
            
class EpacTrainerConfig(EpacConfig):
    
    def __init__(self, args):
        self.taxonomy_fname = args.taxonomy_fname
        self.align_fname = args.align_fname
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
        
        try:
            self.reftree_min_rank = parser.getint("reftree", "min_rank")
            self.reftree_max_seqs_per_leaf = parser.getint("reftree", "max_seqs_per_leaf")
            clades_str = parser.get("reftree", "clades_to_include")
            self.reftree_clades_to_include = self.parse_clades(clades_str)
            clades_str = parser.get("reftree", "clades_to_ignore")
            self.reftree_clades_to_ignore = self.parse_clades(clades_str)
        except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
            pass            
        
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
        
