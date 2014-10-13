#! /usr/bin/env python
try:
    import sys
    import os
    import time
    import glob
    from operator import itemgetter
    from epac.ete2 import Tree, SeqGroup
    from epac.argparse import ArgumentParser
    from epac.config import EpacConfig
    from epac.raxml_util import RaxmlWrapper, FileUtils
    from epac.json_util import RefJsonParser, RefJsonChecker, EpaJsonParser
    from epac.taxonomy_util import Taxonomy,GGTaxonomyFile
    from epac.classify_util import TaxClassifyHelper
except ImportError, e:
    print("Some packages are missing, please re-downloand EPA-classifier")
    print e
    sys.exit()


class LeaveOneTest:
    def __init__(self, config, args):
        self.cfg = config
        self.method = args.method
        self.minlw = args.min_lhw
        self.jplace_fname = args.jplace_fname
        self.ranktest = args.ranktest
        self.output_fname = args.output_dir + "/" + args.output_name

        self.tmp_refaln = config.tmp_fname("%NAME%.refaln")
        self.reftree_lbl_fname = config.tmp_fname("%NAME%_lbl.tre")
        self.reftree_tax_fname = config.tmp_fname("%NAME%_tax.tre")
        self.optmod_fname = self.cfg.tmp_fname("%NAME%.opt")
        self.reftree_fname = self.cfg.tmp_fname("ref_%NAME%.tre")

        try:
            self.refjson = RefJsonParser(config.refjson_fname, ver="1.2")
        except ValueError:
            print("Invalid json file format!")
            sys.exit()
        #validate input json format 
        self.refjson.validate()
        self.rate = self.refjson.get_rate()
        self.node_height = self.refjson.get_node_height()
        self.origin_taxonomy = self.refjson.get_origin_taxonomy()
        self.bid_taxonomy_map = self.refjson.get_bid_tanomomy_map()
        self.tax_tree = self.refjson.get_tax_tree()

        reftree_str = self.refjson.get_raxml_readable_tree()
        reftree = Tree(reftree_str)
        self.reftree_size = len(reftree.get_leaves())

        # IMPORTANT: set EPA heuristic rate based on tree size!                
        self.cfg.resolve_auto_settings(self.reftree_size)
        # If we're loading the pre-optimized model, we MUST set the same rate het. mode as in the ref file        
        if self.cfg.epa_load_optmod:
            self.cfg.raxml_model = self.refjson.get_ratehet_model()
        
        self.classify_helper = TaxClassifyHelper(self.cfg, self.bid_taxonomy_map)

        self.TAXONOMY_RANKS_COUNT = 7
        self.mislabels = []
        self.mislabels_cnt = [0] * self.TAXONOMY_RANKS_COUNT
        self.rank_mislabels = []
        self.rank_mislabels_cnt = [0] * self.TAXONOMY_RANKS_COUNT
        self.misrank_conf_map = {}

    def cleanup(self):
        FileUtils.remove_if_exists(self.tmp_refaln)

    def classify_seq(self, placement):
        edges = placement["p"]
        if len(edges) > 0:
            return self.classify_helper.classify_seq(edges, self.method, self.minlw)
        else:
            print "ERROR: no placements! something is definitely wrong!"

    def rank_level_name(self, rank_level):
        return { 0: "Kingdom",
                 1: "Phylum",
                 2: "Class",
                 3: "Order",
                 4: "Family",
                 5: "Genus",
                 6: "Species"
                }[rank_level]
        
    def check_seq_tax_labels(self, seq_name, orig_ranks, ranks, lws):
        mislabel_lvl = -1
        min_len = min(len(orig_ranks),len(ranks))
        for rank_lvl in range(min_len):
            if ranks[rank_lvl] != Taxonomy.EMPTY_RANK and ranks[rank_lvl] != orig_ranks[rank_lvl]:
                mislabel_lvl = rank_lvl
                break

        if mislabel_lvl >= 0:
            mis_rec = {}
            mis_rec['name'] = seq_name.lstrip(EpacConfig.REF_SEQ_PREFIX)
            mis_rec['level'] = mislabel_lvl
            mis_rec['inv_level'] = -1 * mislabel_lvl  # just for sorting
            mis_rec['orig_ranks'] = orig_ranks
            mis_rec['ranks'] = ranks
            mis_rec['lws'] = lws
            mis_rec['conf'] = lws[mislabel_lvl]
            self.mislabels.append(mis_rec)

            for i in range(mislabel_lvl, self.TAXONOMY_RANKS_COUNT):
                self.mislabels_cnt[i] += 1
            
            return mis_rec
        else:
            return None

    def check_rank_tax_labels(self, rank_name, orig_ranks, ranks, lws):
        mislabel_lvl = -1
        min_len = min(len(orig_ranks),len(ranks))
        for rank_lvl in range(min_len):
            if ranks[rank_lvl] != Taxonomy.EMPTY_RANK and ranks[rank_lvl] != orig_ranks[rank_lvl]:
                mislabel_lvl = rank_lvl
                break

        if mislabel_lvl >= 0:
            mis_rec = {}
            mis_rec['name'] = rank_name
            mis_rec['level'] = mislabel_lvl
            mis_rec['inv_level'] = -1 * mislabel_lvl  # just for sorting
            mis_rec['orig_ranks'] = orig_ranks
            mis_rec['ranks'] = ranks
            mis_rec['lws'] = lws
            mis_rec['conf'] = lws[mislabel_lvl]
            self.rank_mislabels.append(mis_rec)

            for i in range(mislabel_lvl, self.TAXONOMY_RANKS_COUNT):
                self.rank_mislabels_cnt[i] += 1
                
            return mis_rec
        else:
            return None                

    def mis_rec_to_string_old(self, mis_rec):
        lvl = mis_rec['level']
        output = mis_rec['name'] + "\t"
        output += "%s\t%s\t%s\t%.3f\n" % (self.rank_level_name(lvl), 
            mis_rec['orig_ranks'][lvl], mis_rec['ranks'][lvl], mis_rec['lws'][lvl])
        output += ";".join(mis_rec['orig_ranks']) + "\n"
        output += ";".join(mis_rec['ranks']) + "\n"
        output += "\t".join(["%.3f" % conf for conf in mis_rec['lws']]) + "\n"
        return output

    def mis_rec_to_string(self, mis_rec):
        lvl = mis_rec['level']
        output = mis_rec['name'] + "\t"
        output += "%s\t%s\t%s\t%.3f\t" % (self.rank_level_name(lvl), 
            mis_rec['orig_ranks'][lvl], mis_rec['ranks'][lvl], mis_rec['lws'][lvl])
        output += Taxonomy.lineage_str(mis_rec['orig_ranks']) + "\t"
        output += Taxonomy.lineage_str(mis_rec['ranks']) + "\t"
        output += ";".join(["%.3f" % conf for conf in mis_rec['lws']])
        if 'rank_conf' in mis_rec:
            output += "\t%.3f" % mis_rec['rank_conf']
        return output

    def sort_mislabels(self):
        self.mislabels = sorted(self.mislabels, key=itemgetter('inv_level', 'conf'), reverse=True)
        if self.ranktest:
            self.rank_mislabels = sorted(self.rank_mislabels, key=itemgetter('inv_level', 'conf'), reverse=True)
    
    def write_mislabels(self):
        with open("%s.mis" % self.output_fname, "w") as fo_all:
            fields = ["SeqID", "MislabeledLevel", "OriginalLabel", "ProposedLabel", "Confidence", "OriginalTaxonomyPath", "ProposedTaxonomyPath", "PerRankConfidence"]
            if self.ranktest:
                fields += ["HigherRankMisplacedConfidence"]
            header = ";" + "\t".join(fields) + "\n"
            fo_all.write(header)
            if self.cfg.verbose:
                print "Mislabeled sequences:\n"
                print header 
            for mis_rec in self.mislabels:
                output = self.mis_rec_to_string(mis_rec)  + "\n"
                fo_all.write(output)
                if self.cfg.verbose:
                    print(output) 

        if self.ranktest:
            with open("%s.misrank" % self.output_fname, "w") as fo_all:
                fields = ["RankID", "MislabeledLevel", "OriginalLabel", "ProposedLabel", "Confidence", "OriginalTaxonomyPath", "ProposedTaxonomyPath", "PerRankConfidence"]
                header = ";" + "\t".join(fields)  + "\n"
                fo_all.write(header)
                if self.cfg.verbose:
                    print "\nMislabeled higher ranks:\n"
                    print header 
                for mis_rec in self.rank_mislabels:
                    output = self.mis_rec_to_string(mis_rec) + "\n"
                    fo_all.write(output)
                    if self.cfg.verbose:
                        print(output) 

        print "Mislabels counts by ranks:"        
        with open("%s.stats" % self.output_fname, "w") as fo_stat:
            for i in range(self.TAXONOMY_RANKS_COUNT):
                rname = self.rank_level_name(i).ljust(10)
                output = "%s:\t%d" % (rname, self.mislabels_cnt[i])
                if self.ranktest:
                    output += "\t%d" % self.rank_mislabels_cnt[i]
                fo_stat.write(output + "\n")
                print(output) 

    def write_sorted_map(self, fname, out_map):
        total = sum(out_map.itervalues())
        
        sorted_map = [(k,v) for v,k in reversed(sorted(
                 [(v,k) for k,v in out_map.items()]
                 ))
              ]       

        with open(fname, "w") as fout:
            for key, value in sorted_map:
                fout.write("%s: %d %f\n" % (key, value, float(value) / total))
        
    def run_leave_seq_out_test(self):
        job_name = self.cfg.subst_name("epa_%NAME%")
        if self.jplace_fname:
            jp = EpaJsonParser(self.jplace_fname)
        else:        
            jp = self.raxml.run_epa(job_name, self.refalign_fname, self.reftree_fname, self.optmod_fname, mode="l1o_seq")

        placements = jp.get_placement()
        seq_count = 0
        for place in placements:
            seq_name = place["n"][0]
            
            # get original taxonomic label
            nodes = self.tax_tree.get_leaves_by_name(seq_name)
            if len(nodes) != 1:
                print "FATAL ERROR: Sequence %s is not found in the taxonomic tree, or is present more than once!" % seq_name
                sys.exit()
            seq_node = nodes[0]
            orig_ranks = Taxonomy.split_rank_uid(seq_node.up.name)

            # get EPA tax label
            ranks, lws = self.classify_seq(place)
            # check if they match
            mis_rec = self.check_seq_tax_labels(seq_name, orig_ranks, ranks, lws)
            # cross-check with higher rank mislabels
            if self.ranktest and mis_rec:
                rank_conf = 0
                for lvl in range(2,len(orig_ranks)):
                    tax_path = Taxonomy.get_rank_uid(orig_ranks, lvl)
                    if tax_path in self.misrank_conf_map:
                        rank_conf = max(rank_conf, self.misrank_conf_map[tax_path])
                mis_rec['rank_conf'] = rank_conf
            seq_count += 1

        return seq_count    

    def run_leave_subtree_out_test(self):
        job_name = self.cfg.subst_name("tree_epa_%NAME%")
#        if self.jplace_fname:
#            jp = EpaJsonParser(self.jplace_fname)
#        else:        

        #create file with subtrees
        rank_tips = {}
        rank_parent = {}
        for node in self.tax_tree.traverse("postorder"):
            if node.is_leaf() or node.is_root():
                continue
            tax_path = node.name
            ranks = Taxonomy.split_rank_uid(tax_path)
            rank_lvl = Taxonomy.lowest_assigned_rank_level(ranks)
            if rank_lvl < 2:
                continue
                
            parent_ranks = Taxonomy.split_rank_uid(node.up.name)
            parent_lvl = Taxonomy.lowest_assigned_rank_level(parent_ranks)
            if parent_lvl < 1:
                continue
            
            rank_seqs = node.get_leaf_names()
            rank_size = len(rank_seqs)
            if rank_size < 2 or rank_size > self.reftree_size-4:
                continue

#            print rank_lvl, "\t", tax_path, "\t", rank_seqs, "\n"
            rank_tips[tax_path] = node.get_leaf_names()
            rank_parent[tax_path] = parent_ranks
                
        subtree_list = rank_tips.items()
        
        if len(subtree_list) == 0:
            return 0
            
        subtree_list_file = self.cfg.subst_name("treelist_%NAME%.txt")
        with open(subtree_list_file, "w") as fout:
            for rank_name, tips in subtree_list:
                fout.write("%s\n" % " ".join(tips))
        
        jp_list = self.raxml.run_epa(job_name, self.refalign_fname, self.reftree_fname, self.optmod_fname, 
            mode="l1o_subtree", subtree_fname=subtree_list_file)

        subtree_count = 0
        for jp in jp_list:
            placements = jp.get_placement()
            for place in placements:
                ranks, lws = self.classify_seq(place)
                tax_path = subtree_list[subtree_count][0]
                orig_ranks = Taxonomy.split_rank_uid(tax_path)
                rank_level = Taxonomy.lowest_assigned_rank_level(orig_ranks)
                rank_name = GGTaxonomyFile.add_rank_prefix(orig_ranks[rank_level], rank_level)
                parent_ranks = rank_parent[tax_path]
#                print orig_ranks, "\n", parent_ranks, "\n", ranks, "\n"
                mis_rec = self.check_rank_tax_labels(rank_name, parent_ranks, ranks, lws)
                if mis_rec:
                    self.misrank_conf_map[tax_path] = mis_rec['conf']
                subtree_count += 1

        return subtree_count    

    def run_test(self, raxml_mode = True):
        self.raxml = RaxmlWrapper(config)

        print "Total sequences: %d\n" % self.reftree_size

        self.refjson.get_raxml_readable_tree(self.reftree_fname)
        self.refalign_fname = self.refjson.get_alignment(self.tmp_refaln)        
        self.refjson.get_binary_model(self.optmod_fname)

        if self.ranktest:
            subtree_count = self.run_leave_subtree_out_test()

        seq_count = self.run_leave_seq_out_test()

        self.sort_mislabels()
        self.write_mislabels()
        print "\nPercentage of mislabeled sequences: %.2f %%" % (float(self.mislabels_cnt[self.TAXONOMY_RANKS_COUNT-1]) / seq_count * 100)

        if not self.cfg.debug:
            FileUtils.remove_if_exists(self.reftree_fname)
            FileUtils.remove_if_exists(self.optmod_fname)
            FileUtils.remove_if_exists(self.refalign_fname)

def print_options():
    print("usage: python find_mislabels.py -r example/reference.json -t 0.5 -v")
    print("Options:")
    print("    -r reference                   Specify the reference alignment and  taxonomy  in  json  format.\n")
    print("    -t min likelihood weight       A value between 0 and 1, the minimal sum of likelihood weight of")
    print("                                   an assignment to a specific rank. This value represents a confi-")
    print("                                   dence measure of the assignment,  assignments  below  this value")
    print("                                   will be discarded. Default: 0 to output all possbile assignments.\n")
    print("    -o outputdir                   Specify the directory for output files.\n")
    print("    -n name                        Specify output files prefix.\n")
    print("    -m method                      Assignment method 1 or 2")
    print("                                   1: Max sum likelihood (default)")
    print("                                   2: Max likelihood placement\n ")
    print("    -T numthread                   Specify the number of CPUs.\n")
    print("    -j jplacefile                  RAxML EPA placement file to process\n")
    print("    -v                             Print the results on screen.\n")

def parse_args():
    parser = ArgumentParser(description="Find putative mislabeled/misplaced sequences in a taxonomy.",
    epilog="Example: python find_mislabels.py -r example/reference.json -t 0.5 -v")
    parser.add_argument("-r", dest="ref_fname",
            help="""Specify the reference alignment and taxonomy in json format.""")
    parser.add_argument("-t", dest="min_lhw", type=float, default=0.,
            help="""A value between 0 and 1, the minimal sum of likelihood weight of
                    an assignment to a specific rank. This value represents a confidence 
                    measure of the assignment, assignments below this value will be discarded. 
                    Default: 0 to output all possbile assignments.""")
    parser.add_argument("-o", dest="output_dir", default=".",
            help="""Output directory""")
    parser.add_argument("-n", dest="output_name", default="result",
            help="""Query name, will be used as prefix for output file names (default: result)""")
    parser.add_argument("-m", dest="method", default="1",
            help="""Assignment method 1 or 2
                    1: Max sum likelihood (default)
                    2: Max likelihood placement""")
    parser.add_argument("-ranktest", dest="ranktest", action="store_true",
            help="""Test for misplaced higher ranks.""")
    parser.add_argument("-T", dest="num_threads", type=int, default=None,
            help="""Specify the number of CPUs.  Default: 2""")
    parser.add_argument("-v", dest="verbose", action="store_true",
            help="""Print additional info messages to the console.""")
    parser.add_argument("-debug", dest="debug", action="store_true",
            help="""Debug mode, intermediate files will not be cleaned up.""")
    parser.add_argument("-j", dest="jplace_fname", default=None,
            help="""Do not call RAxML EPA, use existing .jplace file as input instead.""")
    parser.add_argument("-c", dest="config_fname", default=None,
            help="Config file name.")
    parser.add_argument("-tmpdir", dest="temp_dir", default=None,
            help="""Directory for temporary files.""")
    args = parser.parse_args()
    return args


def check_args(args):    
    if not args.ref_fname:
        print("Must specify the reference in json format!\n")
        print_options()
        sys.exit()
    
    if not os.path.exists(args.ref_fname):
        print("Input reference json file does not exists: %s" % args.ref_fname)
        sys.exit()
    
    if not os.path.exists(args.output_dir):
        print("Output directory does not exists: %s" % args.output_dir)
        sys.exit()

    if args.jplace_fname and not os.path.exists(args.jplace_fname):
        print("EPA placement file does not exists: %s" % args.jplace_fname)
        sys.exit()

    if args.min_lhw < 0 or args.min_lhw > 1.0:
         args.min_lhw = 0.0
    
    if not (args.method == "1" or args.method == "2"):
        args.method == "1"
        
def print_run_info(config, args):
    print("Mislabels search is running with the following parameters:")
    print(" Reference:......................%s" % args.ref_fname)
    if args.jplace_fname:
        print(" EPA jplace file:................%s" % args.jplace_fname)
    print(" Min likelihood weight:..........%f" % args.min_lhw)
    print(" Assignment method:..............%s" % args.method)
    print(" Number of threads:..............%d" % config.num_threads)
    print("Result will be written to:")
    print(args.output_dir)
    print("")


if __name__ == "__main__":
    if len(sys.argv) == 1: 
        print_options()
        sys.exit()

    args = parse_args()
    check_args(args)
    config = EpacConfig(args)
    print_run_info(config, args)
   
    t = LeaveOneTest(config, args)
    t.run_test()
    if not config.debug:
        t.cleanup()

