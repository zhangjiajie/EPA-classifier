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
    from epac.msa import muscle, hmmer
    from epac.erlang import erlang
    from epac.taxonomy_util import Taxonomy,GGTaxonomyFile
except ImportError, e:
    print("Some packages are missing, please re-downloand EPA-classifier")
    print e
    sys.exit()


class LeaveOneTest:
    def __init__(self, config, args):
        self.cfg = config
        try:
            self.refjson = RefJsonParser(config.refjson_fname)
        except ValueError:
            print("Invalid json file format!")
            sys.exit()
        #validate input json format 
        self.refjson.validate()
        self.rate = self.refjson.get_rate()
        self.node_height = self.refjson.get_node_height()
        self.erlang = erlang()
        self.tmp_refaln = config.tmp_fname("%NAME%.refaln")
        self.reftree_lbl_fname = config.tmp_fname("%NAME%_lbl.tre")
        self.reftree_tax_fname = config.tmp_fname("%NAME%_tax.tre")
        self.seqs = None
        self.TAXONOMY_RANKS_COUNT = 7

        self.output_fname = args.output_dir + "/" + args.output_name
        self.method = args.method
        self.minlw = args.min_lhw
        self.jplace_fname = args.jplace_fname
        self.mislabels = []
        self.mislabels_cnt = [0] * self.TAXONOMY_RANKS_COUNT


    def cleanup(self):
        FileUtils.remove_if_exists(self.tmp_refaln)

    def get_lowest_assigned_rank_level(self, ranks):
        rank_level = len(ranks)-1
        while ranks[rank_level] == Taxonomy.EMPTY_RANK:
            rank_level -= 1
        return rank_level

    def print_ranks(self, rks, confs, minlw = 0.0):
        ss = ""
        css = ""
        for i in range(len(rks)):
            conf = confs[i]
            if conf == confs[0] and confs[0] >=0.99:
                conf = 1.0
            if conf >= minlw:
                ss = ss + rks[i] + ";"
                css = css + "{0:.3f}".format(conf) + ";"
            else:
                break
        if ss == "":
            return None
        else:
            return ss[:-1] + "\t" + css[:-1]

    def run_epa_once(self, reftree, seq_name):
        reftree_fname = self.cfg.tmp_fname("ref_%NAME%_" + seq_name + ".tre")
        job_name = self.cfg.subst_name("epa_%NAME%_" + seq_name)

        reftree.write(outfile=reftree_fname)

        epa_result = self.raxml.run_epa(job_name, self.refalign_fname, reftree_fname, self.optmod_fname)
        self.reftree_lbl_str = epa_result.get_std_newick_tree()        
        placements = epa_result.get_placement()
        
        if len(placements) != 1:
            print "ERROR: 1 placement expected, got: %d" % len(placements)
            return

        # update branchid-taxonomy mapping to account for possible changes in branch numbering
        outgr = self.reftree_outgroup.copy()
        leaf_nodes = outgr.get_leaves_by_name(seq_name)
        if len(leaf_nodes) > 0:
            leaf_nodes[0].delete()
        self.restore_rooting(outgr)
        self.label_reftree_with_ranks()
        self.build_branch_rank_map()

        ranks, lws = self.classify_seq(seq_name, placements[0])
        self.check_tax_labels(seq_name, ranks, lws)
        
        if not self.cfg.debug:
            self.raxml.cleanup(job_name)
            FileUtils.remove_if_exists(reftree_fname)
    
    def restore_rooting(self, outgr):
        self.reftree_tax = Tree(self.reftree_lbl_str)
        outgr_leaves = outgr.get_leaf_names()
        # check if outgroup consists of a single node - ETE considers it to be root, not leaf
        if not outgr_leaves:
            outgr_root = self.reftree_tax&outgr.name
        elif len(outgr_leaves) == 1:
            outgr_root = self.reftree_tax&outgr_leaves[0]
        else:
            # Even unrooted tree is "implicitely" rooted in ETE representation.
            # If this pseudo-rooting happens to be within the outgroup, it cause problems
            # in the get_common_ancestor() step (since common_ancestor = "root")
            # Workaround: explicitely root the tree outside from outgroup subtree
            for node in self.reftree_tax.iter_leaves():
                if not node.name in outgr_leaves:
                    tmp_root = node.up
                    if not self.reftree_tax == tmp_root:
                        self.reftree_tax.set_outgroup(tmp_root)
                        break
            
            outgr_root = self.reftree_tax.get_common_ancestor(outgr_leaves)

        # we could be so lucky that the RAxML tree is already correctly rooted :)
        if outgr_root != self.reftree_tax:
            self.reftree_tax.set_outgroup(outgr_root)

        if self.cfg.debug:
            #    t.show()
            self.reftree_tax.write(outfile=self.reftree_lbl_fname, format=5)

    def label_reftree_with_ranks(self):
        """labeling self.reftree_tax"""
        for node in self.reftree_tax.traverse("postorder"):
            if node.is_leaf():
                seq_ranks = self.origin_taxonomy[node.name]
                rank_level = self.get_lowest_assigned_rank_level(seq_ranks)
                node.add_feature("rank_level", rank_level)
                node.add_feature("ranks", seq_ranks)
                node.name += "__" + seq_ranks[rank_level]
#                print node.name + " -- " + ";".join(node.ranks) + " -- " + str(node.rank_level)
            else:
                if len(node.children) != 2:
                    print "FATAL ERROR: tree is not bifurcating!"
                    sys.exit()
                lchild = node.children[0]
                rchild = node.children[1]
                rank_level = min(lchild.rank_level, rchild.rank_level)
#                if hasattr(node, "B"):
#                    print node.B + " ---- " + str(rank_level) + " --- " + lchild.ranks[rank_level] + " --- " + rchild.ranks[rank_level]
                while rank_level >= 0 and lchild.ranks[rank_level] != rchild.ranks[rank_level]:
                    rank_level -= 1
                node.add_feature("rank_level", rank_level)
                node_ranks = [Taxonomy.EMPTY_RANK] * 7
                if rank_level >= 0:
                    node_ranks[0:rank_level+1] = lchild.ranks[0:rank_level+1]
                    node.name = lchild.ranks[rank_level]
                else:
                    node.name = "Undefined"
#                    print ";".join(lchild.ranks) + " -- " + ";".join(rchild.ranks) + " node_ranks: " + ";".join(node_ranks)
                    if hasattr(node, "B") and self.cfg.verbose:
                        print "INFO: no taxonomic annotation for branch %s (reason: children belong to different kingdoms)" % node.B

                node.add_feature("ranks", node_ranks)

        if self.cfg.debug:
            #    t.show()
            self.reftree_tax.write(outfile=self.reftree_tax_fname, format=3)

    def build_branch_rank_map(self):
        self.bid_taxonomy_map = {}
        for node in self.reftree_tax.traverse("postorder"):
            if not node.is_root() and hasattr(node, "B"):                
                parent = node.up                
                self.bid_taxonomy_map[node.B] = parent.ranks
#                print "%s => %s" % (node.B, parent.ranks)
#            elif self.cfg.verbose:
#                print "INFO: EPA branch label missing, mapping to taxon skipped (%s)" % node.name

    
    def classify_seq(self, placement):
        edges = placement["p"]
        if len(edges) > 0:
            if self.method == "1":
                ranks, lws = self.assign_taxonomy_maxsum(edges, self.minlw)
            else:
                ranks, lws = self.assign_taxonomy(edges)
            return ranks, lws
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
        
    def check_tax_labels(self, seq_name, ranks, lws):
        orig_ranks = self.origin_taxonomy[seq_name]
        mislabel_lvl = -1
        for rank_lvl in range(self.TAXONOMY_RANKS_COUNT):
            if ranks[rank_lvl] != Taxonomy.EMPTY_RANK and ranks[rank_lvl] != orig_ranks[rank_lvl]:
                mislabel_lvl = rank_lvl
                break

        if mislabel_lvl >= 0:
            mis_rec = {}
            mis_rec['name'] = seq_name.lstrip(EpacConfig.REF_SEQ_PREFIX)
            mis_rec['level'] = mislabel_lvl
            mis_rec['inv_level'] = -1 * mislabel_lvl  # just for sorting
            mis_rec['orig_ranks'] = GGTaxonomyFile.strip_prefix(orig_ranks)
            mis_rec['ranks'] = GGTaxonomyFile.strip_prefix(ranks)
            mis_rec['lws'] = lws
            mis_rec['conf'] = lws[mislabel_lvl]
            self.mislabels.append(mis_rec)

            for i in range(mislabel_lvl, self.TAXONOMY_RANKS_COUNT):
                self.mislabels_cnt[i] += 1

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
        return output

    def sort_mislabels(self):
        self.mislabels = sorted(self.mislabels, key=itemgetter('inv_level', 'conf'), reverse=True)
    
    def write_mislabels(self):
        with open("%s.mis" % self.output_fname, "w") as fo_all:
	    fields = ["SeqID", "MislabeledRank", "OriginalLabel", "ProposedLabel", "Confidence", "OriginalTaxonomyPath", "ProposedTaxonomyPath", "PerRankConfidence"]
	    fo_all.write(";" + "\t".join(fields) + "\n")
            for mis_rec in self.mislabels:
                output = self.mis_rec_to_string(mis_rec)            
                fo_all.write(output + "\n")
                if self.cfg.verbose:
                    print(output) 

        print "Mislabels counts by ranks:"        
        with open("%s.stats" % self.output_fname, "w") as fo_stat:
            for i in range(self.TAXONOMY_RANKS_COUNT):
                output = "%s: %d" % (self.rank_level_name(i), self.mislabels_cnt[i])            
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

    def assign_taxonomy_maxsum(self, edges, minlw = 0.):
        """this function sums up all LH-weights for each rank and takes the rank with the max. sum """
        # in EPA result, each placement(=branch) has a "weight"
        # since we are interested in taxonomic placement, we do not care about branch vs. branch comparisons,
        # but only consider rank vs. rank (e. g. G1 S1 vs. G1 S2 vs. G1)
        # Thus we accumulate weights for each rank, there are to measures:
        # "own" weight  = sum of weight of all placements EXACTLY to this rank (e.g. for G1: G1 only)
        # "total" rank  = own rank + own rank of all children (for G1: G1 or G1 S1 or G1 S2)
        rw_own = {}
        rw_total = {}
        rb = {}
        
        for edge in edges:
            br_id = str(edge[0])
            lweight = edge[2]
            lowest_rank = None

            if lweight == 0.:
                continue
            
            # accumulate weight for the current sequence                
            ranks = self.bid_taxonomy_map[br_id]
            for i in range(len(ranks)):
                rank = ranks[i]
                if rank != Taxonomy.EMPTY_RANK:
                    rw_total[rank] = rw_total.get(rank, 0) + lweight
                    lowest_rank = rank
                    if not rank in rb:
                        rb[rank] = br_id
                else:
                    break

            if lowest_rank:
                rw_own[lowest_rank] = rw_own.get(lowest_rank, 0) + lweight
                rb[lowest_rank] = br_id
            elif self.cfg.verbose:
                print "WARNING: no annotation for branch ", br_id
            

        # we assign the sequence to a rank, which has the max "own" weight AND 
        # whose "total" weight is greater than a confidence threshold
        max_rw = 0.
        s_r = None
        for r in rw_own.iterkeys():
            if rw_own[r] > max_rw and rw_total[r] >= minlw:
                s_r = r
                max_rw = rw_own[r] 
        if not s_r:
            s_r = max(rw_total.iterkeys(), key=(lambda key: rw_total[key]))

        a_br_id = rb[s_r]
        a_ranks = self.bid_taxonomy_map[a_br_id]

        # "total" weight is considered as confidence value for now
        a_conf = [0.] * len(a_ranks)
        for i in range(len(a_conf)):
            rank = a_ranks[i]
            if rank != Taxonomy.EMPTY_RANK:
                a_conf[i] = rw_total[rank]

        return a_ranks, a_conf

    def run_raxml_leave_test(self):
        job_name = self.cfg.subst_name("epa_%NAME%")
        if self.jplace_fname:
            jp = EpaJsonParser(self.jplace_fname)
        else:        
            jp = self.raxml.run_epa(job_name, self.refalign_fname, self.reftree_fname, self.optmod_fname, leave_one_test=True)

        placements = jp.get_placement()
        seq_count = 0
        for place in placements:
#            print "Placement # %d" % (seq_count + 1)
            seq_name = place["n"][0]
            ranks, lws = self.classify_seq(place)
            self.check_tax_labels(seq_name, ranks, lws)
            seq_count += 1

        return seq_count    

    def run_test(self, raxml_mode = True):
        self.raxml = RaxmlWrapper(config)
        self.refalign_fname = self.refjson.get_alignment(fout = self.tmp_refaln)        
        self.optmod_fname = self.cfg.tmp_fname("%NAME%.opt")
        self.refjson.get_binary_model(self.optmod_fname)
        self.origin_taxonomy = self.refjson.get_origin_taxonomy()
        self.orig_bid_taxonomy_map = self.refjson.get_bid_tanomomy_map()
        self.reftree_fname = self.cfg.tmp_fname("ref_%NAME%.tre")

        reftree_str = self.refjson.get_raxml_readable_tree()
        reftree = Tree(reftree_str)

        self.reftree_size = len(reftree.get_leaves())

        # IMPORTANT: set EPA heuristic rate based on tree size!                
        self.cfg.resolve_auto_settings(self.reftree_size)
        # If we're loading the pre-optimized model, we MUST set the same rate het. mode as in the ref file        
        if self.cfg.epa_load_optmod:
            self.cfg.raxml_model = self.refjson.get_ratehet_model()

        print "Total sequences: %d\n" % self.reftree_size

#        self.family_stats_fname = "%s_family.stats" % self.output_fname
#        self.cnt_family = {}
        
        if raxml_mode:
            self.refjson.get_raxml_readable_tree(self.reftree_fname)
            self.bid_taxonomy_map = self.orig_bid_taxonomy_map
            seq_count = self.run_raxml_leave_test()
            if not self.cfg.debug:
                FileUtils.remove_if_exists(self.reftree_fname)
        else:
            self.reftree_outgroup = self.refjson.get_outgroup()
            seq_count = 0
            for node in reftree.iter_leaves(): #reftree.get_leaves_by_name('r_EU338490|S002165559'):
                tmp_tree = reftree.copy() 
                leave_node = tmp_tree.get_leaves_by_name(node.name)[0]
                leave_node.delete()
                self.run_epa_once(tmp_tree, node.name)
                seq_count += 1
    #            if seq_count > 30:            
    #                break

        self.sort_mislabels()
        self.write_mislabels()
        print "\nPercentage of mislabeled sequences: %.2f %%" % (float(self.mislabels_cnt[self.TAXONOMY_RANKS_COUNT-1]) / seq_count * 100)

#        self.write_sorted_map(self.family_stats_fname, self.cnt_family)

        if not self.cfg.debug:
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
    parser = ArgumentParser(description="Find putative mislabeled/misplaced sequences in a taxonomy.")
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
    args = parser.parse_args()
    return args


def check_args(args):    
    if not args.ref_fname:
        print("Must specify the reference in json format!\n")
        print_options()
        sys.exit()
    
    if not os.path.exists(args.ref_fname):
        print("Input reference json file does not exists: %s" % args.ref_fname)
        print_options()
        sys.exit()
    
    if not os.path.exists(args.output_dir):
        print("Output directory does not exists: %s" % args.output_dir)
        print_options()
        sys.exit()

    if args.jplace_fname and not os.path.exists(args.jplace_fname):
        print("EPA placement file does not exists: %s" % args.jplace_fname)
        print_options()
        sys.exit()

    if args.min_lhw < 0 or args.min_lhw > 1.0:
         args.min_lhw = 0.0
    
    if not (args.method == "1" or args.method == "2"):
        args.method == "1"
        
def print_run_info(config, args):
    print("EPA-classifier running with the following parameters:")
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

