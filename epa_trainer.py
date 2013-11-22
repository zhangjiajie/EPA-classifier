#!/usr/bin/env python

import sys
import os
import shutil
from epac.ete2 import Tree, SeqGroup
from epac.argparse import ArgumentParser,RawTextHelpFormatter
from epac.config import EpacTrainerConfig
from epac.raxml_util import RaxmlWrapper, FileUtils
from epac.taxonomy_util import Taxonomy, GGTaxonomyFile, TaxTreeBuilder
from epac.json_util import RefJsonBuilder
from epac.erlang import tree_param 
from epac.msa import hmmer

class RefTreeBuilder:
    # this prefix will be added to every sequence name in reference to prevent 
    # name clashes with query sequences, which are coded with numbers
    REF_SEQ_PREFIX="r_"

    def __init__(self, config): 
        self.cfg = config
        self.mfresolv_job_name = "mfresolv"
        self.epalbl_job_name = "epalbl"
        self.optmod_job_name = "optmod"
        self.raxml_wrapper = RaxmlWrapper(config)
        self.taxonomy = GGTaxonomyFile(config.taxonomy_fname, RefTreeBuilder.REF_SEQ_PREFIX)
        
        self.outgr_fname = self.cfg.tmp_fname("%NAME%_outgr.tre")
        self.reftree_mfu_fname = self.cfg.tmp_fname("%NAME%_mfu.tre")
        self.reftree_bfu_fname = self.cfg.tmp_fname("%NAME%_bfu.tre")
        self.optmod_fname = self.cfg.tmp_fname("%NAME%.opt")
        self.lblalign_fname = self.cfg.tmp_fname("%NAME%_lblq.fa")
        self.reftree_lbl_fname = self.cfg.tmp_fname("%NAME%_lbl.tre")
        self.reftree_tax_fname = self.cfg.tmp_fname("%NAME%_tax.tre")
        self.brmap_fname = self.cfg.tmp_fname("%NAME%_map.txt")

    def build_multif_tree(self):
        c = self.cfg
        tb = TaxTreeBuilder(c, self.taxonomy)
        (t, ids) = tb.build(c.reftree_min_rank, c.reftree_max_seqs_per_leaf, c.reftree_clades_to_include, c.reftree_clades_to_ignore)
        self.reftree_ids = frozenset(ids)
        self.reftree_size = len(ids)
        self.reftree_multif = t

        if self.cfg.debug:
            refseq_fname = self.cfg.tmp_fname("%NAME%_seq_ids.txt")
            # list of sequence ids which comprise the reference tree
            with open(refseq_fname, "w") as f:
                for sid in ids:
                    f.write("%s\n" % sid)

            # original tree with taxonomic ranks as internal node labels
            reftax_fname = self.cfg.tmp_fname("%NAME%_mfu_tax.tre")
            t.write(outfile=reftax_fname, format=8)
        #    t.show()

    def export_ref_alignment(self):
        """This function transforms the input alignment in the following way:
           1. Filter out sequences which are not part of the reference tree
           2. Add sequence name prefix (r_)"""
        self.refalign_fname = self.cfg.tmp_fname("%NAME%_matrix.afa")
        with open(self.cfg.align_fname, "r") as fin:
            with open(self.refalign_fname, "w") as fout:
                while True:
                    line = fin.readline()
                    if not line: 
                        break
                    sid = RefTreeBuilder.REF_SEQ_PREFIX + line.strip()[1:]
                    if sid in self.reftree_ids:
                        fout.write(">" +  sid + "\n") 
                        fout.write(fin.readline())
                    else:
                        fin.readline()    
    
    def export_ref_taxonomy(self):
        self.taxonomy_map = {}
        
        for sid, ranks in self.taxonomy.iteritems():
            if sid in self.reftree_ids:
                self.taxonomy_map[sid] = ranks
            
        if self.cfg.debug:
            tax_fname = self.cfg.tmp_fname("%NAME%_tax.txt")
            with open(tax_fname, "w") as fout:
                for sid, ranks in self.taxonomy_map.iteritems():
                    ranks_str = self.taxonomy.lineage_str(sid, True) 
                    fout.write(sid + "\t" + ranks_str + "\n")   

    def save_rooting(self):
        rt = self.reftree_multif

        # remove unifurcation at the root
        if len(rt.children) == 1:
            rt = rt.children[0]

        if len(rt.children) > 1:
            outgr = rt.children[0]    
            outgr_size = len(outgr.get_descendants())
            for child in rt.children:
                if child != outgr:
                    child_size = len(child.get_descendants())
                    if child_size < outgr_size:
                        outgr = child
                        outgr_size = child_size
        else:
            print "Invalid tree: unifurcation at the root node!"
            sys.exit()

        
        self.reftree_outgroup = outgr
        outgr.write(outfile=self.outgr_fname, format=9)
        if self.cfg.verbose:
            print "Outgroup for rooting was saved to: " + self.outgr_fname + ", outgroup size: " + str(outgr_size+1)
        
        # now we can safely unroot the tree and remove internal node labels to make it suitable for raxml
        rt.write(outfile=self.reftree_mfu_fname, format=9)

    # RAxML call to convert multifurcating tree to the strictly bifurcating one
    def resolve_multif(self):
        print "\nReducing the alignment: \n"
        self.reduced_refalign_fname = self.raxml_wrapper.reduce_alignment(self.refalign_fname)
        
        # hack: use CAT for large trees, since GAMMA has numerical problems with them
        if self.reftree_size > 10000 and self.cfg.raxml_model != "GTRCAT":
            orig_raxml_model = self.cfg.raxml_model
            self.cfg.raxml_model = "GTRCAT"
        else:
            orig_raxml_model = ""

        print "\nResolving multifurcation: \n"
        raxml_params = ["-s", self.reduced_refalign_fname, "-g", self.reftree_mfu_fname]
        self.raxml_wrapper.run(self.mfresolv_job_name, raxml_params)
        if self.raxml_wrapper.result_exists(self.mfresolv_job_name):        
            self.raxml_wrapper.copy_result_tree(self.mfresolv_job_name, self.reftree_bfu_fname)
            if not self.cfg.debug:
                self.raxml_wrapper.cleanup(self.mfresolv_job_name)
        
            # restore the original model (s. above)
            if orig_raxml_model:
                self.cfg.raxml_model = orig_raxml_model
            
            # RAxML call to optimize model parameters and write them down to the binary model file
            print "\nOptimizing model parameters under GTR+GAMMA: \n"
            raxml_params = ["-f", "e", "-s", self.reduced_refalign_fname, "-t", self.reftree_bfu_fname]
            self.raxml_wrapper.run(self.optmod_job_name, raxml_params)
            if self.raxml_wrapper.result_exists(self.optmod_job_name):
                self.raxml_wrapper.copy_optmod_params(self.optmod_job_name, self.optmod_fname)
                if not self.cfg.debug:
                    self.raxml_wrapper.cleanup(self.optmod_job_name)
            else:
                print "RAxML run failed (model optimization), please examine the log for details: %s" \
                        % self.raxml_wrapper.make_raxml_fname("output", self.optmod_fname)
                sys.exit()  
        else:
            print "RAxML run failed (mutlifurcation resolution), please examine the log for details: %s" \
                    % self.raxml_wrapper.make_raxml_fname("output", self.mfresolv_job_name)
            sys.exit()  
            
    def load_reduced_refalign(self):
        formats = ["fasta", "phylip_relaxed"]
        for fmt in formats:
            try:
                self.reduced_refalign_seqs = SeqGroup(sequences=self.reduced_refalign_fname, format = fmt)
                break
            except:
                pass
        if self.reduced_refalign_seqs == None:
            print("FATAL ERROR: Invalid input file format in %s! (load_reduced_refalign)" % self.reduced_refalign_fname)
            sys.exit()
    
    # dummy EPA run to label the branches of the reference tree, which we need to build a mapping to tax ranks    
    def epa_branch_labeling(self):
        # create alignment with dummy query seq
        seqlen = len(self.reduced_refalign_seqs.get_seqbyid(0))
        self.reduced_refalign_seqs.write(format="fasta", outfile=self.lblalign_fname)
        
        with open(self.lblalign_fname, "a") as fout:
            fout.write(">" + "DUMMY131313" + "\n")        
            fout.write("A"*seqlen + "\n")        
        
        epa_result = self.raxml_wrapper.run_epa(self.epalbl_job_name, self.lblalign_fname, self.reftree_bfu_fname, self.optmod_fname)
        self.reftree_lbl_str = epa_result.get_std_newick_tree()
        if self.raxml_wrapper.epa_result_exists(self.epalbl_job_name):        
            if not self.cfg.debug:
                self.raxml_wrapper.cleanup(self.epalbl_job_name)
        else:
            print "RAxML EPA run failed, please examine the log for details: %s" \
                    % self.raxml_wrapper.make_raxml_fname("output", self.epalbl_job_name)
            sys.exit()        

    def restore_rooting(self):
        self.reftree_tax = Tree(self.reftree_lbl_str)
        outgr = self.reftree_outgroup
        outgr_leafs = outgr.get_leaf_names()
        # check if outgroup consists of a single node - ETE considers it to be root, not leaf
        if not outgr_leafs:
            outgr_leafs = [outgr.name]
        if len(outgr_leafs) > 1:
            outgr_root = self.reftree_tax.get_common_ancestor(outgr_leafs)
        else:
            outgr_root = self.reftree_tax&outgr_leafs[0]
        self.reftree_tax.set_outgroup(outgr_root)

        if self.cfg.debug:
            #    t.show()
            self.reftree_tax.write(outfile=self.reftree_lbl_fname, format=5)

    def label_reftree_with_ranks(self):
        """labeling self.reftree_tax"""
        for node in self.reftree_tax.traverse("postorder"):
            if node.is_leaf():
                # assume that seqs are always assigned to species level (7th rank, 0-based index = 6)           
                seq_ranks = self.taxonomy.get_seq_ranks(node.name)
                rank_level = self.taxonomy.lowest_assigned_rank_level(node.name)
                node.add_feature("rank_level", rank_level)
                node.add_feature("ranks", seq_ranks)
                node.name += "__" + seq_ranks[rank_level]
     #           print node.B
            else:
                if len(node.children) != 2:
                    print "FATAL ERROR: tree is not bifurcating!"
                    sys.exit()
                lchild = node.children[0]
                rchild = node.children[1]
                rank_level = min(lchild.rank_level, rchild.rank_level)
    #            print str(rank_level) + "---" + lchild.ranks[rank_level] + " --- " + lchild.ranks[rank_level]
                while rank_level >= 0 and (lchild.ranks[rank_level] != rchild.ranks[rank_level]):
                    rank_level -= 1
                if (rank_level >= 0):
                    node.add_feature("rank_level", rank_level)
                    node.name = lchild.ranks[rank_level]
                    node_ranks = [Taxonomy.EMPTY_RANK] * 7
                    node_ranks[0:rank_level+1] = lchild.ranks[0:rank_level+1]
    #                print "Node rank: " + str(rank_level)
    #                print node_ranks
                    node.add_feature("ranks", node_ranks)
    #                print node_ranks
                else:
                    print "FATAL ERROR: sequences belong to different kingdoms! (or something went really wrong...)"
                    sys.exit()

        if self.cfg.debug:
            #    t.show()
            self.reftree_tax.write(outfile=self.reftree_tax_fname, format=3)

    def build_branch_rank_map(self):
        self.bid_ranks_map = {}
        for node in self.reftree_tax.traverse("postorder"):
            if not node.is_root() and hasattr(node, "B"):                
                parent = node.up                
                self.bid_ranks_map[node.B] = parent.ranks
            elif self.cfg.verbose:
                print "INFO: EPA branch label missing, mapping to taxon skipped (%s)" % node.name
    
    def write_branch_rank_map(self):
        with open(self.brmap_fname, "w") as fbrmap:    
            for node in self.reftree_tax.traverse("postorder"):
                if not node.is_root() and hasattr(node, "B"):                
                    fbrmap.write(node.B + "\t" + ";".join(self.bid_ranks_map[node.B]) + "\n")
    
    def calc_node_heights(self):
        """Calculate node heights on the reference tree (used to define branch-length cutoff during classification step)
           Algorithm is as follows:
           Tip node or node resolved to Species level: height = 1 
           Inner node resolved to Genus or above:      height = min(left_height, right_height) + 1 
         """
        nh_map = {}
        dummy_added = False
        for node in self.reftree_tax.traverse("postorder"):
            if not node.is_root():
                if not hasattr(node, "B"):                
                    # In a rooted tree, there is always one more node/branch than in unrooted one
                    # That's why one branch will be always not EPA-labelled after the rooting
                    if not dummy_added: 
                        node.B = "DDD"
                        dummy_added = True
                        species_rank = Taxonomy.EMPTY_RANK
                    else:
                        print "FATAL ERROR: More than one tree branch without EPA label (calc_node_heights)"
                        sys.exit()
                else:
                    species_rank = self.bid_ranks_map[node.B][6]
                bid = node.B
                if node.is_leaf() or species_rank != Taxonomy.EMPTY_RANK:
                    nh_map[bid] = 1
                else:
                    lchild = node.children[0]
                    rchild = node.children[1]
                    nh_map[bid] = min(nh_map[lchild.B], nh_map[rchild.B]) + 1

        # remove heights for dummy nodes, since there won't be any placements on them
        if dummy_added:
            del nh_map["DDD"]
            
        self.node_height_map = nh_map

    def __get_all_rank_names(self, root):
        rnames = set([])
        for node in root.traverse("postorder"):
            ranks = node.ranks
            for rk in ranks:
                rnames.add(rk)
        return rnames

    def mono_index(self):
        """This method will calculate monophyly index by looking at the left and right hand side of the tree"""
        children = self.reftree_tax.children
        if len(children) == 1:
            while len(children) == 1:
                children = children[0].children 
        if len(children) == 2:
            left = children[0]
            right =children[1]
            lset = self.__get_all_rank_names(left)
            rset = self.__get_all_rank_names(right)
            iset = lset & rset
            return iset
        else:
            print("Error: input tree not birfurcating")
            return set([])

    def write_json(self):
        jw = RefJsonBuilder()

        jw.set_taxonomy(self.bid_ranks_map)
        jw.set_tree(self.reftree_lbl_str)

        seqs = self.reduced_refalign_seqs.get_entries()    
        jw.set_sequences(seqs)
        
        # this stupid workaround is needed because RAxML outputs the reduced
        # alignment in relaxed PHYLIP format, which is not supported by HMMER
        refalign_fasta = self.cfg.tmp_fname("%NAME%_ref_reduced.fa")
        self.reduced_refalign_seqs.write(outfile=refalign_fasta)

        print "Building the HMMER profile...\n"

        try:
            hmm = hmmer(self.cfg, refalign_fasta)
            fprofile = hmm.build_hmm_profile()
            jw.set_hmm_profile(fprofile)
        finally:
            FileUtils.remove_if_exists(refalign_fasta)
            FileUtils.remove_if_exists(fprofile)
 
        orig_tax = self.taxonomy_map
        jw.set_origin_taxonomy(orig_tax)
        
        print "Calculating the speciation rate...\n"
        tp = tree_param(tree = self.reftree_lbl_str, origin_taxonomy = orig_tax)
        jw.set_rate(tp.get_speciation_rate_fast())
        jw.set_nodes_height(self.node_height_map)
        
        jw.set_binary_model(self.optmod_fname)
        
        print "Writing down the reference file...\n"
        jw.dump(self.cfg.refjson_fname)

    def cleanup(self):
        FileUtils.remove_if_exists(self.outgr_fname)
        FileUtils.remove_if_exists(self.reftree_mfu_fname)
        FileUtils.remove_if_exists(self.reftree_bfu_fname)
        FileUtils.remove_if_exists(self.optmod_fname)
        FileUtils.remove_if_exists(self.lblalign_fname)
        FileUtils.remove_if_exists(self.outgr_fname)
        FileUtils.remove_if_exists(self.reduced_refalign_fname)
        FileUtils.remove_if_exists(self.refalign_fname)

    # top-level function to build a reference tree    
    def build_ref_tree(self):
        print "\n=> Building a multifurcating tree from taxonomy: " + self.cfg.taxonomy_fname + "...\n"
        self.build_multif_tree()
#        sys.exit()        
        print "\n==> Building the reference alignment " + "...\n"
        self.export_ref_alignment()
        self.export_ref_taxonomy()
        print "\n===> Saving the outgroup for later re-rooting " + "...\n"
        self.save_rooting()
        print "\n====> RAxML call: resolve multifurcation " + "...\n"
        self.resolve_multif()
        self.load_reduced_refalign()
        print "\n=====> RAxML-EPA call: labeling the branches " + "...\n"
        self.epa_branch_labeling()
        print "\n======> Re-rooting the reference tree" + "...\n"
        self.restore_rooting()
        print "\n=======> Labeling the reference tree with taxonomic ranks" + "...\n"
        self.label_reftree_with_ranks()
        print "\n========> Building the mapping between EPA branch labels and ranks" + "...\n"
        self.build_branch_rank_map()
        self.calc_node_heights()
        
        print "\n=========> Checking branch labels ...\n"
        print "shared rank names before trainning: " + repr(self.taxonomy.get_common_ranks())
        print "shared rank names after  trainning: " + repr(self.mono_index())
        
        print "\n=========> Saving the reference JSON file" + "...\n"
        self.write_json()
        print "\n***********  Done!  **********\n"

def parse_args():
    parser = ArgumentParser(description="Build a reference tree for EPA taxonomic placement.",
    epilog="Example: ./epa_trainer.py -t example/training_tax.txt -s example/training_seq.fa -r example/ref.json",
    formatter_class=RawTextHelpFormatter)
    parser.add_argument("-t", dest="taxonomy_fname",
            help="""Reference taxonomy file.""")
    parser.add_argument("-s", dest="align_fname",
            help="""Reference alignment file. Sequences must be aligned, their IDs must correspond to those
in taxonomy file.""")
    parser.add_argument("-r", dest="ref_fname",
            help="""Reference output file. It will contain reference alignment, phylogenetic tree and other
information needed for taxonomic placement of query sequences.""")
    parser.add_argument("-T", dest="num_threads", type=int, default=None,
            help="""Specify the number of CPUs.  Default: 2""")            
    parser.add_argument("-v", dest="verbose", action="store_true",
            help="""Print additional info messages to the console.""")
    parser.add_argument("-debug", dest="debug", action="store_true",
            help="""Debug mode, intermediate files will not be cleaned up.""")
    parser.add_argument("-c", dest="config_fname", default=None,
            help="""Config file name.""")
    
    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()
    
    return args
 
def check_args(args):
    #check if taxonomy file exists
    if not os.path.isfile(args.taxonomy_fname):
        print "ERROR: Taxonomy file not found: %s" % args.tax_fname
        sys.exit()

    #check if alignment file exists
    if not os.path.isfile(args.align_fname):
        print "ERROR: Alignment file not found: %s" % args.align_fname
        sys.exit()

    #check if reference tree file already exists
    if os.path.isfile(args.ref_fname):
        print "ERROR: Reference tree file already exists: %s" % args.ref_fname
        print "Please delete it explicitely if you want to overwrite."
        sys.exit()

# -------
# MAIN
# -------
if __name__ == "__main__":
    args = parse_args()
    check_args(args)
    config = EpacTrainerConfig(args)
    builder = RefTreeBuilder(config)
    builder.build_ref_tree()
    if not args.debug:
        builder.cleanup()
