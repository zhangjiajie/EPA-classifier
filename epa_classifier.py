#! /usr/bin/env python
try:
    import sys
    import os
    import json
    import time
    import glob
    from epac.ete2 import Tree, SeqGroup
    from epac.argparse import ArgumentParser
    from epac.epa_util import epa
    from epac.json_util import jsonparser, json_checker
    from epac.msa import muscle, hmmer
    from epac.erlang import erlang
    from epac.taxonomy_util import TaxonomyUtils
except ImportError:
    print("Some packages are missing, please re-downloand EPA-classifier")
    sys.exit()


class magic:
    def __init__(self, refjson, query, verbose = True, numcpu = "2"):
        self.numcpus = numcpu
        self.v = verbose
        self.refjson = jsonparser(refjson)
        self.refjson.validate()
        self.bid_taxonomy_map = self.refjson.get_bid_tanomomy_map()
        self.reftree = self.refjson.get_reftree()
        self.rate = self.refjson.get_rate()
        self.node_height = self.refjson.get_node_height()
        self.query = query
        self.basepath = os.path.dirname(os.path.abspath(__file__))
        self.erlang = erlang()
        self.tmppath = self.basepath + "/epac/tmp"
        self.name = str(time.time())
        self.tmp_refaln = self.tmppath + "/" + self.name + ".refaln"
        self.epa_alignment = self.tmppath + "/" + self.name + ".afa"
        self.hmmprofile = self.tmppath + "/" + self.name + ".hmmprofile"
        self.tmpquery = self.tmppath + "/" + self.name + ".tmpquery"
        self.min_confidence=0.2
        self.seqs = None


    def cleanup(self):
        self.remove(self.tmp_refaln)
        self.remove(self.epa_alignment)
        self.remove(self.hmmprofile)
        self.remove(self.tmpquery)
        reduced = glob.glob(self.tmppath + "/*.reduced")
        for f in reduced:
            self.remove(f)
        reduced = glob.glob(self.tmppath + "/*.jplace")
        for f in reduced:
            self.remove(f)


    def remove(self, filename):
        if os.path.exists(filename):
            os.remove(filename)


    def align_to_refenence(self, noalign, minp = 0.9):
        self.refjson.get_hmm_profile(self.hmmprofile)
        refaln = self.refjson.get_alignment(fout = self.tmp_refaln)
        hm = hmmer(refalign = refaln , query = self.tmpquery, refprofile = self.hmmprofile, discard = noalign, seqs = self.seqs, minp = minp)
        self.epa_alignment = hm.align()


    def merge_alignment(self, query_seqs):
        refaln = self.refjson.get_alignment_list()
        with open(self.epa_alignment, "w") as fout:
            for seq in refaln:
                fout.write(">" + seq[0] + "\n" + seq[1] + "\n")
            for name, seq, comment, sid in query_seqs.iter_entries():
                fout.write(">" + str(sid) + "\n" + seq + "\n")


    def checkinput(self, noalignpath, minp = 0.9):
        self.seqs = SeqGroup(sequences=self.query, format = "fasta")
        self.seqs.write(format="fasta_internal", outfile=self.tmpquery)
        print("Checking if query sequences are aligned ...")
        entries = self.seqs.get_entries()
        seql = len(entries[0][1])
        aligned = True
        for entri in entries[1:]:
            l = len(entri[1])
            if not seql == l:
                aligned = False
                break
        
        if aligned and len(self.seqs) > 1:
            print("Query sequences are aligned")
            refalnl = self.refjson.get_alignment_length()
            if refalnl == seql:
                print("Merging query alignment with reference alignment")
                self.merge_alignment(self.seqs)
            else:
                print("Merging query alignment with reference alignment using MUSCLE")
                require_muscle()
                refaln = self.refjson.get_alignment(fout = self.tmp_refaln)
                m = muscle()
                self.epa_alignment = m.merge(refaln, self.tmpquery)
        else:
            print("Query sequences are not aligned")
            print("Align query sequences to the reference alignment using HMMER")
            require_hmmer()
            self.align_to_refenence(noalignpath, minp = minp)
        
        print("Running EPA ......")


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


    def classify(self, fout = None, fnoalign = None, method = "1", minlw = 0.0, pv = 0.02, minp = 0.9):
        self.checkinput(fnoalign, minp = minp)
        EPA = epa()
        placements = EPA.run(reftree = self.refjson.get_raxml_readable_tree(), alignment = self.epa_alignment, num_thread = self.numcpus).get_placement()
        EPA.clean()
        if fout!=None:
            fo = open(fout, "w")
        
        for place in placements:
            taxa_name = place["n"][0]
            origin_taxa_name = self.seqs.get_name(int(taxa_name))
            edges = place["p"]
            edges = self.erlang_filter(edges, p = pv)
            if len(edges) > 0:
                if method == "1":
                    ranks, lws = self.assign_taxonomy_maxsum(edges)
                else:
                    ranks, lws = self.assign_taxonomy(edges)
                
                isnovo = self.novelty_check(place_edge = str(edges[0][0]), ranks =ranks, lws = lws, minlw = minlw)
                rankout = self.print_ranks(ranks, lws, minlw)
                if rankout == None:
                    output = origin_taxa_name+ "\t\t\t?"
                else:
                    if isnovo: 
                        output = origin_taxa_name+ "\t" + self.print_ranks(ranks, lws, minlw) + "\t*"
                    else:
                        output = origin_taxa_name+ "\t" + self.print_ranks(ranks, lws, minlw) + "\to"
                if self.v:
                    print(output) 
                if fout!=None:
                    fo.write(output + "\n")
            else:
                output = origin_taxa_name+ "\t\t\t?"
                if self.v:
                    print(output) 
                if fout!=None:
                    fo.write(output + "\n")
        
        if fout!=None:
            fo.close()


    def erlang_filter(self, edges, p = 0.02):
        newedges = []
        for edge in edges:
            edge_nr = str(edge[0])
            pendant_length = edge[4]
            pv = self.erlang.one_tail_test(rate = self.rate, k = int(self.node_height[edge_nr]), x = pendant_length)
            if pv >= p:
                newedges.append(edge)
        return newedges


    def novelty_check(self, place_edge, ranks, lws, minlw):
        """If the taxonomic assignment is not assigned to the genus level, 
        we need to check if it is due to the incomplete reference taxonomy or 
        it is likely to be something new:
        
        1. If the final ranks are assinged because of lw cut, that means with samller lw
        the ranks can be further assinged to lowers. This indicate the undetermined ranks 
        in the assignment is not due to the incomplete reference taxonomy, so the query 
        sequence is likely to be something new.
        
        2. Otherwise We check all leaf nodes' immediate lower rank below this ml placement point, 
        if they are not empty, output all ranks and indicate this could be novelty.
        """
        
        lowrank = 0
        for i in range(len(ranks)):
            if i < 6:
                """above genus level"""
                rk = ranks[i]
                lw = lws[i]
                if rk == "-":
                    break
                else:
                    lowrank = lowrank + 1
                    if lw >=0 and lw < minlw:
                        return True
        
        if lowrank >= 5 and not ranks[lowrank] == "-":
            return False
        else:
            placenode = self.reftree.search_nodes(B = place_edge)[0]
            if placenode.is_leaf():
                return False
            else:
                leafnodes = placenode.get_leaves()
                flag = True
                for leaf in leafnodes:
                    br_num = leaf.B
                    branks = self.bid_taxonomy_map[br_num]
                    if branks[lowrank] == "-":
                        flag = False
                        break
                        
                return flag


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


    def assign_taxonomy_maxsum(self, edges):
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
            
            # accumulate weight for the current sequence                
            ranks = self.bid_taxonomy_map[br_id]
            for i in range(len(ranks)):
                rank = ranks[i]
                if rank != TaxonomyUtils.EMPTY_RANK:
                    rw_total[rank] = rw_total.get(rank, 0) + lweight
                    lowest_rank = rank
                    if not rank in rb:
                        rb[rank] = br_id
                else:
                    break
            
            rw_own[lowest_rank] = rw_own.get(lowest_rank, 0) + lweight
            rb[lowest_rank] = br_id

        # we assign the sequence to a rank, which has the max "own" weight AND 
        # whose "total" weight is greater than a confidence threshold
        max_rw = 0.
        for r in rw_own.iterkeys():
            if rw_own[r] > max_rw and rw_total[r] >= self.min_confidence:
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
            if rank != TaxonomyUtils.EMPTY_RANK:
                a_conf[i] = rw_total[rank]

        return a_ranks, a_conf


def print_options():
    print("usage: python epa_classifier.py -r example/reference.json -q example/query.fa -t 0.5 -v")
    print("Options:")
    print("    -r reference                   Specify the reference alignment and  taxonomy  in  json  format.\n")
    print("    -q query sequence              Specify the query seqeunces file.")
    print("                                   If the query sequences are aligned to  the  reference  alignment")
    print("                                   already, the epa_classifier will  classify the  queries  to  the")
    print("                                   lowest rank possible. If the  query sequences  are aligned,  but")
    print("                                   not to the reference, then a profile alignment will be perfermed")
    print("                                   to merge the two alignments. If  the  query  sequences  are  not")
    print("                                   aligned, then HMMER will be used to align  the  queries  to  the")
    print("                                   reference alignment. \n")
    print("    -t min likelihood weight       A value between 0 and 1, the minimal sum of likelihood weight of")
    print("                                   an assignment to a specific rank. This value represents a confi-")
    print("                                   dence measure of the assignment,  assignments  below  this value")
    print("                                   will be discarded. Default: 0 to output all possbile assignments.\n")
    print("    -o outputfile                  Specify the file name for output.\n")
    print("    -p p-value                     P-value for Erlang test.  Default: 0.02\n")
    print("    -minalign min-aligned%         Minimal percent of aligned sites.  Default: 0.9 (90%)\n")
    print("    -m method                      Assignment method 1 or 2")
    print("                                   1: Max sum likelihood (default)")
    print("                                   2: Max likelihood placement\n ")
    print("    -T numthread                   Specify the number of CPUs.\n")
    print("    -v                             Print the results on screen.\n")


def require_muscle():
    basepath = os.path.dirname(os.path.abspath(__file__))
    if not os.path.exists(basepath + "/epac/bin/muscle"):
        print("The pipeline uses MUSCLE to merge alignment,")
        print("please download the programm from:")
        print("http://www.drive5.com/muscle/downloads.htm")
        print("Rename the executable to usearch and put it to epac/bin/  \n")
        sys.exit() 


def require_hmmer():
    basepath = os.path.dirname(os.path.abspath(__file__))
    if not os.path.exists(basepath + "/epac/bin/hmmbuild") or not os.path.exists(basepath + "/epac/bin/hmmalign"):
        print("The pipeline uses HAMMER to align the query seqeunces,")
        print("please download the programm from:")
        print("http://hmmer.janelia.org/")
        print("Copy the executables hmmbuild and hmmalign to bin/  \n")
        sys.exit()


def parse_args():
    parser = ArgumentParser(description="Assign taxonomy ranks to query sequences")
    parser.add_argument("-r", dest="ref_fname",
            help="""Specify the reference alignment and  taxonomy in json  format.""")
    parser.add_argument("-q", dest="query_fname",
            help="""Specify the query seqeunces file.\n
                    If the query sequences are aligned to  the  reference  alignment
                    already, the epa_classifier will  classify the  queries  to  the
                    lowest rank possible. If the  query sequences  are aligned,  but
                    not to the reference, then a profile alignment will be performed
                    to merge the two alignments. If  the  query  sequences  are  not
                    aligned, then HMMER will be used to align  the  queries  to  the
                    reference alignment.""")
    parser.add_argument("-t", dest="min_lhw", type=float, default=0.,
            help="""A value between 0 and 1, the minimal sum of likelihood weight of
                    an assignment to a specific rank. This value represents a confidence 
                    measure of the assignment,  assignments  below  this value will be discarded. 
                    Default: 0 to output all possbile assignments.""")
    parser.add_argument("-o", dest="output_fname", default="",
            help="""Query name, will be used as a name of results folder (default: <reftree>_YYYYMMDD_HHMM)""")
    parser.add_argument("-p", dest="p_value", type=float, default=0.02,
            help="""P-value for Erlang test.  Default: 0.02""")
    parser.add_argument("-minalign", dest="minalign", type=float, default=0.9,
            help="""Minimal percent of the sites aligned to the reference alignment.  Default: 0.9""")
    parser.add_argument("-m", dest="method", default="1",
            help="""Assignment method 1 or 2
                    1: Max sum likelihood (default)
                    2: Max likelihood placement""")
    parser.add_argument("-T", dest="num_threads", type=int, default=2,
            help="""Specify the number of CPUs.  Default: 2""")
    parser.add_argument("-v", dest="verbose", action="store_true",
            help="""Print the results on screen.""")
    parser.add_argument("-debug", dest="debug", action="store_true",
            help="""debug model, intermediate files will not be cleaned up.""")
    parser.add_argument("-a", dest="align_only", action="store_true",
            help="""Alignment only: Do not perform classification, just build the combined alignment (RS+QS) (default: OFF)""")
    parser.add_argument("-j", dest="jplace_fname",
            help="""Do not call RAxML EPA, use existing .jplace file as input instead.""")
    parser.add_argument("-c", dest="config_fname", default="epac.cfg",
            help="config file name (default: epac.cfg)")
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
    
    if not args.query_fname:
        print("The query can not be empty!\n")
        print_options()
        sys.exit()
    
    if not os.path.exists(args.query_fname):
        print("Input query file does not exists: %s" % args.query_fname)
        print_options()
        sys.exit()
    
    if args.output_fname == "":
        args.output_fname = args.query_fname + ".assignment.txt"
        args.noalign_fname = args.query_fname + ".noalign"
    else:
        args.noalign_fname = args.output_fname + ".noalign"
    
    if args.min_lhw < 0 or args.min_lhw > 1.0:
         args.min_lhw = 0.0
    
    if args.p_value < 0:
        args.p_value = 0
    
    if not (args.method == "1" or args.method == "2"):
        args.method == "1"
    
    print("EPA-classifier running with the following parameters:")
    print(" Reference:.....................%s" % args.ref_fname)
    print(" Query:.........................%s" % args.query_fname)
    print(" Min likelihood weight:.........%f" % args.min_lhw)
    print(" Assignment method:.............%s" % args.method)
    print(" P-value for Erlang test:.......%f" % args.p_value)
    print(" Number of threads:.............%d" % args.num_threads)
    print("Result will be written to:")
    print(args.output_fname)
    print("Sequences which can not be aligned will be written to:")
    print(args.noalign_fname)
    print("")


if __name__ == "__main__":
    if len(sys.argv) == 1: 
        print_options()
        sys.exit()

    args = parse_args()
    check_args(args)
   
    m = magic(refjson = args.ref_fname, query = args.query_fname, verbose = args.verbose, numcpu = args.num_threads)
    m.classify(fout = args.output_fname, fnoalign = args.noalign_fname, method = args.method, minlw = args.min_lhw, pv = args.p_value, minp = args.minalign)
    if not args.debug:
        m.cleanup()


    
