#! /usr/bin/env python
import sys
import os
import json
import operator
from coretype.seqgroup import SeqGroup
from coretype.tree import Tree
from subprocess import call
from taxonomy_parser import TreeBuilder
from epa_util import epa, raxml
from json_util import jsonparser
from msa import hmmer
from erlang import tree_param 


class sequence:
    """Provide sequence id and ranks"""
    def __init__(self, seq_id, tax_id="0"):
        """
           ranks[6] = species 
           ranks[5] = genus
           ranks[4] = family
           ranks[3] = order
           ranks[2] = class
           ranks[1] = phylum 
           ranks[0] = kindom
        """  
        self.ranks = [""]*7
        self.curr_ranks_idx = 0
        self.seq_id = seq_id
        self.tax_id = tax_id
    
    
    def correct_rank_greengene(self):
        if self.ranks[0] == "k__":
            self.ranks[0:] = ["-"] * 7
            return 
        if self.ranks[1] == "p__":
            self.ranks[1:] = ["-"] * 6
            return 
        if self.ranks[2] == "c__":
            self.ranks[2:] = ["-"] * 5
            return 
        if self.ranks[3] == "o__":
            self.ranks[3:] = ["-"] * 4
            return
        if self.ranks[4] == "f__":
            self.ranks[4:] = ["-"] * 3
            return
        if self.ranks[5] == "g__":
            self.ranks[5:] = ["-"] * 2
            return
        if self.ranks[6] == "s__":
            self.ranks[6] = "-"
            return
    
    
    def __str__(self):
        allranks = ""
        for r in self.ranks:
            allranks =  allranks + ";" + r
        allranks = allranks[1:]
        allranks = self.seq_id + "  " + allranks
        return allranks

    """Not used"""
    def __assignRank(self, rank, name):
        name = name.replace(" ", "_")
        if rank == "species":
            self.ranks[0] = name
        elif rank == "genus":
            self.ranks[1] = name
        elif rank == "family":
            self.ranks[2] = name
        elif rank == "order":
            self.ranks[3] = name
        elif rank == "class":
            self.ranks[4] = name
        elif rank == "phylum":
            self.ranks[5] = name
        elif rank == "kindom":
            self.ranks[6] = name
    
    """This function is likely not needed"""
    def __checkRank(self, rank_id_list, id_rank_map, id_name_map):
        for i, rank in enumerate(self.ranks):
            if rank == "":
                print("checking ranks for " + self.seq_id)
                if i == 0:
                    for rid in rank_id_list:
                        rkname = id_rank_map[rid]
                        if rkname == "forma":
                            name = id_name_map[rid]
                            self.ranks[0] = name
                            break
                        elif rkname == "varietas":
                            name = id_name_map[rid]
                            self.ranks[0] = name
                            break
                        elif rkname == "subspecies":
                            name = id_name_map[rid]
                            self.ranks[0] = name
                            break
                        
                elif i == 1:
                    for rid in rank_id_list:
                        rkname = id_rank_map[rid]
                        if rkname == "subgenus":
                            name = id_name_map[rid]
                            self.ranks[1] = name
                            break
                        
                elif i == 2:
                    for rid in rank_id_list:
                        rkname = id_rank_map[rid]
                        if rkname == "subtribe":
                            name = id_name_map[rid]
                            self.ranks[2] = name
                            break
                        elif rkname == "tribe":
                            name = id_name_map[rid]
                            self.ranks[2] = name
                            break
                        elif rkname == "subfamily":
                            name = id_name_map[rid]
                            self.ranks[2] = name
                            break
                        elif rkname == "superfamily":
                            name = id_name_map[rid]
                            self.ranks[2] = name
                            break
                        
                elif i == 3:
                    for rid in rank_id_list:
                        rkname = id_rank_map[rid]
                        if rkname == "parvorder":
                            name = id_name_map[rid]
                            self.ranks[3] = name
                            break
                        elif rkname == "infraorder":
                            name = id_name_map[rid]
                            self.ranks[3] = name
                            break
                        elif rkname == "suborder":
                            name = id_name_map[rid]
                            self.ranks[3] = name
                            break
                        elif rkname == "superorder":
                            name = id_name_map[rid]
                            self.ranks[3] = name
                            break
                elif i == 4:
                    for rid in rank_id_list:
                        rkname = id_rank_map[rid]
                        if rkname == "infraclass":
                            name = id_name_map[rid]
                            self.ranks[4] = name
                            break
                        elif rkname == "subclass":
                            name = id_name_map[rid]
                            self.ranks[4] = name
                            break
                        elif rkname == "superclass":
                            name = id_name_map[rid]
                            self.ranks[4] = name
                            break
                elif i == 5:
                    for rid in rank_id_list:
                        rkname = id_rank_map[rid]
                        if rkname == "subphylum":
                            name = id_name_map[rid]
                            self.ranks[5] = name
                            break
                        elif rkname == "superphylum":
                            name = id_name_map[rid]
                            self.ranks[5] = name
                            break
    
    """This function is likely not needed"""
    def findMyRanksByTree(self, tax_tree):
        nodes = tax_tree.search_nodes(name=self.tax_id)
        if len(nodes)!=0:
            tax = nodes[0]
            self.__assignRank(tax.R, tax.N)
            while tax.up!=None:
                tax = tax.up
                self.__assignRank(tax.R, tax.N)
    
    """This function is likely not needed"""
    def findMyRanksByDB(self, id_pid_map, id_rank_map, id_name_map):
        rank_id_list = []
        currid = self.tax_id
        rank_id_list.append(currid)
        while currid!="1":
            pid = id_pid_map.get(currid, "1")
            currid = pid
            rank_id_list.append(currid)
        for rid in rank_id_list:
            rk = id_rank_map[rid]
            name = id_name_map[rid]
            self.__assignRank(rk, name) 
        self.__checkRank(rank_id_list, id_rank_map, id_name_map)



class seq_db:
    def __init__(self):
        self.tax_tree = ""
        self.name_taxid_map = {}
        self.seq_list = []
    
    def __str__(self):
        str1 = ""
        for seq in self.seq_list:
            str1 = str1 + seq.__str__() + "\n"
        return str1
    
    def init_db_greengene(self, sfin):
        with open(sfin) as fo:
            lines = fo.readlines()
        for line in lines:
            toks = line.strip().split("\t")
            seqid = toks[0].strip()
            ranks = toks[1].split(";")
            for i in range(len(ranks)):
                ranks[i]= ranks[i].strip()
            if len(ranks) == 7:
                seq = sequence(seq_id = seqid)
                seq.ranks = ranks
                seq.correct_rank_greengene()
                self.seq_list.append(seq)
            else:
                print("incomplete rank for seq: " + seqid)
                sys.exit()
    
    def get_origin_taxonomy_map(self):
        self.orign_taxonomy = {}
        for seq in self.seq_list:
            self.orign_taxonomy[seq.seq_id] = seq.ranks
        return self.orign_taxonomy
        
    
    """NCBI taxonomy support"""
    def init_db_from_file(self, sfin):
        fin = open(sfin)
        lines = fin.readlines()
        fin.close()
        for line in lines:
            line = line.strip()
            toks = line.split("|")
            if len(toks) == 7:
                seq = sequence(seq_id = toks[0])
                seq.ranks = toks[1:]
                seq.ranks.reverse()
                self.seq_list.append(seq)
            else:
                print("incomplete rank for seq" + toks[0])
                
    """NCBI taxonomy support"""
    def build_with_tax_tree(self, stax_tree, sname_taxid):
        self.tax_tree = Tree(stax_tree, format = 1)
        fname_taxid = open(sname_taxid)
        lines = fname_taxid.readlines()
        fname_taxid.close()
        for line in lines:
            line = line.strip()
            toks = line.split();
            self.name_taxid_map[toks[0]] = toks[1]
        
        for seqname in self.name_taxid_map.keys():
            taxid = self.name_taxid_map[seqname]
            seq = sequence(seq_id = seqname, tax_id = taxid)
            seq.findMyRanksByTree(self.tax_tree)
            self.seq_list.append(seq)
    
    """NCBI taxonomy support"""
    def build_with_db(self, s_id_name, s_id_pid_rank, sname_taxid):
        fname_taxid = open(sname_taxid)
        lines = fname_taxid.readlines()
        fname_taxid.close()
        for line in lines:
            line = line.strip()
            toks = line.split()
            self.name_taxid_map[toks[0]] = toks[1]
        
        id_pid_map = {}
        id_rank_map = {}
        id_name_map = {}
        f_id_pid_rank = open(s_id_pid_rank)
        lines = f_id_pid_rank.readlines()
        f_id_pid_rank.close()
        for line in lines:
            line = line.strip()
            toks = line.split("|")
            id_pid_map[toks[0]] = toks[1]
            id_rank_map[toks[0]] = toks[2]
        
        f_id_name = open(s_id_name)
        lines = f_id_name.readlines()
        f_id_name.close()
        for line in lines:
            line = line.strip()
            toks = line.split("|")
            id_name_map[toks[0]] = toks[1]
        
        for seqname in self.name_taxid_map.keys():
            taxid = self.name_taxid_map[seqname]
            seq = sequence(seq_id = seqname, tax_id = taxid)
            seq.findMyRanksByDB(id_pid_map, id_rank_map, id_name_map)
            self.seq_list.append(seq)
    
    def to_file(self, sfout):
        fout = open(sfout, "w")
        for seq in self.seq_list:
            fout.write(seq.__str__() + "\n")
        fout.close()
        
    def get_seq_by_name(self, seq_name):
        for seq in self.seq_list:
            if seq.seq_id == seq_name:
                return seq
        return None    
    
    def get_all_rank_names(self):
        """A lit of tubes (rank_name, rank_number)"""
        rk_name_set = set()
        for seq in self.seq_list:
            rks = seq.ranks
            for i, rk in enumerate(rks):
                if rk!="" and rk!="-":
                    rk_name_set.add((rk, i))
        return rk_name_set
    
    def get_non_species_rank_names(self):
        """A lit of tubes (rank_name, rank_number)"""
        rk_name_set = set()
        for seq in self.seq_list:
            rks = seq.ranks
            for i, rk in enumerate(rks):
                if rk!="" and rk!="-":
                    if i!=6:
                        rk_name_set.add((rk, i))
        return rk_name_set
    
    def rank_stas(self, seq_list):
        """find the most abundant ranks for sequence object list"""
        rank_stas = [["",0]] * 7
        numseq = float(len(seq_list))
        for i in range(7):
            rank_num_map = {}
            for seq in seq_list:
                rkname = seq.ranks[i]
                if rkname in rank_num_map:
                    rank_num_map[rkname] = rank_num_map[rkname] + 1
                else:
                    rank_num_map[rkname] = 1
            sorted_rank_num_map = sorted(rank_num_map.iteritems(), key=operator.itemgetter(1), reverse = True)
            maxrank = list(sorted_rank_num_map[0])
            maxrank[1] = float(maxrank[1]) / numseq
            rank_stas[i] = maxrank
        return rank_stas
    
    def rank_abundance(self, seq_list, rank_name, rank_num):
        num_seqs = float(len(seq_list))
        cnt = 0.0
        for seq in seq_list:
            ranks = seq.ranks
            if ranks[rank_num] == rank_name:
                cnt = cnt + 1.0
        return cnt/num_seqs



class phylogeny_annotator:
    """sphylogeny: birfurcating phlogenetic tree;  s_seq_db: greengene taxonomic flat file  """
    def __init__(self, sphylogeny, s_seq_db, t=0.99):
        self.tree_input = sphylogeny
        self.taxonomy_file = s_seq_db
        self.threshold = t
        self.root = Tree(sphylogeny, format=1)
        self.seqs = seq_db()
        self.seqs.init_db_greengene(s_seq_db)
        self.max_rank = 6 
        rks = self.seqs.get_all_rank_names()
        self.all_rank_names = []
        for rk in rks:
            self.all_rank_names.append(rk[0])
        
        self.nid_freq_map = {} # nid : [[r0name, f0],[r1name, f1],...,[r5name,f5]]
        self.nid_assigned_map = {} # nid : True or False, indicate if this node rank has been fully determined
        self.nid_ranks_map = {} #  nid : [r0name, r1name, ... ,r5name]
        self.nid_ranknum_map = {} #  nid : final_rank_num
    
    
    def __get_child_ranks(self, internal_node, rank_num):
        """input:internal node, rank_num; output: rankname frequency map """
        leaves = internal_node.get_leaves()
        rname_cnt_map = {}
        for leaf in leaves:
                seq = self.seqs.get_seq_by_name(leaf.name)
                rank_name = seq.ranks[rank_num]
                if rank_name!="-":
                    if rank_name in rname_cnt_map:
                        rname_cnt_map[rank_name] = rname_cnt_map[rank_name] + 1
                    else:
                        rname_cnt_map[rank_name] = 1 
        return rname_cnt_map
    
    
    def __sum_rank_num(self):
        """ideally every branch should be assainged to the lowest rank possbile, that is the higher ranknum the better"""
        s = 0
        for nid in self.nid_ranknum_map.keys():
            freq_table = self.nid_freq_map[nid]
            for ft in freq_table:
                s = s + ft[1]
        return s
    
    
    def __count_miss_labled(self):
        cnt = 0
        leave = self.root.get_leaves()
        for leaf in leave:
            oriranks = self.seqs.get_seq_by_name(leaf.name).ranks
            #oriranks.reverse()
            if self.nid_ranks_map[leaf.nid] != oriranks:
                print(leaf.name)
                print("Correct:" + str(self.nid_ranks_map[leaf.nid]))
                print("Misslab:" + str(oriranks))
                cnt = cnt + 1
        return cnt
        
    
    def __mono_index(self):
        """This method will calculate monophyly index by looking at the left and right hand side of the tree"""
        children = self.root.children
        if len(children) == 1:
            while len(children) == 1:
                children = children[0].children 
        if len(children) == 2:
            left = children[0]
            right =children[1]
            lset = self.__get_all_rank_names(left)
            rset = self.__get_all_rank_names(right)
            iset = lset & rset
            return len(iset)
        else:
            print("Fatel error: input tree not birfurcating")
            return None
    
    
    def __get_all_rank_names(self, node):
        nset = set()
        for n in node.traverse(strategy = "preorder"):
            ranks = self.nid_ranks_map[n.nid]
            for rk in ranks:
                if rk!="-":
                    nset.add(rk)
        return nset
    
    
    def __assign_all_descendent_node_rank(self, node, rank_num, rank_name, correct_errors = False):
        """Assign a rank name to all descendent nodes, and correct prior ranks errors if required"""
        descent_nodes = node.get_descendants()
        descent_nodes.append(node)
        find_error = False
        for nodei in descent_nodes:
            if nodei.is_leaf():
                ranks = self.nid_ranks_map[nodei.nid]
                ranks[rank_num] = rank_name
                self.nid_ranks_map[nodei.nid] = ranks
                
                if correct_errors:
                    seq = self.seqs.get_seq_by_name(nodei.name)
                    if seq.ranks[rank_num] != rank_name:
                        nodei.add_feature("is_correct", "No")
                        find_error = True
            else:
                ranks = self.nid_ranks_map[nodei.nid]
                ranks[rank_num] = rank_name
                self.nid_ranks_map[nodei.nid] = ranks
        
        if find_error and correct_errors:
            #recalculate all frequency vectors
            seq_util = seq_db()
            for nodei in node.traverse(strategy = "preorder"):
                leaves = nodei.get_leaves()
                seqs = []
                for leaf in leaves:
                    seqs.append(self.seqs.get_seq_by_name(leaf.name))
                    freq_table = seq_util.rank_stas(seqs)
                    self.nid_freq_map[nodei.nid] = freq_table
    
    
    def annotate_all_branches_td(self):
        self.nid_freq_map = {} # nid : [[r0name, f0],[r1name, f1],...,[r5name,f5]]
        self.nid_assigned_map = {} # nid : True or False, indicate if this node rank has been fully determined
        self.nid_ranks_map = {} #  nid : [r0name, r1name, ... ,r5name]
        self.nid_ranknum_map = {} #  nid : final_rank_num
        all_leaves = self.root.get_leaves()
        for leaf in all_leaves:
            leaf.add_feature("is_correct", "yes")
        
        """traversal the tree to: 1. add node id: nid 
                                  2. calculate the frequence table for each node/branch
                                  3. init ranks with "-"
        """
        seq_util = seq_db()
        i = 0
        for node in self.root.traverse(strategy = "preorder"):
            i = i + 1
            node.add_feature("nid", i)
            self.nid_assigned_map[i] = False
            ranks = ["-"] * 7
            self.nid_ranks_map[i] = ranks
            leaves = node.get_leaves()
            seqs = []
            for leaf in leaves:
                seqs.append(self.seqs.get_seq_by_name(leaf.name))
            freq_table = seq_util.rank_stas(seqs)
            self.nid_freq_map[i] = freq_table
        
        """traversal the tree preorder to assign ranks for each nodes"""
        for node in self.root.traverse(strategy = "preorder"):
            freq_table = self.nid_freq_map[node.nid]
            ranks = self.nid_ranks_map[node.nid]
            assigning_rank_idx = 0
            if node.is_root():
                self.nid_ranknum_map[node.nid] = -1
            else:
                next_rank_idx = self.nid_ranknum_map[node.up.nid] + 1
                flag = True
                while flag:
                    if next_rank_idx < 7:
                        rk_freq = freq_table[next_rank_idx]
                        if rk_freq[1] == 1.0:
                            self.__assign_all_descendent_node_rank(node, next_rank_idx, rk_freq[0])
                            next_rank_idx = next_rank_idx + 1
                        else:
                            childs = node.get_children()
                            lchild = childs[0]
                            rchild = childs[1]
                            lfreq_table = self.nid_freq_map[lchild.nid]
                            rfreq_table = self.nid_freq_map[rchild.nid]
                            lrk_freq = lfreq_table[next_rank_idx]
                            rrk_freq = rfreq_table[next_rank_idx]
                            #!!!!!core algorithm!!!!!!
                            if rk_freq[1] < lrk_freq[1] and rk_freq[1] < rrk_freq[1]:
                                flag = False
                            else: #TODO: checking all possibilties 
                                if lrk_freq[0] == rrk_freq[0] and rk_freq[1]>self.threshold:
                                    self.__assign_all_descendent_node_rank(node, next_rank_idx, rk_freq[0])
                                    next_rank_idx = next_rank_idx + 1
                                else:
                                    flag = False
                    else:
                        #assign taxonomy to species level
                        assigning_rank_idx = 6
                        rk_freq = freq_table[assigning_rank_idx]
                        flag = False
                        ranks[assigning_rank_idx] = rk_freq[0]
                self.nid_ranknum_map[node.nid] = next_rank_idx - 1
            self.nid_assigned_map[node.nid] = True
            
        return self.__mono_index()
    
    
    def show_tree_with_rank(self):
        allnodes = self.root.get_descendants()
        for node in allnodes:
            ranks = self.nid_ranks_map[node.nid]
            node.add_face(TextFace(str(ranks)), column=0, position = "branch-right")
            if node.is_leaf():
                seq = self.seqs.get_seq_by_name(node.name)
                rk = seq.ranks
                node.add_face(TextFace(str(rk)), column=0, position = "branch-right")
        self.root.show()
    
    
    def rooting_by_outgroup_names(self, outgroup_names):
        all_leaves = self.root.get_leaves()
        sog_names = set(outgroup_names)
        ca1 = self.root
        #Traversal all nodes to find the common ancestor of the input outgroup_names
        for node in self.root.traverse():
            currleaves = node.get_leaves()
            currlnames = []
            for lv in currleaves:
                currlnames.append(lv.name)
            scurrnames = set(currlnames)
            if scurrnames == sog_names:
                ca1 = node
            break
        #Check if the found ca is the root, if yes, find the complmentary names of the tree
        if ca1!=self.root:
            self.root.set_outgroup(ca1)
        else:
            restnodes = []
            for leaf in all_leaves:
                if leaf.name not in sog_names:
                    restnodes.append(leaf.name)
            srestnodes = set(restnodes)
            for node in self.root.traverse():
                currleaves = node.get_leaves()
                currlnames = []
                for lv in currleaves:
                    currlnames.append(lv.name)
                scurrnames = set(currlnames)
                if scurrnames == srestnodes:
                    self.root.set_outgroup(node)
                    break
    
    
    def annotate(self):
        """Auto root function"""
        #find all bipartations:
        list_bipar = []
        #all_leaves = self.root.get_leaves()
        for node in self.root.traverse("postorder"):
            if not node.is_root():
                leaves = node.get_leaves()
                leave_names = []
                for leaf in leaves:
                    leave_names.append(leaf.name) 
                list_bipar.append(leave_names)
        
        #find the root:
        min_mindex= 9999999999
        best_bipar = None
        for bipar in list_bipar:
            #Search the current tree to find the partitions:
            Node0 = self.root.search_nodes(name = bipar[0])[0]
            if len(bipar) == 1:
                self.root.set_outgroup(Node0)
            else:
                self.rooting_by_outgroup_names(bipar)
            mindex = self.annotate_all_branches_td()
            if mindex < min_mindex:
                min_mindex = mindex
                best_bipar = bipar
        
        if len(best_bipar) == 1:
            Node0 = self.root.search_nodes(name = best_bipar[0])[0]
            self.root.set_outgroup(Node0)
        else:
            self.rooting_by_outgroup_names(best_bipar)
        
        mindex = self.annotate_all_branches_td()
        print("Root with mono Index: " +  repr(mindex))
     
        
    def correct_leaf_ranks(self):
        """This method should be only called when the taxonomy is not congrunt with the phylogeny"""
        leaves = self.root.get_leaves()
        for leaf in leaves:
            if not leaf.is_root():
                father = leaf.up
                lranks = self.nid_ranks_map[leaf.nid]
                franks = self.nid_ranks_map[father.nid]
                lsp = lranks[6]
                for i, rk in enumerate(franks):
                    lranks[i] = rk
                lranks[6] = lsp 
                self.nid_ranks_map[leaf.nid] = lranks



class trainning:
    """This class parse the temp epa json file and generated a json file that has annotated taxonomy"""
    def __init__(self, temp_epa_json, ref_taxonomy, ref_sequences):
        self.jdata = json.load(open(temp_epa_json))
        self.jdata.pop("placements", None)
        self.jdata.pop("fields", None)
        self.jdata.pop("version", None)
        self.ref_taxonomy = ref_taxonomy
        self.ref_sequences = ref_sequences
        self.tree = self.jdata["tree"]
        self.tree = self.tree.replace("{", "[&&NHX:B=")
        self.tree = self.tree.replace("}", "]")
        self.pa = phylogeny_annotator(self.tree, self.ref_taxonomy)
        self.seqgroup = SeqGroup(sequences = self.ref_sequences)
        
    def trainning2json(self, fout, outgrop_names = ""):
        self.pa.annotate()
        bid_ranks_map = {}
        nid_ranks_map = self.pa.nid_ranks_map
        for node in self.pa.root.traverse("postorder"):
            ranks = nid_ranks_map.get(node.nid, ["-", "-", "-", "-", "-","-"])
            if hasattr(node, "B"):
                bid = node.B
                bid_ranks_map[bid] = ranks
        self.jdata["taxonomy"] = bid_ranks_map
        self.jdata["tree"] = self.tree
        self.jdata["raxmltree"] = Tree(self.tree, format=1).write(format=5)
        seqids = self.pa.root.get_leaf_names()
        seqs = []
        for entri in self.seqgroup.iter_entries():
            name = entri[0]
            if name in seqids:
                seqs.append(entri)
            
        self.jdata["sequences"] = seqs
        hmm = hmmer(refalign = self.ref_sequences)
        fprofile = hmm.build_hmm_profile()
        with open(fprofile) as fp:
            lines = fp.readlines()
        self.jdata["hmm_profile"] = lines
        os.remove(fprofile)
        self.jdata["origin_taxonomy"] = self.pa.seqs.get_origin_taxonomy_map()
        
        tp = tree_param(tree = self.tree, origin_taxonomy = self.jdata["origin_taxonomy"])
        self.jdata["rate"] = tp.get_speciation_rate()
        self.jdata["node_height"] = tp.get_nodesheight()
        self.jdata["version"] = "1.0"
        self.jdata["author"] = "Jiajie Zhang"
        
        with open(fout, "w") as fo:
            json.dump(self.jdata, fo, indent=4, sort_keys=True)



class leave_one_test:
    """This class will do a leave one test on every sequence in the alignment"""
    def __init__(self, salignment, stree, model):
        self.salign = salignment
        self.stre = stree
        self.tree = Tree(stree)
        self.alignment = salignment
        self.model = model
        
    
    def prune_one_tax(self, taxa_name):
        self.tree = Tree(self.stre)
        target_node = self.tree.search_nodes(name = taxa_name)[0]
        parent = target_node.up
        if parent.is_root():
            return -1
        else:
            sister_leaf_names = parent.get_leaf_names()
            sister_leaf_names.remove(taxa_name)
            
        target_node.delete(preserve_branch_length=True)
        epa_tree = self.tree.write(format = 5)
        EPA = epa()
        jp = EPA.run(reftree = epa_tree, alignment = self.salign, model = self.model)

        return self.check_placement(jp, sister_leaf_names)
    
    def check_placement(self, jp, location_leaf_names):
        newtree = Tree(jp.get_std_newick_tree())
        place_branch = str(jp.get_placement()[0]["p"][0][0])
        place_node = newtree.search_nodes(B=place_branch)[0]
        
        loca_leaf_nodes = []
        for n in location_leaf_names:
            nodes = newtree.search_nodes(name = n)
            loca_leaf_nodes.append(nodes[0])
        
        if len(loca_leaf_nodes) == 1:
            old_node = loca_leaf_nodes[0]
        else:
            old_node = loca_leaf_nodes[0].get_common_ancestor(loca_leaf_nodes[1:])
        
        if old_node.is_root():
            print("Wrong root in EPA tree")
            return -1
        else:
            dis = old_node.get_distance(target= place_node, topology_only=True)
            print("Placement distance:" + repr(dis))
            return dis
            
def testlo(fjson):
    jp = jsonparser(fjson)
    align = "tempaln"
    jp.get_alignment(align)
    
    t = jp.get_reftree()
    names = t.get_leaf_names()
    dissum = 0.0
    for name in names:
        lot = leave_one_test(align, jp.get_raxml_readable_tree())
        dis = lot.prune_one_tax(name)
        if dis >= 0:
            dissum = dissum + dis
    
    print dissum/float(len(names))

def testlo2(fjson):
    #jp = jsonparser(fjson)
    align = "/home/zhangje/GIT/EPA-classifier/example/testdata/clostridia_otu97/clostridia_otu97.fasta"
    t = Tree("/home/zhangje/GIT/EPA-classifier/example/testdata/clostridia_otu97/tax.bf.tre")
    EPA = epa()
    RA = raxml()
    fmodel = RA.get_model_parameters(bftree = "/home/zhangje/GIT/EPA-classifier/example/testdata/clostridia_otu97/tax.bf.tre", alignment = align)
    #RA.clean()
    
    names = t.get_leaf_names()
    dissum = 0.0
    for name in names:
        lot = leave_one_test(align, "/home/zhangje/GIT/EPA-classifier/example/testdata/clostridia_otu97/tax.bf.tre", fmodel)
        dis = lot.prune_one_tax(name)
        if dis >= 0:
            dissum = dissum + dis
    
    print dissum/float(len(names))


if __name__ == "__main__":
    if len(sys.argv) < 3: 
        print("usage: ./taxonomy_training.py <bifurcating_phylogenetic.tre> <greengene flat taxonomy file>")
        sys.exit()
    t = phylogeny_annotator(sphylogeny = sys.argv[1], s_seq_db = sys.argv[2])
    t.annotate()
    t.show_tree_with_rank()
    
