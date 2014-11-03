#! /usr/bin/env python
from epac.taxonomy_util import Taxonomy
    
class TaxClassifyHelper:
    def __init__(self, cfg, bid_taxonomy_map):
        self.cfg = cfg
        self.bid_taxonomy_map = bid_taxonomy_map

    def classify_seq(self, edges, method = "1", minlw = 0.):
        if len(edges) > 0:
            if method == "1":
                ranks, lws = self.assign_taxonomy_maxsum(edges, minlw)
            else:
                ranks, lws = self.assign_taxonomy_maxlh(edges)
            return ranks, lws
        else:
            return None, None      
     
    def assign_taxonomy_maxlh(self, edges):
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

    def assign_taxonomy_maxsum(self, edges, minlw):
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
                rank_id = Taxonomy.get_rank_uid(ranks, i)
                if rank != Taxonomy.EMPTY_RANK:
                    rw_total[rank_id] = rw_total.get(rank_id, 0) + lweight
                    lowest_rank = rank_id
                    if not rank_id in rb:
                        rb[rank_id] = br_id
                else:
                    break

            if lowest_rank:
                rw_own[lowest_rank] = rw_own.get(lowest_rank, 0) + lweight
                rb[lowest_rank] = br_id
            elif self.cfg.verbose:
                print "WARNING: no annotation for branch ", br_id
            
        # if all branches have empty ranks only, just return this placement
        if len(rw_total) == 0:
            return ranks, [1.] * len(ranks)
        
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
                rank_id = Taxonomy.get_rank_uid(a_ranks, i)
                a_conf[i] = rw_total[rank_id]

        return a_ranks, a_conf
    
