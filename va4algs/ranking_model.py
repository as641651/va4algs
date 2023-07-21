from algorithm_ranking import RankVariantsDFGTr
from algorithm_ranking import RankVariantsSort2

class RankingModel:
    def __init__(self,name,method=2):
        self.name = name
        self.method = method
    
    def train(self):
        pass
    
    def methodology2(self,alg_measurements,alg_list):
        rv = RankVariantsDFGTr(alg_measurements, alg_list)
        ranks =  rv.rank_variants_reliable()[0].iloc[:,:2]
        h0_ = rv.graph.get_separable_arrangement()
        return ranks,h0_
    
    def methodology3(self,alg_measurements,alg_list):
        rv_ = RankVariantsDFGTr(alg_measurements, alg_list)
        #ranks_ =  rv_.rank_variants_reliable()[0].iloc[:,:2]
        ranks_ = rv_.rank_variants()
        h0_ = rv_.graph.get_separable_arrangement()
        rv = RankVariantsSort2(alg_measurements,h0_)
        ranks = rv.rank_variants()
        return ranks, h0_
        
    def get_ranks(self,alg_measurements, alg_list=None):
        if not alg_list:
            alg_list = list(alg_measurements.keys())
        ## TODO ML model: decide methodolog2 or methodology3 and cutoff points
        if self.method == 2:
            ranks,h0_ = self.methodology2(alg_measurements,alg_list)  
        elif self.method==3:  
            ranks,h0_ = self.methodology3(alg_measurements,alg_list)
        else:
            print("Ranking method not implemented")
            exit(-1)  
        cutoffs = [0,0]
        return ranks, cutoffs, h0_