from algorithm_ranking import RankVariantsDFGTr

class RankingModel:
    def __init__(self):
        pass
    
    def train(self):
        pass
    
    def methodology2(self,alg_measurements,alg_list):
        rv = RankVariantsDFGTr(alg_measurements, alg_list)
        ranks =  rv.rank_variants_reliable()[0].iloc[:,:2]
        h0_ = rv.graph.get_separable_arrangement()
        return ranks,h0_
    
    def methodology3(self,alg_measurements,alg_list):
        pass
        
    def get_ranks(self,alg_measurements, alg_list=None):
        if not alg_list:
            alg_list = list(alg_measurements.keys())
        ## TODO ML model: decide methodolog2 or methodology3 and cutoff points
        ranks,h0_ = self.methodology2(alg_measurements,alg_list)    
        cutoffs = [0,0]
        return ranks, cutoffs, h0_