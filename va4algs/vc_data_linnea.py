from variants_compare import VariantsCompare
from pm4py.objects.conversion.log import converter as log_converter
from .ranking_data_linnea import RankingDataLinnea

class VCDataLinnea:
    def __init__(self,rdl:RankingDataLinnea):
        self.rdl = rdl
        self.dml = self.rdl.dml
        self.thread_str = self.rdl.thread_str
        
        ret = self.rdl.load()
        assert ret != -1, "First rank data"
        
        self.da = self.rdl.data_anomalies
        
        ## Regular VC: Rank 0 vs none
        self.vc = {}
        
        ## VC type 1: compare only Min flops
        self.vc1 = {}
        
        ## VC type 2: compare rank 0 algs vs min flops algs
        self.vc2 = {}
        
        ## VC type 3: compare algs 0 vs rest upto relflops cutoff
        self.vc3 = {}
        
        
          
    def get_vc(self,op_str):
        
        try:
            return 1, self.vc[op_str]
        except KeyError:
            pass
        
        ranks = self.rdl.data_ranks[op_str]
        rank0_algs = set(ranks[ranks.iloc[:,1]==0]['case:concept:name'])
        worst_algs = set(ranks[ranks.iloc[:,1]!=0]['case:concept:name'])
        
        if len(worst_algs) == 0:
            print("All algs are ranked 0")
            return -1, None
        
        vc = self._get_vc(op_str,rank0_algs,worst_algs)
        self.vc[op_str] = vc
        
        return 1,vc
        
    
    def get_vc1(self,op_str):
        
        try:
            return 1, self.vc1[op_str]
        except KeyError:
            pass
        
        rf_cutoff = self.da[self.da['op_str']==op_str]['rel-flops-cutoff'].values[0]
        if rf_cutoff > 0.0:
            print("No min flop alg in best variants")
            return -1,None
        
        adj_risk = self.da[self.da['op_str']==op_str]['adj_risk'].values[0]
        if adj_risk == 0:
            print("All min flop algs in best variants")
            return -1,None
        
        ranks = self.rdl.data_ranks[op_str]
        
        rel0_algs = set(ranks[ranks['case:rel-flops']==0]['case:concept:name'])
        rank0_algs = set(ranks[ranks.iloc[:,1]==0]['case:concept:name'])
        
        best_algs = rank0_algs.intersection(rel0_algs)
        worst_algs = rel0_algs - rank0_algs
        
        vc = self._get_vc(op_str,best_algs,worst_algs)
        
        self.vc1[op_str] = vc
        
        return 1,vc
    
    
    def get_vc2(self,op_str):
        
        try:
            return 1, self.vc2[op_str]
        except KeyError:
            pass
        
        
        adj_risk = self.da[self.da['op_str']==op_str]['adj_risk'].values[0]
        if adj_risk == 0:
            print("All min flop algs in best variants")
            return -1,None
        
        
        ranks = self.rdl.data_ranks[op_str]
        
        rel0_algs = set(ranks[ranks['case:rel-flops']==0]['case:concept:name'])
        rank0_algs = set(ranks[ranks.iloc[:,1]==0]['case:concept:name'])
        
        worst_algs = rel0_algs - rank0_algs
        
        #print(len(worst_algs),len(rel0_algs))
        
        vc = self._get_vc(op_str,rank0_algs,worst_algs)
        
        self.vc2[op_str] = vc
        
        return 1,vc
        
                
    
    def get_vc3(self,op_str):
                
        try:
            return 1, self.vc3[op_str]
        except KeyError:
            pass
        
        
        ranks = self.rdl.data_ranks[op_str]
        rank0_algs = set(ranks[ranks.iloc[:,1]==0]['case:concept:name']) 
        
        rel_flops_cutoff = self.da[self.da['op_str']==op_str]['rel-flops-cutoff'].values[0]   
        worst_algs = set(ranks[ranks['case:rel-flops']<=rel_flops_cutoff]['case:concept:name']) - rank0_algs
        
        if len(worst_algs) == 0:
            print("All algs upto {} FLOPs are in the best algs".format(rel_flops_cutoff))
            return -1, None
        
        #print(len(worst_algs),len(rank0_algs))
        vc = self._get_vc(op_str,rank0_algs,worst_algs)
        self.vc3[op_str] = vc
        
        return 1,vc
        
        
    
    
    def _get_vc(self,op_str,best_a,worst_a):
        ml = self.dml.mls[self.thread_str][op_str]
        
        et = self._get_et(ml)
        
        xes_log = log_converter.apply(et)
        
        activity_key = 'concept:name'
        vc = VariantsCompare(xes_log,list(best_a),list(worst_a),activity_key=activity_key)
        
        return vc
    
    def _get_et(self,ml):
        dc = ml.data_collector
        et = ml.filter_table(dc.get_meta_table())
        et['concept:name'] = et['concept:name'].apply(lambda row: self._clean_concept_eq(row))
        et['concept:name'] = et['concept:name'].apply(lambda row: self._clean_concept_remove_LAPACK(row))
        return et
        
    
    
    def _clean_concept_eq(self, name):
        splits = name.split('=')
        if len(splits) > 1:
            return splits[-1].strip()
        return splits[0].strip()


    def _clean_concept_remove_LAPACK(self,name):
        splits = name.split('LAPACK.')
        if len(splits) > 1:
            return splits[-1].strip()
        return splits[0].strip()
        