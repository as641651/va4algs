import os
import pandas as pd
import pickle

from am4pa.linnea import DataManagerLinnea
from variants_compare import VariantsCompare
from algorithm_ranking import MeasurementsVisualizer
from pm4py.objects.conversion.log import converter as log_converter

from .ranking_model import RankingModel


class RankingDataLinnea:
    def __init__(self, dml:DataManagerLinnea, rm:RankingModel, thread_str):
        self.dml = dml
        self.rm = rm
        self.thread_str = thread_str
        
        self.data_vcs_flops = {}
        #self.data_vcs_nflops = {}
        self.data_kernels = None
        self.data_relations = None
        self.data_best_kseq = None
        self.data_worst_kseq = None
        self.data_ext = None
        self.data_ranks = {}
        self.data_h0 = {}
        self.obj_path = os.path.join(dml.lc.local_dir,'ranking-data','rdl_{}.pkl'.format(thread_str))
        
        
    def rank3way(self):
        
        data_nodes = []
        data_edges = []
        data_ext = []
        dbest_a = []
        dworst_a = []
        
        for op_str, ml in self.dml.mls[self.thread_str].items():
            
            #collect data
            ml.case_durations_manager.clear_case_durations()
            for i in self.dml.measurements_data[self.thread_str][op_str]:
                ml.collect_measurements(i)
                
            ranks,cutoffs,h0_ = self.rm.get_ranks(ml.get_alg_measurements())
            best_algs = ranks[ranks.iloc[:,1]<=cutoffs[0]]['case:concept:name'].tolist()
            worst_algs = ranks[ranks.iloc[:,1]>cutoffs[1]]['case:concept:name'].tolist()

            dc = ml.data_collector
            et = ml.filter_table(dc.get_meta_table())
            et['concept:name'] = et['concept:name'].apply(lambda row: self._clean_concept_eq(row))
            et['concept:name'] = et['concept:name'].apply(lambda row: self._clean_concept_remove_LAPACK(row))

            for alg in best_algs:
                dbest_a.append((et[et['case:concept:name']==alg]['concept:name'].apply(lambda x: x.split('_')[0]).tolist()))

            for alg in worst_algs:
                dworst_a.append((et[et['case:concept:name']==alg]['concept:name'].apply(lambda x: x.split('_')[0]).tolist()))

            xes_log = log_converter.apply(et)

            activity_key = 'concept:name'
            vc = VariantsCompare(xes_log,best_algs,worst_algs,activity_key=activity_key)
            dn, de = vc.get_diff_data()
            
            dn['operands'] = op_str
            de['operands'] = op_str

            ct = ml.filter_table(dc.get_case_table())
            min_flop = ct['case:flops'].min()
            ct['case:rel-flops'] = ct.apply(lambda row: (row['case:flops'] - min_flop) / min_flop, axis=1)
            et = et.merge(ct, on='case:concept:name')
            et['kernel'] = et.apply(lambda x: x['concept:name'].split('_')[0], axis=1)
            ext = et[['kernel', 'concept:flops', 'case:rel-flops']]
            ext = ext.drop_duplicates().reset_index(drop=True)
            
            #ext = et.drop_duplicates(subset=['concept:name'])[['concept:name', 'concept:flops']]
            data_nodes.append(dn)
            data_edges.append(de)
            data_ext.append(ext)
            
            self.data_vcs_flops[op_str] = vc
            self.data_ranks[op_str] = ranks
            self.data_h0[op_str] = h0_

        self.data_kernels = pd.concat(data_nodes).reset_index(drop=True)
        self.data_relations = pd.concat(data_edges).reset_index(drop=True)
        
        def get_flops(str_):
            if not '@@' in str_:
                return float(str_.split('_')[1])
            return 0
        self.data_kernels['flops'] = self.data_kernels.apply(lambda x: get_flops(x['node']), axis=1)
        self.data_kernels['kernel'] = self.data_kernels.apply(lambda x: x['node'].split('_')[0], axis=1)
        self.data_relations['flopsA'] = self.data_relations.apply(lambda x: get_flops(x['nodeA']), axis=1)
        self.data_relations['kernelA'] = self.data_relations.apply(lambda x: x['nodeA'].split('_')[0], axis=1)
        self.data_relations['flopsB'] = self.data_relations.apply(lambda x: get_flops(x['nodeB']), axis=1)
        self.data_relations['kernelB'] = self.data_relations.apply(lambda x: x['nodeB'].split('_')[0], axis=1)
        
        self.data_best_kseq =  dbest_a
        self.data_worst_kseq = dworst_a
        self.data_ext = pd.concat(data_ext).reset_index(drop=True)
        
        
    def _clean_concept_eq(self, name):
        splits = name.split('=')
        if len(splits) > 1:
            return splits[-1].strip()
        return splits[0].strip()
    
    def _clean_concept_remove_cost(self,name):
        splits = name.split('_')
        if len(splits) > 1:
            return splits[0].strip()
        return splits[0].strip()

    def _clean_concept_remove_LAPACK(self,name):
        splits = name.split('LAPACK.')
        if len(splits) > 1:
            return splits[-1].strip()
        return splits[0].strip()
    
    def visualize_box_plots(self, op_str, scale=0.8, tick_size=16):
        ml = self.dml.mls[self.thread_str][op_str]
        mv = MeasurementsVisualizer(ml.get_alg_measurements(), self.data_h0[op_str])
        fig = mv.show_measurements_boxplots(scale=scale,tick_size=tick_size)
        return fig
    
    
    def __getstate__(self):
        state = {
            'kernels':self.data_kernels,
            'relations':self.data_relations,
            'ranks':self.data_ranks,
            'h0':self.data_h0,
            'vcs_f':self.data_vcs_flops,
            'best_k':self.data_best_kseq,
            'worst_k':self.data_worst_kseq,
            'ext':self.data_ext
        }
        return state
    
    def __setstate__(self,state):
        self.data_kernels = state['kernels']
        self.data_relations = state['relations']
        self.data_ranks = state['ranks']
        self.data_h0 = state['h0']
        self.data_vcs_flops = state['vcs_f']
        self.data_best_kseq = state['best_k']
        self.data_worst_kseq = state['worst_k']
        self.data_ext = state['ext']
        
    def save(self):
        if not os.path.exists(os.path.dirname(self.obj_path)):
            os.makedirs(os.path.dirname(self.obj_path))
            
        with open(self.obj_path,"wb") as f:
            pickle.dump(self,f)
            
    def load(self):
        if not os.path.exists(self.obj_path):
            return -1
        
        with open(self.obj_path,"rb") as f:
            self.__setstate__(pickle.load(f).__getstate__())
            
        for op_str, ml in self.dml.mls[self.thread_str].items():
            
            #collect data
            ml.case_durations_manager.clear_case_durations()
            for i in self.dml.measurements_data[self.thread_str][op_str]:
                ml.collect_measurements(i)
        
    