import pandas as pd
from am4pa.linnea import LinneaConfig, MeasurementsLinnea
from algorithm_ranking import RankVariantsDFGTr, MeasurementsVisualizer
from pm4py.objects.conversion.log import converter as log_converter
from variants_compare import VariantsCompare
from typing import List

from .operands_manager import OperandsManager

class DataExtractorLinnea:
    def __init__(self, linnea_config:LinneaConfig):
        self.lc = linnea_config
        self.om = OperandsManager() 
        self.mls = {}
        self.backend = linnea_config.backend
        self.run_ids = set()
        self.Ranker = RankVariantsDFGTr
        self.bRankReliable = True
        self.q_max = 75
        self.q_min = 25
        
        self.data_vcs = {}
        self.data_kernels = None
        self.data_relations = None
        self.data_best_kseq = None
        self.data_worst_kseq = None
        self.data_ext = None
        self.data_ranks = {}
        
    def add_operands(self, ops:List, bSync=True):
        self.om.add_operands(ops)
        for op in ops:
            self.mls[self.om.op_str(op)] = MeasurementsLinnea(self.lc, op)
            if bSync:
                self._sync_data(self.mls[self.om.op_str(op)])
    
    def _sync_data(self, ml):
        if self.backend:
            if self.lc.bm.check_if_dir_exists(ml.runner.operands_dir):
                self.om.set_generated(ml.op_sizes)
                ret = ml.data_collector.get_runtimes_table()
                if isinstance(ret, pd.DataFrame):
                    ml.measured_once = True
        else:
            #TODO
            pass
        
    def generate_variants(self):
        for k,v in self.mls.items():
            if not self.om.is_generated(k.split("_")):
                print("Generating Variants for {}".format(k))
                v.generate_variants()
                self.om.set_generated(k.split("_"))
                
    def gather_all_variants(self):
        for k,v in self.mls.items():
            v.gather_all_variants()
            
    def filter_on_flops(self, rel_flops=1.2):
        for k,v in self.mls.items():
            v.filter_on_flops(rel_flops=rel_flops)
            
    def filter_on_rel_duration(self, rel_d=1.2):
        for k,v in self.mls.items():
            if not v.measured_once:
                v.measure_once()
            v.filter_on_flops_rel_duration(rel_duration=rel_d)
            
            
    def check_measured(self, run_id):
        for k,v in self.mls.items():
            df = v.data_collector.get_runtimes_competing_table(run_id)
            if isinstance(df, pd.DataFrame):
                self.om.set_measured(k.split("_"))    
                self.run_ids.add(run_id)
                
    def measure_variants(self, reps=5, run_id=0):
        if not run_id in self.run_ids:
            self.reset_measured()
            self.run_ids.add(run_id)
        
        for k,v in self.mls.items():
            if not self.om.is_measured(k.split("_")): 
                print("Measuring Variants for {}".format(k))
                v.measure(reps, run_id)
                self.om.set_measured(k.split("_")) 
                
    def reset_measured(self):
        self.om.df['measured'] = False
        
        
    def _clean_concept(self, name):
        splits = name.split('=')
        if len(splits) > 1:
            return splits[-1].strip()
        return splits[0].strip()

    def _rank_variants(self, alg_measurements, alg_list):
        rv = self.Ranker(alg_measurements, alg_list)
        return rv.rank_variants(q_max = self.q_max, q_min = self.q_min)

    def _rank_variants_reliable(self, alg_measurements, alg_list):
        rv = self.Ranker(alg_measurements, alg_list)
        return rv.rank_variants_reliable()[0].iloc[:,:2]
    
    
    def visualize_box_plots(self, op_str, scale=0.8, tick_size=16):
        ml = self.mls[op_str]
        mv = MeasurementsVisualizer(ml.get_alg_measurements(), ml.h0)
        fig = mv.show_measurements_boxplots(scale=scale,tick_size=tick_size)
        return fig
    
    def prepare_data_for_analysis(self, cutoff=0, clear=False):
        #vc = {}
        data_nodes = []
        data_edges = []
        data_ext = []
        dbest_a = []
        dworst_a = []
        
        for op_str, ml in self.mls.items():
            
            #collect data
            if clear:
                ml.data_collector.delete_local_data()
            ml.case_durations_manager.clear_case_durations()
            for i in self.run_ids:
                ml.collect_measurements(i)
                
            #dfs.append(rank_variants_dfg_reliable(m.get_alg_measurements(), m.h0))
            if self.bRankReliable:
                ranks = self._rank_variants_reliable(ml.get_alg_measurements(), ml.h0)
            else:
                ranks = self._rank_variants(ml.get_alg_measurements(), ml.h0)
            best_algs = ranks[ranks.iloc[:,1]<=cutoff]['case:concept:name'].tolist()
            worst_algs = ranks[ranks.iloc[:,1]>cutoff]['case:concept:name'].tolist()

            dc = ml.data_collector
            et = ml.filter_table(dc.get_meta_table())
            et['concept:name'] = et['concept:name'].apply(lambda row: self._clean_concept(row))

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

            
            data_nodes.append(dn)
            data_edges.append(de)
            data_ext.append(ext)
            
            self.data_vcs[op_str] = vc
            self.data_ranks[op_str] = ranks

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
        
