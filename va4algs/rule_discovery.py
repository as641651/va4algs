from enum import Enum
import pandas as pd
from va4algs import RankingDataLinnea
from va4algs import VCDataLinnea

class AnalysisType(Enum):
    TYPE1 = 1
    TYPE2 = 2
    TYPE3 = 3
    TYPE4 = 4
    TYPE5 = 5

class RulesDiscovery:
    def __init__(self,rdl:RankingDataLinnea):
        self.rdl = rdl
        self.dml = rdl.dml
        self.thred_str = rdl.thread_str
        
        self.vc = VCDataLinnea(self.rdl)
        self.da = self.rdl.data_anomalies
        
        self.data_kernels = {}
        self.data_relations = {}
        
    def get_analysis_data(self, a_type:AnalysisType):
        try:       
            return self.data_kernels[a_type], self.data_relations[a_type]
        except KeyError:
            pass
            
        
        nodes = []
        edges = []
        
        if a_type == AnalysisType.TYPE1:
            ## only min flops
            op_strs = self.da[self.da['adj_risk']>0.0]['op_str'].tolist()
            vc_f = self.vc.get_vc1
        elif a_type == AnalysisType.TYPE2:
            ## rank 0 vs min flops for rel flops > 0
            op_strs = self.da[self.da['rel-flops-cutoff']>0.0]['op_str'].tolist()
            vc_f = self.vc.get_vc2
        elif a_type == AnalysisType.TYPE3: 
            ## only rel-flops cut off > 0
            op_strs = self.da[self.da['rel-flops-cutoff']>0.0]['op_str'].tolist()
            vc_f = self.vc.get_vc3
        elif a_type == AnalysisType.TYPE4:  
            ## vc3 with adj_risk > 0
            op_strs = self.da[self.da['adj_risk']>0.0]['op_str'].tolist()
            vc_f = self.vc.get_vc3
        elif a_type == AnalysisType.TYPE5:
            ## All ops best vs rest
            op_strs = self.da['op_str'].tolist()
            vc_f = self.vc.get_vc     
        else:
            print("Choose a relavant analysis type: 1, 2 or 3")
            return    
            
        
        for op_str in op_strs:
            
            ret, vc = vc_f(op_str)
            if ret == 1: 
                dn, de = vc.get_diff_data()
                dn['operands'] = op_str
                de['operands'] = op_str

                nodes.append(dn)
                edges.append(de)
            
        if nodes:
            self.data_kernels[a_type] = pd.concat(nodes).reset_index(drop=True)
            self.data_relations[a_type] = pd.concat(edges).reset_index(drop=True)
            self.data_kernels[a_type], self.data_relations[a_type] = self._prepare_kernel_relations_data(self.data_kernels[a_type],
                                                                                            self.data_relations[a_type])
        
            return self.data_kernels[a_type], self.data_relations[a_type]
        else:
            return None,None

    
    def discover_kernel_rules(self,a_type:AnalysisType):
        df_k, _ = self.get_analysis_data(a_type)
        
        sa = df_k.groupby(["kernel", "class"])['class'].count().unstack(fill_value=0)
        sa = self._add_missing_class_data(sa)
        sa = self._compute_class_statistics(sa)
        s2 = sa[(sa['good (%)'] == 1.0) | (sa['bad (%)'] == 1.0)].reset_index()
        
        return s2
    
    def discover_relations_rules(self, a_type):
        _, df_r = self.get_analysis_data(a_type)
        
        sa = df_r.groupby(["kernelA", "kernelB", "class"])['class'].count().unstack(fill_value=0)
        sa = self._add_missing_class_data(sa)
        sa = self._compute_class_statistics(sa)
        
        s2 = sa[(sa['good (%)'] == 1.0) | (sa['bad (%)'] == 1.0)].reset_index()
#         if filter_kernel_rules:
#             bl = self.rules['kernel'][1] + self.rules['kernel'][-1]
#             return s2[(~s2['kernelA'].isin(bl)) & (~s2['kernelB'].isin(bl))].reset_index(drop=True)
        return s2
    

    def _prepare_kernel_relations_data(self,data_kernels,data_relations):
        def get_flops(str_):
            if not '@@' in str_:
                return float(str_.split('_')[1])
            return 0
        data_kernels['flops'] = data_kernels.apply(lambda x: get_flops(x['node']), axis=1)
        data_kernels['kernel'] = data_kernels.apply(lambda x: x['node'].split('_')[0], axis=1)
        data_relations['flopsA'] = data_relations.apply(lambda x: get_flops(x['nodeA']), axis=1)
        data_relations['kernelA'] = data_relations.apply(lambda x: x['nodeA'].split('_')[0], axis=1)
        data_relations['flopsB'] = data_relations.apply(lambda x: get_flops(x['nodeB']), axis=1)
        data_relations['kernelB'] = data_relations.apply(lambda x: x['nodeB'].split('_')[0], axis=1)
        
        return data_kernels, data_relations
    
    def _add_missing_class_data(self, df):
        if not 1 in df.columns:
            df[1] = 0
        if not -1 in df.columns:
            df[-1] = 0
        return df
    
    def _compute_class_statistics(self, df):
        df['total'] = df.apply(lambda x: x[-1]+x[0]+x[1], axis=1)
        df['good (%)'] = df.apply(lambda x: x[1]/float(x['total']), axis=1).round(2)
        df['bad (%)'] = df.apply(lambda x: x[-1]/float(x['total']), axis=1).round(2)
        df['selection score'] = df.apply(lambda x: x[1]/float(x[1]+ x[-1]), axis=1).fillna(-1).round(2)
        return df.reset_index()
        
            