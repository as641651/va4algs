import pandas as pd
from typing import List
from sklearn.tree import DecisionTreeClassifier
from .utils import get_tree_gviz, gini_impurity, plot_clf_regions

class DataAnalyzer3Way:
    def __init__(self,df_kernels:pd.DataFrame, 
                 df_relations:pd.DataFrame,
                 df_extra = None):
        
        self._df_k = df_kernels
        self._df_r = df_relations
        self._df_ext = df_extra
        self.kernel_stats = None
        self.relations_stats = None
        
        self.classifiers = {}
        
        self.update_data()
        
        self.rules = {}
        self.rules['kernel'] = {}
        self.rules['kernel'][1] = []
        self.rules['kernel'][-1] = []
        
        self.rules['relation'] = {}
        self.rules['relation'][1] = []
        self.rules['relation'][-1] = []
        
        
    def update_data(self, df_kernels:pd.DataFrame=None, 
                    df_relations:pd.DataFrame=None, 
                    df_extra:pd.DataFrame=None):
        
        if df_kernels:
            self._df_k = df_kernels
        if df_relations:
            self._df_r = df_relations
        if df_extra:
            self._df_ext = df_extra
            
        self._prepare_kernel_stats()
        self._prepare_relations_stats()
        
    def discover_kernel_rules(self):
        # Explainable if rel flops > 0
        # Surprising otherwise
        s = self.kernel_stats
        s2 = s[(s['good (%)'] == 1.0) | (s['bad (%)'] == 1.0)].reset_index()
        self.rules['kernel'][1] =  s2[s2['good (%)'] == 1.0]['kernel'].to_list()
        self.rules['kernel'][-1] =  s2[s2['bad (%)'] == 1.0]['kernel'].to_list()
        return s2
    
    def discover_relations_rules(self, filter_kernel_rules = True):
        s = self.relations_stats
        s2 = s[(s['good (%)'] == 1.0) | (s['bad (%)'] == 1.0)].reset_index()
        if filter_kernel_rules:
            bl = self.rules['kernel'][1] + self.rules['kernel'][-1]
            return s2[(~s2['kernelA'].isin(bl)) & (~s2['kernelB'].isin(bl))].reset_index(drop=True)
        return s2
    
    def discover_kernel_anomalies(self):
        s = self.kernel_stats
        return s[(s['case:rel-flops'] > 0.0) & (s['selection score'] == 1.0)]
    
    def get_kernel_stats(self, kernel):
        return self.kernel_stats[self.kernel_stats['kernel']==kernel]
    
    def get_relation_stats(self, kernelA, kernelB):
        s = self.relations_stats
        return s[(s['kernelA']==kernelA) & (s['kernelB']==kernelB)]
        
    
    def discover_impurities(self):
        s = self._df_r
        s4 = s.groupby(["kernelA", "flopsA", "kernelB", "flopsB"])['class'].apply(list).reset_index()
        s4['impurity'] = s4.apply(lambda x: gini_impurity(x['class']), axis=1)
        return s4[s4['impurity'] != 0.0]
    
    def filter_relations_not_in_rules(self):
        s = self.relations_stats
        s = s[(s['good (%)'] != 1.0) & (s['bad (%)'] != 1.0)].reset_index().sort_values(by=['selection score'])
        return s[s['selection score'] != -1.0].reset_index()
    
    
    def filter_relations_between_selection_scores(self, upper=1.0, lower=0.0):
        s = self.relations_stats
        
        return s[(s['selection score'] < upper) &  (s['selection score'] > lower) ]\
        .reset_index().sort_values(by=['selection score'])
        
    def filter_relations_on_selection_scores(self, val1, val2=None):
        s = self.relations_stats
        
        if val2 == None:
            val2 = val1
        
        return s[(s['selection score'] == val1) |  (s['selection score'] == val2) ]\
        .reset_index().sort_values(by=['selection score'])
    
    def train_dtree_classifiers(self, kernelsA:List, kernelsB:List):
        
        for kernelA, kernelB in zip(kernelsA, kernelsB):
            k_str = self.get_cls_key(kernelA, kernelB)
            X = self.get_relation_data(kernelA, kernelB)
            clf = DecisionTreeClassifier(random_state=1234)
            model = clf.fit(X.iloc[:,:-1], X.iloc[:,-1])
            self.classifiers[k_str] = clf
            
    
    def get_cls_key(self, kernelA, kernelB):
        return "{}/{}".format(kernelA, kernelB)
    
    
    def get_relation_data(self, kernelA, kernelB):
        df = self._df_r
        x = df[(df['kernelA'] == kernelA) & (df['kernelB'] == kernelB)]
        x['flopsA'] = x['flopsA']*10**-6
        x['flopsB'] = x['flopsB']*10**-6
        return x[['flopsA', 'flopsB', 'class']].rename(columns={
            'flopsA':kernelA, 'flopsB':kernelB
        })
    
    def get_regions_plot(self, kernelA, kernelB):
        X = self.get_relation_data(kernelA, kernelB)
        fig = plot_clf_regions(X.values[:,:-1], X.values[:,-1], 
                               X.columns[0], X.columns[1],
                              self.classifiers[self.get_cls_key(kernelA, kernelB)])
        return fig
    
    def get_dtree(self, kernelA, kernelB):
        g = get_tree_gviz(self.classifiers[self.get_cls_key(kernelA, kernelB)],
                          kernelA, kernelB)
        return g
    
    def find_operands_on_cl(self, cl:int, kernelA:str, kernelB:str=None):
        df = self._df_r
        
        if kernelB:
            condition = (df['kernelA'] == kernelA) & \
                    (df['kernelB'] == kernelB) & \
                    (df['class'] == cl)
        else:
            condition = (df['kernelA'] == kernelA) | \
                    (df['kernelB'] == kernelA) & \
                    (df['class'] == cl)
    
        return list(set(df[condition]['operands'].to_list()))
    
    def find_operands_on_flops(self, kernelA:str, kernelB:str,
                           ineqA:str, flopsA:float,
                           ineqB:str, flopsB:float):
        df = self._df_r
        condition = None
        
        if ineqA == "<" and ineqB == "<":
            condition = (df['kernelA'] == kernelA) & \
                    (df['flopsA'] < flopsA*10**6) & \
                    (df['kernelB'] == kernelB) & \
                    (df['flopsB'] < flopsB*10**6)
        elif ineqA == ">" and ineqB == ">":
            condition = (df['kernelA'] == kernelA) & \
                    (df['flopsA'] > flopsA*10**6) & \
                    (df['kernelB'] == kernelB) & \
                    (df['flopsB'] > flopsB*10**6)
        elif ineqA == "<" and ineqB == ">":
            condition = (df['kernelA'] == kernelA) & \
                    (df['flopsA'] < flopsA*10**6) & \
                    (df['kernelB'] == kernelB) & \
                    (df['flopsB'] > flopsB*10**6)          
        elif ineqA == ">" and ineqB == "<":
            condition = (df['kernelA'] == kernelA) & \
                    (df['flopsA'] > flopsA*10**6) & \
                    (df['kernelB'] == kernelB) & \
                    (df['flopsB'] < flopsB*10**6)
        
        return df[condition]['operands'].to_list()
    
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
    
    def _prepare_kernel_stats(self):
        sa = self._df_k.groupby(["kernel", "class"])['class'].count().unstack(fill_value=0)
        sa = self._add_missing_class_data(sa)
        sa = self._compute_class_statistics(sa)
        
        if isinstance(self._df_ext, pd.DataFrame):
            sb = pd.DataFrame(self._df_ext.groupby(['kernel'])['case:rel-flops'].min()).reset_index()
            self.kernel_stats = sb.merge(sa, on='kernel')
        else:
            self.kernel_stats = sa
        
    def _prepare_relations_stats(self):    
        sa = self._df_r.groupby(["kernelA", "kernelB", "class"])['class'].count().unstack(fill_value=0)
        sa = self._add_missing_class_data(sa)
        self.relations_stats = self._compute_class_statistics(sa)
        
    