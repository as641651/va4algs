
from am4pa.linnea import MeasurementsMaganer
import os
import json
from .operands_sampler_neighbour import OperandsSamplerNeighbour 



class SmartMeasurementManager:
    def __init__(self,rdl, adj_risk_thresh=0.05):
        self.rdl = rdl
        self.dml = rdl.dml
        self.ranking_method = self.rdl.rm.method
        
        ret = self.rdl.load()
        assert ret != -1, "First rank data"
        print("Ranking data from {} has been loaded.".format(self.rdl.obj_path))
        
        self.history = {'focus':[],'dirty':False}
        self.history_file = os.path.join(self.dml.lc.local_dir,
                                         'history_{}.json'.format(self.rdl.obj_path.split('/')[-1].split('.pkl')[0]))
        if os.path.exists(self.history_file):
            with open(self.history_file, 'r') as jf:
                self.history = json.load(jf)
            
        
        self.im = None
        self.focus = []
        
        
        self.num_ops = len(self.rdl.data_anomalies['op_str'].tolist()[0].split('_'))
        self.osn = OperandsSamplerNeighbour(self.num_ops,delta=50)
        self.mm = MeasurementsMaganer(self.dml,self.osn,self.rdl.thread_str)
        
        self.adj_risk_thresh = adj_risk_thresh
        _ = self.filter_interesting_operands()
        
    
    def filter_interesting_operands(self):
        df = self.rdl.data_anomalies
        self.im = df[df['adj_risk']>self.adj_risk_thresh].sort_values(['rel-flops-cutoff','adj_risk'], ascending=[False,False])
            
            
        focus_ops = self.im[self.im['adj_risk']>self.adj_risk_thresh]['op_str'].tolist()
        
        self.focus = []
        for op_str in focus_ops:
            if op_str not in self.history['focus']:
                self.focus.append(op_str)
        
        return self.im
        
    def _smart_generate_variants(self):
        for op_str in self.focus:
            print("Focus: {}".format(op_str))
            self.osn.set_focus(list(map(int,op_str.split('_'))))
            self.mm.generate_variants_sampler(5)
           
    def generate_measure(self,reps,run_id,bSlrum):
        if self.history['dirty']:
            print("Call rank data before further measurements")
            return -1
        
        if not self.focus:
            print("Nothing to focus")
            return -1
            
        self._smart_generate_variants()
        self.mm.measure_variants(reps,run_id=run_id,bSlrum=bSlrum)
        self.history['dirty'] = True
        self.dml._update_json(self.history,self.history_file)
            
        return 1
        
    def check_slrum_status(self):
        self.dml.check_completed_slrum_jobs()
        if not self.dml.slrum_running_jobs['r']:
            print("Completed")
            return 1
        else:
            print(self.dml.slrum_running_jobs['r'])
            return -1
        
    def rank_update(self):
        if self.history['dirty']:
            self.rdl.rank3way(update=True)
            for op in self.focus:
                if not op in self.history['focus']:
                    self.history['focus'].append(op)
            self.history['dirty'] = False
            self.dml._update_json(self.history,self.history_file)
            self.rdl.save()
        else:
            print("No new measurements")
            
    def stop_search_path(self):
        if not self.history['dirty']:
            for op in self.focus:
                if not op in self.history['focus']:
                    self.history['focus'].append(op)
            self.dml._update_json(self.history,self.history_file)
            

