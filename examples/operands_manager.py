import pandas as pd
from typing import List

class OperandsManager:
    def __init__(self):
        self.df = pd.DataFrame(columns=['op', 'generated', 'measured'])
        
    def op_str(self, op:List):
        return "_".join(op)
        
    def add_operands(self, ops:List):
        for op in ops:
            op_str = self.op_str(op)
            if not op_str in self.df['op'].values:
                row = {'op':op_str, 'generated': False, 'measured':False}
                self.df = self.df.append(row, ignore_index=True)
                
    def set_generated(self, op:List):
        self.df.loc[self.df['op'] == self.op_str(op), 'generated'] = True
        
    def is_generated(self, op:List):
        return self.df.loc[self.df['op'] == self.op_str(op), 'generated'].values[0]
    
    def set_measured(self, op:List):
        self.df.loc[self.df['op'] == self.op_str(op), 'measured'] = True
        
    def is_measured(self, op:List):
        return self.df.loc[self.df['op'] == self.op_str(op), 'measured'].values[0]