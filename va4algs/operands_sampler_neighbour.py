from am4pa.linnea import OperandsSamplerBase
import random

class OperandsSamplerNeighbour(OperandsSamplerBase):
    def __init__(self,num_operands,delta=50,seed=108):
        super().__init__(num_operands)
        self.seed = seed
        random.seed(self.seed)
        self.delta = delta
        self.focus = None
        
    def set_focus(self,focus):
        self.focus = focus
    
    def sample(self):
        
        if not self.focus:
            return
        
        ops = []
        
        for i in range(self.n):
            
            min_ = max(10,self.focus[i]-self.delta)
            max_ = self.focus[i]+self.delta
            
            ops.append(random.randint(min_,max_))
        return ops