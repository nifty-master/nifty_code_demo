import numpy as np
import torch
import ot
#import pandas as pd
#import matplotlib.pyplot as plt
##data1 = torch.tensor(pd.read_csv("ynew").values)
##data2 = torch.tensor(pd.read_csv("ytest").values) 
data1 = np.loadtxt("ynew", delimiter=',')
data2 = np.loadtxt("ytest",delimiter=',') 
n = data1.shape[0] 
p = data1.shape[1]  
W1 = ot.sliced_wasserstein_distance(data1, data2, seed=0) 