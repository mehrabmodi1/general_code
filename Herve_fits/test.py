# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 11:13:39 2020

@author: Mehrab
"""

import numpy as np
act_fn = "blah/blah/blah/{}/fly{}/dFF_data.mat".format('Gamma', '1');

print(act_fn)
lobe = 'Gamma1';
fly = 1;

folder_fit = "C:/Data/Data/Analysed_data/data_sharing/KC_param_fitting_sandbox/fits_{}_fly{}".format(lobe, fly)

for a in np.arange(5):
    print(a)
    
    
import os as os    
from pathlib import Path, PureWindowsPath
base_path = PureWindowsPath("C:\Data\Data\Analysed_data\data_sharing\KC_long_trace");
base_path = Path(base_path);
obj_list = os.listdir(base_path);
n_objs = len(obj_list);
n_dirs = 0;
for obj_n in range(0, (n_objs - 1) ):
    curr_path = "{}/{}".format(base_path, obj_list[obj_n]);
    if os.path.isdir(curr_path) == 1:
        n_dirs = n_dirs + 1;


        
import matplotlib as plt
test_var = np.random.rand(1601,)
test_var[800:1600,] = np.nan;
#plt.pyplot.imshow(test_var)
breakpoint()
plt.pyplot.plot(test_var)