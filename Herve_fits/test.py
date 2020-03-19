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