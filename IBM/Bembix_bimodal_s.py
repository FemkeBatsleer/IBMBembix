'''
Created on 6 dec. 2017

@author: fbatslee
'''
import Bembix_model_simple
import time
import numpy as np
import sys
start_time=time.clock()

for i in range(int(sys.argv[1]), int(sys.argv[2])):#range(int(sys.argv[1]), int(sys.argv[2])):
    #sample params
    active = np.random.choice([True, False])
    CA_w_env = np.random.choice([True, False])
    CA_day = np.random.choice([True, False])
    response_func = np.random.choice([Bembix_model_simple.response_type2, Bembix_model_simple.response_par])
    bi_CA = np.random.choice([True, False])#defines, when a wasp makes its next nest closeby, if it also takes into account the surrounding nests there
    weight_day = np.random.normal(10,3)
    while weight_day <= 0:
        weight_day = np.random.normal(10,3)
    print(active, CA_w_env, CA_day, response_func, bi_CA, weight_day)
   
    pop5 = Bembix_model_simple.Population(amount_ind=455, amount_days=30, scenario='bimodal',
                                   active=active, CA_w_env=CA_w_env, CA_day=CA_day, response_func=response_func,
                                   bi_CA=bi_CA, weight_day=weight_day)
    pop5.let_it_run()
    pop5.create_output(i, savepath='data')
    
    print('Bimodal population {}'.format(i))

print(time.clock()-start_time, 'seconds')