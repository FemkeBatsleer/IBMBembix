'''
Created on 6 dec. 2017

@author: fbatslee
'''
import Bembix_model_simple
import time
import numpy as np
import sys
start_time=time.clock()

for i in range(int(sys.argv[1]), int(sys.argv[2])):
    #sample params
    active = np.random.choice([True, False])
    CA_w_env = np.random.choice([True, False])
    CA_day = np.random.choice([True, False])
    response_func = np.random.choice([Bembix_model_simple.response_type2, Bembix_model_simple.response_par])
    weight_day = np.random.normal(10,3)
    while weight_day <= 0:
        weight_day = np.random.normal(10,3)
        
    choosers_on_followers = np.random.normal(0.559, 0.1)#np.random.sample()#the relative amount of choosers and followers [0,1)
    print(choosers_on_followers)
   
    pop5 = Bembix_model_simple.Population(amount_ind=455, amount_days=30, choosers_on_followers=choosers_on_followers,
                                   active=active, CA_w_env=CA_w_env, CA_day=CA_day, response_func=response_func,
                                   weight_day=weight_day)
    pop5.let_it_run()
    pop5.create_output(i, savepath='$VSC_DATA/Outputs/CF')
    
    print('CF population {}'.format(i))

print(time.clock()-start_time, 'seconds')