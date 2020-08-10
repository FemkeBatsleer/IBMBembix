'''
Created on 5 dec. 2017

@author: fbatslee
'''
import Bembix_model_simple
import time
import numpy as np
import sys
start_time=time.clock()

for i in range(int(sys.argv[1]), int(sys.argv[2])):
    time_step = time.clock()
    pop4 = Bembix_model_simple.Population(amount_ind=455, amount_days=30, scenario='prev-exp-BU')
    pop4.let_it_run()
    pop4.create_output(i, savepath='$VSC_DATA/Outputs/PE')
    print('PEBU population {}'.format(i))
    print(time.clock()-time_step)
    

print(time.clock()-start_time, 'seconds')