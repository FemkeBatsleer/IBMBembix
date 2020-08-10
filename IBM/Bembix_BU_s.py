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

    pop2 = Bembix_model_simple.Population(amount_ind=455, amount_days=30, scenario='bottom-up')
    pop2.let_it_run()
    pop2.create_output(number=i, savepath='$VSC_DATA/data/Outputs/BU')
    print('BU population {}'.format(i))

print(time.clock()-start_time, 'seconds')