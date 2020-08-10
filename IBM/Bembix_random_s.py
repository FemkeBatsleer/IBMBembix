import Bembix_model_simple
import time
import sys
start_time=time.clock()

for i in range(int(sys.argv[1]), int(sys.argv[2])):
    pop1 = Bembix_model_simple.Population(amount_ind=455, amount_days=30, scenario='random')
    pop1.let_it_run()
    pop1.create_output(number=i, savepath='$VSC_DATA/Outputs/random')
    print('Random population {}'.format(i))

print(time.clock()-start_time, 'seconds')