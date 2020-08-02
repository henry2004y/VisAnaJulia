from spacepy import pybats
#import timeit

import time

t0 = time.time()
#data = pybats.IdlFile("3d_var_region0_0_t00001527_n00005528.out")
data = pybats.IdlFile("../3d_var_region0_0_t00001205_n00037679.out")
#data = pybats.IdlFile("/Users/hyzhou/Documents/Computer/Julia/SWMF/test/z=0_raw_1_t25.60000_n00000258.out")
t1 = time.time()
total = t1-t0
print(total)

#Bx = data['Bx']
