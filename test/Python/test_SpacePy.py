from spacepy import pybats
#import timeit

import time

t0 = time.time()
data = pybats.IdlFile("../3d_var_region0_0_t00001205_n00037679.out")
t1 = time.time()
total = t1-t0
print(total)

#Bx = data['Bx']
