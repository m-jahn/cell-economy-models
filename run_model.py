# Examples to import and run various cellular economy models
# ----------------------------------------------------------
#
# Synechocystis sp. PCC6803 (cyanobacterium, photosynthetic)
# ----------------------------------------------------------

# 1. import libraries and model(s)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plth
from models import synechocystis_steadystate
from models import synechocystis_dynamic



# 2. define initial parameters
sub = 100                                       # initial substrate concentration, CO2/HCO3-
Ki = 5000                                       # light inhibition constant for photosystems
mumax = 0.11                                    # maximum growth rate, used to calculate protein utilization
c_upper = [1,1,1,1,1,1,90,25,25,5,1,1,1]        # optional list of concentration upper bounds
reserve = [0.0, 0.0, 0.0, 0.0, 0.0]             # fraction of enzyme reserve in total

# 3. define light conditions

# (A) light in % max intensity, log decrease
#light = 100.0/1.5**np.array(range(0,12))

# (B) light as step change
#light = np.array([5.0]*25+[50.0]*26)

# (C) light coming in pulses
#light = np.array([3.0]*12+[50.0]*6+[3.0]*12+[50.0]*6+[3.0]*13)

# (D) light as smooth day night cycle
# use sine function to simulate one full day at length 2*pi = 6.283,
# so 2 days is 4*pi, and step width = 4*pi/96,
# since sine(x) is between -1 and 1, we rescale by (sine(x)+1)*50 (0 to 100)
light = np.round((np.sin(np.arange(0, 4*3.1415, 4*3.1415/96))+1)*50)+1

# time as a function of light step number, in hours
time = np.arange(0, len(light)/2, 0.5)


# 3. run model simulations
#    loop through different values of ribosome reserve
for i in [0.0, 0.05, 0.10, 0.15]:
    reserve[4] = i
    result_ss = synechocystis_steadystate.simulate(time, light, sub, c_upper, reserve, mumax, Ki)
    result_ss.table.to_csv('results/synechocystis/daynight/steady_state_RIB' + str(i) + '.csv')
    result_dy = synechocystis_dynamic.simulate(time, light, sub, c_upper, reserve, mumax, Ki, a = result_ss.a, c = result_ss.c)
    result_dy.table.to_csv('result_dynamic_RIB' + str(i) + '.csv')


# 4. results visualization






plt.figure()
# growth rate/ribosomes
plt.subplot(2,2,1)
plt.axis([0,100,0,0.2])
plt.plot(m.time, v['RIB'].value,'g-',label='mu')
# enzyme rates
plt.subplot(2,2,2)
for i in enz:
    plt.plot(m.time, v[i].value,'b--',label='r')

# protein alpha optima and concentration
plt.subplot(2,2,3)
plt.axis([0,100,0,0.32])
for i in pro:
    plt.plot(m.time, a[i].value,'g-',label='a')
    plt.plot(m.time, c[i].value,'r--',label='c')

plt.axis([0,100,0,0.32])

# concentration metabolites
plt.subplot(2,2,4)
for i in met:
    plt.plot(m.time, c[i].value,'r--',label='c')

plt.show()

