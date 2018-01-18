#! /usr/bin/env python

import os

#for n in [1,2,5,10]:
#  for i in range(1,201):
#    print("Iteration number: %d, epsilon: %f, %d" % (i, float(i*0.0090), n))
#    os.system("./clusterer.py data.txt %f %d epsilon_n%d.txt > /dev/null" % (float(i*0.0090), n, n))

for i in range(1,11):
  print("Iteration number with fixed epsilon and N=%d" % i)
  os.system("./clusterer.py data.txt 0.009 %d epsilon_fixed_0.009.txt > /dev/null " % i)
