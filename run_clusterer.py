#! /usr/bin/env python

import os

for i in range(1,1001):
  print("Iteration number: %d" % i)
  os.system("./clusterer.py data.txt %f 1 epsilon.txt > /dev/null " % float(i*0.0015))
