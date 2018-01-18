#! /usr/bin/env python

import numpy as np
import sys
import pprint
from sklearn.cluster import dbscan

pars = sys.argv
if len(pars) != 5:
	print("\tUsage: ./clusterer.py [filename] [epsilon] [min_samples] [output]")
	sys.exit(-1)

pp = pprint.PrettyPrinter(indent=2)

# Notable SO goodie:
# https://stackoverflow.com/questions/3844801/check-if-all-elements-in-a-list-are-identical

# List of pixels:
#   pixel: [x, y, pID, MCl]
pixels = np.loadtxt(pars[1])
sim_particles = {x[3]:x[2] for x in pixels}
core_samples, cl_labels = dbscan(pixels, eps=float(pars[2]), min_samples=int(pars[3]), metric='euclidean')

# List of clusters:
#  cluster: { ClID: [pixel, ...] }
clusters = {}

# Fill a dictionary with clusters data
for pos,label in enumerate(cl_labels):
  # Dictionary of clusters does not contain any entry with corresponding cl_label
  if label not in clusters:
    clusters[label] = [list(pixels[pos])]
  # Label already present, append new pixel to pixel list
  else:
    clusters[label].append(list(pixels[pos]))

# Skim data:
#   Create on the fly the list of mc labels of pixels and test whether they contain
#   all the same labels
bad_clusters = {k:v for k,v in clusters.iteritems() if not [x[3] for x in v][1:] == [x[3] for x in v][:-1]}
hom_clusters = {k:v for k,v in clusters.iteritems() if [x[3] for x in v][1:] == [x[3] for x in v][:-1]}
frg_clusters = {}

# Iter over skimmed, find duplicates
for k1,v1 in hom_clusters.iteritems():
  for k2,v2 in hom_clusters.iteritems():
    if v1[0][3] == v2[0][3] and k1 != k2:
      frg_clusters[k1] = v1

# Remove fragmented from homogeneous, find good clusters
good_clusters = {k:v for k,v in hom_clusters.iteritems() if not k in frg_clusters.keys()}

# pp.pprint(hom_clusters)
# pp.pprint(bad_clusters)
# pp.pprint(frg_clusters)
# pp.pprint(good_clusters)

rec_pions = { k:v for k,v in hom_clusters.iteritems() if v[0][2] == 211 }
rec_heliums = { k:v for k,v in hom_clusters.iteritems() if v[0][2] == 1000020030 }
good_pions = { k:v for k,v in good_clusters.iteritems() if v[0][2] == 211 }
good_heliums  = { k:v for k,v in good_clusters.iteritems() if v[0][2] == 1000020030 }

print("Clusters Found:\n\tTotal:       %d\n\tHomogeneous: %d\n\tBad:         %d\n\tFragmented:  %d\n\tGood:        %d" % (len(clusters), len(hom_clusters), len(bad_clusters), len(frg_clusters), len(good_clusters)))
print("Particles:\n\tTotal MC:      %d\n\tTotal pions:   %d\n\tTotal Heliums: %d" % (len(sim_particles), len({ k:v for k,v in sim_particles.iteritems() if v == 211 }), len({k:v for k,v in sim_particles.iteritems() if v == 1000020030})))
print("\n\tTotal RC:      %d\n\tTotal pions:   %d\n\tTotal Heliums: %d" % (len(clusters), len(rec_pions), len(rec_heliums)))
print("\n\tTotal good RC:      %d\n\tTotal good pions:   %d\n\tTotal good Heliums: %d" % (len(good_clusters), len(good_pions), len(good_heliums)))

with open(pars[4],"a") as of:
  of.write("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%d\n" % 
    (len(sim_particles), len({ k:v for k,v in sim_particles.iteritems() if v == 211 }), len({k:v for k,v in sim_particles.iteritems() if v == 1000020030}), len(clusters), len(good_clusters), len(good_pions), len(good_heliums), float(pars[2]), float(pars[3])))
