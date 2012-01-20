#/usr/bin/python

import json, time
import subprocess, os, sys
import matplotlib.pyplot as plt

data = open(sys.argv[1])

i = 0

for line in data.readlines():
    i = i + 1
    line = line.split("\t")
    points = line[1:]
    x = []
    y = []
    for idx, val in enumerate(points):
    	if idx%2 == 0:
		x.append(val)
	else:
		y.append(val)
    plt.plot(x,y,'b.')
    com = line[0].split(" ")
    plt.plot([com[0]], [com[1]], "rs")
    plt.axis((float(com[0])-2, float(com[0])+2,float(com[1])-2, float(com[1])+2))

    filename = str('%04d' % i) + '.png'
    plt.savefig(filename, dpi=100)

    print 'Wrote file', filename

    #
    # Clear the figure to make way for the next image.
    #
    plt.clf()


sys.exit(0)
