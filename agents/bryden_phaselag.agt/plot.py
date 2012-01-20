import sys, os, re
import matplotlib.pyplot as plt

name = sys.argv[1].split(".")
number = int(name[1])
while True:
	x = []
	a = []
	b = []
	c = []
	d = []
	e = []
	path = "%s.%d.out"%(name[0], number)
	if not os.path.exists(path):
		break
	f = open(path)
	count = 0
	for line in f:
		if count > 4000 and count < 5000:
			line = line.split(" ")
			x.append(float(line[0]))
			a.append(float(line[int(sys.argv[2])]))
			b.append(float(line[int(sys.argv[3])]))
			c.append(float(line[int(sys.argv[4])]))
			d.append(float(line[int(sys.argv[5])]))
			e.append(float(line[int(sys.argv[6])]))
		if count > 5000:
			break
		count = count + 1
	plt.plot(x, a)
	plt.plot(x, b)
	plt.plot(x, c)
	plt.plot(x, d)
	plt.plot(x, e)
	number = number + 1
plt.show()
