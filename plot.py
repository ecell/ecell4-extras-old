import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import numpy
import sys

data = numpy.loadtxt(sys.argv[1])
# data = (data[-90:].T + data[-180:-90].T + data[-270:-180].T) / 3.0

newdata = []
for i in range(9):
    tmp = data[len(data) - 90 + i * 10: len(data) - 90 + (i + 1) * 10].T
    newdata.append(((i + 0.5) * 10, sum(tmp[2]), sum(tmp[3])))
newdata = numpy.array(newdata).T
print newdata

data = data[-90:].T
plt.plot(data[1], data[2], '-', color='#a0ffa0')
plt.plot(data[1], data[3], '-', color='#ffa0a0')
plt.bar(newdata[0] - 3, newdata[1] * 0.1, width=3, color='g', label='A1')
plt.bar(newdata[0], newdata[2] * 0.1, width=3, color='r', label='A2')
plt.xlim(0, 89)
plt.ylim(0, 22.5)
plt.legend(loc='best', shadow=True)
plt.ylabel("The Number of Molecules")
plt.xlabel("Coordinate along x-axis")
plt.savefig('data.eps')
