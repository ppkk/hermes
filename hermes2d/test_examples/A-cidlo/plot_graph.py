# import libraries
import numpy, pylab
from pylab import *

# plot DOF convergence graph
pylab.title("Point value")
pylab.xlabel("relative permitivity")
pylab.ylabel("value")
#axis('equal')
pylab.grid(True)
ax = pylab.gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(20)
data = numpy.loadtxt("pic/normal_calculations.dat")
x = data[:, 0]
y = data[:, 1]
plot(x, y, 's', label="exact")

data = numpy.loadtxt("pic/pgd_calculations.dat")
x = data[:, 0]
for i in range(1, data[0].size):
    y = data[:, i]
    plot(x, y, '-', label=str(i))    

legend()

# finalize
#show()
savefig("pic/dependence.png")
pylab.close()

#
# plot parameters separately 
#
for i in range(1, data[0].size):    
    filename = "pic/parameter00%d.dat" % (i-1)
    picturename = "pic/parameter00%d.png" % (i-1)
    data2 = numpy.loadtxt(filename)
    x2 = data2[:, 0]
    y2 = data2[:, 1]
    plot(x2, y2, '-', label="$p_%d$" % i)
    legend()
    savefig(picturename)
    pylab.close()

#
# plot parameters to one file 
#
picturename = "pic/parameters.png"
for i in range(1, data[0].size):    
    filename = "pic/parameter00%d.dat" % (i-1)
    data2 = numpy.loadtxt(filename)
    x2 = data2[:, 0]
    y2 = data2[:, 1]
    plot(x2, y2, '-', label="$p_%d$" % i)
legend()
savefig(picturename)
pylab.close()


pylab.title("Convergence")
pylab.xlabel("iteration")
pylab.ylabel("function difference")

lines = [np.array(map(double, line.split())) for line in open('pic/convergence.dat')]
i = 1
for line in lines:
    x = range(1, len(line) + 1)
    y = list(line)
    semilogy(x, y, '-', label="step %d" % i)
    i += 1

legend()
savefig("pic/convergence.png")
pylab.close()
    