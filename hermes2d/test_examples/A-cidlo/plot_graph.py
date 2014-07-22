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
