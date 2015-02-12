from numpy import linspace, array, polyval
import matplotlib.pyplot as plt

xi1 = -1.0
xi2 = -1.0/3.0
xi3 = 1.0/3.0
xi4 = 1.0

psi = lambda x, xi, xi_i, xi_j, xi_k: ((x-xi_i)*(x-xi_j)*(x-xi_k))/\
	  ((xi - xi_i)*(xi - xi_j)*(xi - xi_k))
psi1 = lambda x: psi(x, xi1, xi2, xi3, xi4)
psi2 = lambda x: psi(x, xi2, xi1, xi3, xi4)
psi3 = lambda x: psi(x, xi3, xi1, xi2, xi4)
psi4 = lambda x: psi(x, xi4, xi1, xi2, xi3)

p1 = array([ -9.0/16,  9.0/16,   1.0/16, -1.0/16])
p2 = array([ 27.0/16, -9.0/16, -27.0/16,  9.0/16])
p3 = array([-27.0/16, -9.0/16,  27.0/16,  9.0/16])
p4 = array([  9.0/16,  9.0/16,  -1.0/16, -1.0/16])

x = linspace(-1, 1, 20)

plt.plot(x, psi2(x), 'y', x, psi3(x), 'b', x, psi4(x), 'g',x, psi1(x), 'r' )
plt.plot(x, polyval(p1, x), '-', x, polyval(p2, x),'o', x, polyval(p3, x), '*', x, polyval(p4,x),'x')
plt.title('Problem 2')
plt.show()
