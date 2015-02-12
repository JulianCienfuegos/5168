from numpy import log, linspace, array
from numpy.linalg import solve
import matplotlib.pyplot as plt

# Part A
alpha = 1.0/12.0
el1 = lambda x: 2.0*x
el2 = lambda x: 2.0*(1.0 - x)
u = lambda x: log(1+x)/(log(2))-x

x = linspace(0, 1, 20)
x1 = linspace(0, 0.5, 10)
x2 = linspace(0.5, 1, 10)

plt.subplot(2, 1, 1)
plt.plot(x1, alpha*el1(x1),'r', x2, alpha*el2(x2), 'r', x, u(x), 'b')
plt.title('Problem 1a')

# Part B
K = array([[10, -5.5, 0],[-5.5, 12, -6.5],[0, -6.5, 14]])
f = 0.25
F = array([f, f, f])
alpha =  solve(K,F)

el1 = lambda x: alpha[0]*4*x
el2 = lambda x: alpha[1]*4*(x - 0.25) + alpha[0]*4*(0.5 - x)
el3 = lambda x: alpha[1]*4*(0.75 - x) + alpha[2]*4*(x - 0.5)
el4 = lambda x: alpha[2]*4*(1 - x)
u   = lambda x: log(1+x)/(log(2))-x

x  = linspace(0,1,20)
x1 = linspace(0, 0.25, 10)
x2 = linspace(0.25, 0.5, 10)
x3 = linspace(0.5, 0.75, 10)
x4 = linspace(0.75, 1, 10)

plt.subplot(2, 1, 2)
plt.plot(x1, el1(x1), 'r', x2, el2(x2), 'r', x3, el3(x3), 'r', x4, el4(x4), 'r', x, u(x), 'b')
plt.title('Problem 1b')
plt.show()


