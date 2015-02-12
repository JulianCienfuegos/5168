from numpy import log, linspace
import matplotlib.pyplot as plt

alpha = 1.0/12.0
el1 = lambda x: 2.0*x
el2 = lambda x: 2.0*(1.0 - x)
u = lambda x: log(1+x)/(log(2))-x

x = linspace(0, 1, 20)
x1 = linspace(0, 0.5, 10)
x2 = linspace(0.5, 1, 10)

plt.plot(x1, alpha*el1(x1),'r', x2, alpha*el2(x2), 'r', x, u(x), 'b')
plt.title('Problem 1a')
plt.show()
