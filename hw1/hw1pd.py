from scipy import sin, cos, sinh
from scipy.integrate import quad as di
from numpy import pi, linspace, array
from numpy.linalg import solve as sls
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math

# Test Functions.
p1  = lambda x: x*(1 - x)
p2  = lambda x: x*(1 - x)*(1.0/2.0 - x)
p3  = lambda x: x*(1 - x)*(1.0/3.0 - x)*(2.0/3.0 - x)
dp1 = lambda x: 1 - 2*x
dp2 = lambda x: 3*x**2 - 3*x + 1.0/2.0
dp3 = lambda x: 2.0/9.0 - 22.0*x/9.0 + 6*x**2 - 4*x**3

p11 = lambda x: p1(x)**2
p22 = lambda x: p2(x)**2
p33 = lambda x: p3(x)**2

p12 = lambda x: p1(x) * p2(x)
p13 = lambda x: p1(x) * p3(x)
p23 = lambda x: p2(x) * p3(x)

dp11 = lambda x: dp1(x)**2
dp22 = lambda x: dp2(x)**2
dp33 = lambda x: dp3(x)**2

dp12 = lambda x: dp1(x) * dp2(x)
dp13 = lambda x: dp1(x) * dp3(x)
dp23 = lambda x: dp2(x) * dp3(x)

f1 = lambda x: x*p1(x)
f2 = lambda x: x*p2(x)
f3 = lambda x: x*p3(x)

def dInt(fun):
    return di(fun, m, M)[0]
    
# Parameters 
m = 0
M = 1
X = linspace(m, M, 50)

# Real Solution
u = lambda x: x - sinh(x)/sinh(1)

def plot_func(title, img_name, f, g, X, err):
    prec = 10000000
    plt.figure()
    plt.plot(X, f(X),'r' , X, g(X),'g')
    plt.title(title)
    plt.xlabel('X')
    plt.ylabel('u')
    red_patch = mpatches.Patch(color='red', label='True Solution')
    green_patch = mpatches.Patch(color='green', label='Approximate Solution')
    plt.legend(handles=[red_patch, green_patch], loc = 1, prop={'size':10})
    ax=plt.gca()
    col_labels=['% Error']
    row_labels=['x=0.25','x=0.5', 'x=0.75']
    table_vals=[[math.ceil(prec*err[0])/prec],[math.ceil(prec*err[1])/prec],[math.ceil(prec*err[2])/prec]]
    # the rectangle is where I want to place the table
    the_table = plt.table(cellText=table_vals,
                      colWidths = [0.15],
                      rowLabels=row_labels,
                      colLabels=col_labels,
                      loc='lower right')
    plt.savefig(img_name)
    
# N = 1
i11 = lambda x: p11(x) + dp11(x)
alpha = dInt(f1)/dInt(i11)
u_approx = lambda x: alpha*p1(x)

err = [abs(u(0.25) - u_approx(0.25))/u(0.25), 
       abs(u(0.5) - u_approx(0.5))/u(0.5), 
       abs(u(0.75) - u_approx(0.75))/u(0.75)]
plot_func('1 Polynomial Basis Function', 'one_poly.png', u, u_approx, X, err)

# N = 2
i12 = lambda x: p12(x) + dp12(x)
i22 = lambda x: p22(x) + dp22(x)
K = array([[dInt(i11), dInt(i12)],
           [dInt(i12), dInt(i22)]])
F = array([dInt(f1),
           dInt(f2)])
            
alpha = sls(K, F)
u_approx = lambda x: alpha[0]*p1(x) + alpha[1]*p2(x)
err = [abs(u(0.25) - u_approx(0.25))/u(0.25),
       abs(u(0.5) - u_approx(0.5))/u(0.5),
       abs(u(0.75) - u_approx(0.75))/u(0.75)]
plot_func('2 Polynomial Basis Functions', 'two_poly.png', u, u_approx, X, err)


# N = 3
i13 = lambda x: p13(x) + dp13(x)
i23 = lambda x: p23(x) + dp23(x)
i33 = lambda x: p33(x) + dp33(x)
K = array([[dInt(i11), dInt(i12), dInt(i13)],
           [dInt(i12), dInt(i22), dInt(i23)],
           [dInt(i13), dInt(i23), dInt(i33)]])

F = array([dInt(f1),
           dInt(f2),
           dInt(f3)])
                       
alpha = sls(K, F)
u_approx = lambda x: alpha[0]*p1(x) + alpha[1]*p2(x) + alpha[2]*p3(x)
err = [abs(u(0.25) - u_approx(0.25))/u(0.25),
       abs(u(0.5) - u_approx(0.5))/u(0.5),
       abs(u(0.75) - u_approx(0.75))/u(0.75)]
plot_func('Three Polynomial Basis Functions', 'three_poly.png', u, u_approx, X, err)
    
    
