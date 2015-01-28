from scipy import sin, cos, sinh
from scipy.integrate import quad as di
from numpy import pi, linspace, array
from numpy.linalg import solve as sls
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math

# Test Functions.
phi    = lambda x, c: sin(c*pi*x)
d_phi  = lambda x, c: c*pi*cos(c*pi*x)
phi2   = lambda x, c1, c2: phi(x, c1) * phi(x, c2)
d_phi2 = lambda x, c1, c2: d_phi(x, c1) * d_phi(x, c2)
f      = lambda x: x

# Parameters 
m = 0
M = 1
X    = linspace(m, M, 20)

# Real Solution
u = lambda x: x - sinh(x)/sinh(1)

# These functions generate matrix elements.
def K_(a, b, i, j):
    f = lambda x: phi2(x, i, j) + d_phi2(x, i, j)
    return di(f, a, b)[0]    

def F_(a, b, i):
    F = lambda x: f(x)*phi(x,i)
    return di(F, a, b)[0]

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
    col_labels=['Abs Error']
    row_labels=['x=0.25','x=0.5', 'x=0.75']
    table_vals=[[math.ceil(prec*err[0])/prec],[math.ceil(prec*err[1])/1000000],[math.ceil(prec*err[2])/prec]]
    # the rectangle is where I want to place the table
    the_table = plt.table(cellText=table_vals,
                      colWidths = [0.15],
                      rowLabels=row_labels,
                      colLabels=col_labels,
                      loc='lower right')
    plt.savefig(img_name)
    
# N = 1
alpha = F_(m, M, 1)/K_(m, M, 1, 1)
u_approx = lambda x: alpha*phi(x,1)
err = [abs(u(0.25) - u_approx(0.25)), abs(u(0.5) - u_approx(0.5)), abs(u(0.75) - u_approx(0.75))]
plot_func('One Basis Func', 'one.png', u, u_approx, X, err)

# N = 2
K2 = array([[K_(m, M, 1, 1), K_(m, M, 1, 2)],
            [K_(m, M, 2, 1), K_(m, M, 2, 2)]])
F2 = array([F_(m, M, 1),
            F_(m, M, 2)])
            
alpha = sls(K2, F2)
u_app = lambda x: alpha[0]*phi(x, 1) + alpha[1]*phi(x, 2)
err = [abs(u(0.25) - u_app(0.25)), abs(u(0.5) - u_app(0.5)), abs(u(0.75) - u_app(0.75))]
plot_func('Two Basis Functions', 'two.png', u, u_app, X, err)

# N = 3
K3 = array([[K_(m, M, 1, 1), K_(m, M, 1, 2), K_(m, M, 1, 3)],
            [K_(m, M, 2, 1), K_(m, M, 2, 2), K_(m, M, 2, 3)],
            [K_(m, M, 3, 1), K_(m, M, 3, 2), K_(m, M, 3, 3)]])

F3 = array([F_(m, M, 1),
            F_(m, M, 2),
            F_(m, M, 3)])
            
alpha = sls(K3, F3)
u_a = lambda x: alpha[0]*phi(x, 1) + alpha[1]*phi(x, 2) + alpha[2]*phi(x,3)
err = [abs(u(0.25) - u_a(0.25)), abs(u(0.5) - u_a(0.5)), abs(u(0.75) - u_a(0.75))]
plot_func('Three Basis Functions', 'three.png', u, u_a, X, err)


    
    
