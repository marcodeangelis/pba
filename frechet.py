"""
Created on Wed Apr 15 2021
@author: Marco De Angelis
University of Liverpool
github.com/marcodeangelis
MIT License

# Code the computation of Frechet bounds with intrusive arithmetic.

To see how this code can be used in practice, check the Readme.md file in this repository. 



"""


from matplotlib import pyplot
import numpy


def interval(lo,hi):
    return [lo, hi]

def pbox(leftcdf,rightcdf,steps=50, unbounded=True, interval=False):
    a_l = leftcdf
    a_r = rightcdf
    p = numpy.linspace(0.001,0.999,steps+1)
    if interval:
        p_l = list(p).copy()
        p_r = list(p).copy()
        x_l = [a_l] * (steps+1)
        x_r = [a_r] * (steps+1)
    else:
        p_r = list(p).copy()
        p_l = list(p[1:])+list([p[-1]])
        x_l = list(a_l.ppf(p))
        x_r = list(a_r.ppf(p))
    
    l = [x_l,p_l]
    r = [x_r,p_r]
    return l,r
    
def plot(X,ax=None):
    def stepdata(x,y): # x,y must be python lists
        xx,yy = x*2, y*2
        xx.sort()
        yy.sort()
        return xx, [0.]+yy[:-1]
    if ax is None:
        fig,ax = pyplot.subplots(figsize=(12,9))
    x,y=stepdata(X[0][0],X[0][1])
    ax.plot(x,y)
    x,y=stepdata(X[1][0],X[1][1])
    ax.plot(x,y)



# ----- Arithmetic operations start here ----- #

# X[0][0]: x-values left cdf bound 
# X[0][1]: p-values left cdf bound 

# X[1][0]: x-values right cdf bound 
# X[1][1]: p-values right cdf bound 

def frechet_sum(X,Y):
    p_Z_hi = X[0][1]
    p_Z_lo = X[1][1]
    X_lo, Y_lo = X[1][0][:-1], Y[1][0][:-1]
    X_hi, Y_hi = X[0][0][1:],  Y[0][0][1:]
    N = len(X_lo)
    Z_lo, Z_hi = [], []
    for i in range(N):
        z_lo, z_hi = [], []
        for j in range(i,N):
            z_lo.append(X_lo[j]+Y_lo[i-j+N-1])
        if i == 0:
            z_hi.append(X_hi[0]+Y_hi[0])
        else:
            for j in range(i):
                z_hi.append(X_hi[j]+Y_hi[i-j])
        Z_lo.append(min(z_lo))
        Z_hi.append(max(z_hi))

    Z = (([Z_hi[0]]+Z_hi,p_Z_hi),(Z_lo+[Z_lo[-1]],p_Z_lo))
    return Z

def frechet_minus(X,Y):
    p_Z_hi = X[0][1]
    p_Z_lo = X[1][1]
    X_lo, Y_lo = X[1][0][:-1], Y[1][0][:-1]
    X_hi, Y_hi = X[0][0][1:],  Y[0][0][1:]
    N = len(X_lo)
    Z_lo, Z_hi = [], []
    for i in range(N):
        z_lo , z_hi = [], []
        for j in range(i,N):
            z_lo.append(X_lo[j]-Y_hi[j-i])
        if i == 0:
            z_hi.append(X_hi[0]-Y_lo[N-1])
        else:
            for j in range(i):
                z_hi.append(X_hi[j]-Y_lo[j-i+N-1])
        Z_lo.append(min(z_lo))
        Z_hi.append(max(z_hi))

    Z = (([Z_hi[0]]+Z_hi,p_Z_hi),(Z_lo+[Z_lo[-1]],p_Z_lo))
    return Z

def frechet_times(X,Y):
    p_Z_hi = X[0][1]
    p_Z_lo = X[1][1]
    X_lo, Y_lo = X[1][0][:-1], Y[1][0][:-1]
    X_hi, Y_hi = X[0][0][1:],  Y[0][0][1:]
    N = len(X_lo)
    Z_lo, Z_hi = [], []
    for i in range(N):
        z_lo , z_hi = [], []
        for j in range(i,N):
            z_lo.append(X_lo[j]*Y_lo[i-j+N-1])
        if i == 0:
            z_hi.append(X_hi[0]*Y_hi[0])
        else:
            for j in range(i):
                z_hi.append(X_hi[j]*Y_hi[i-j])
        Z_lo.append(min(z_lo))
        Z_hi.append(max(z_hi))

    Z = (([Z_hi[0]]+Z_hi,p_Z_hi),(Z_lo+[Z_lo[-1]],p_Z_lo))
    return Z

def frechet_div(X,Y):
    p_Z_hi = X[0][1]
    p_Z_lo = X[1][1]
    X_lo, Y_lo = X[1][0][:-1], Y[1][0][:-1]
    X_hi, Y_hi = X[0][0][1:],  Y[0][0][1:]
    N = len(X_lo)
    Z_lo, Z_hi = [], []
    for i in range(N):
        z_lo , z_hi = [], []
        for j in range(i,N):
            z_lo.append(X_lo[j]/Y_hi[j-i])
        if i == 0:
            z_hi.append(X_hi[0]/Y_lo[N-1])
        else:
            for j in range(i):
                z_hi.append(X_hi[j]/Y_lo[j-i+N-1])
        Z_lo.append(min(z_lo))
        Z_hi.append(max(z_hi))

    Z = (([Z_hi[0]]+Z_hi,p_Z_hi),(Z_lo+[Z_lo[-1]],p_Z_lo))
    return Z