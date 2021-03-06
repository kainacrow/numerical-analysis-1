#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 19:37:52 2018

@author: Kaina
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 12:48:09 2018

@author: Kaina
"""
import numpy as np
import matplotlib.pyplot as plt

f = lambda x,t: x - math.sin(t)
x0 = 0

x = lambda t: -0.5*((np.e**t))+(0.5*((np.sin(t))+(np.cos(t))))

def calculateError(tArray, xArray):
    e = 0.0
    for i in range(len(tArray)):
        if(np.abs(x(tArray[i]) - xArray[i]) >e):
            e = np.abs(x(tArray[i])-xArray[i])
    return e

N = 10
delta = np.pi/N
tValues = np.arange(0, np.pi, delta)
xValues = x(tValues)

x0,x1 = 0.0, 0.0
        
def trapezoid(x0, N, delta):
    xnplus1 = lambda xn, n: xn*((1+(0.5*delta))/(1-(0.5*delta)))+((0.5*delta)/(0.5*delta)-1)*(np.sin(delta*(n+1))+np.sin(delta*(delta*n)))
    tValues = [0.0]
    xValues = [x0]
    for n in np.arange(delta, N*delta+delta, delta):
        tValues.append(n)
        xValues.append(xnplus1(xValues[-1], tValues[-1]))
    return xValues[-1]

x1 = trapezoid(x0, delta, delta/N)
    
def AB2(x0, x1, N, delta):
    xnplus1 = lambda xn, tn, xnminus1, tnminus1: xn + delta/2.0 * (3*(f(xn, tn))-(f(xnminus1, tnminus1)))
    tValues = [0.0, delta]
    xValues = [x0, x1]
    for n in np.arange(2*delta, np.pi+delta, delta):
        tValues.append(n)
        xValues.append(xnplus1(xValues[-1], tValues[-1], xValues[-2], tValues[-2]))
    return [np.array(tValues), np.array(xValues)]
    
ts,xs = AB2(x0, x1, N, delta)
calculateError(ts, xs)
#print(delta)
#print(calculateError(ts, xs))
#for i in range(len(ts)):
#    print((ts[i], xs[i]))
    
## first order method??
    
Nlist = [int(np.sqrt(10)**(i+2)) for i in range(11)]
#print Nlist
errorList = []

for N in Nlist:
    delta = np.pi/N
    x1 = trapezoid(x0, delta, delta/N)
    ts, xs = AB2(x0, x1, N, delta)
    error = calculateError(ts, xs)
    errorList.append(error)
  
print errorList

#for i in range(len(Nlist)):
#    print((Nlist[i], errorList[i]))

ln = lambda x: np.log(x) 
print np.polyfit(ln(Nlist[2:]), ln(errorList[2:]),1)    
    
plt.figure(1,figsize=(20, 15))
plt.subplot(222)
plt.title('Trapezoid Method')
plt.plot(Nlist, errorList,label = 'N vs Error')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.grid(True)
plt.show()
#    