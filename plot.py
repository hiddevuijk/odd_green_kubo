import numpy as np
import matplotlib.pyplot as plt
from sys import exit

dt = 1e-1

t = np.loadtxt("t.dat") * dt
cxx = np.loadtxt("cxx.dat")
cyy = np.loadtxt("cyy.dat")
cxy = np.loadtxt("cxy.dat")
cyx = np.loadtxt("cyx.dat")

cvxx= np.loadtxt("cvxx.dat")
cvyy= np.loadtxt("cvyy.dat")
cvxy= np.loadtxt("cvxy.dat")
cvyx= np.loadtxt("cvyx.dat")

cvxvx= np.loadtxt("cvxvx.dat")
cvyvy= np.loadtxt("cvyvy.dat")
cvxvy= np.loadtxt("cvxvy.dat")
cvyvx= np.loadtxt("cvyvx.dat")



#plt.scatter(t[1:], cxx[1:] / (2 * t[1:]), label="cxx")
#plt.scatter(t[1:], cyy[1:] / (2 * t[1:]), label="cyy")
#plt.scatter(t[1:], cxy[1:] / (2 * t[1:]))
#plt.scatter(t[1:], cyx[1:] / (2 * t[1:]))

plt.scatter(t, cvxx, label="cvxx")
#plt.scatter(t, cvyy, label="cvyy")
plt.scatter(t, cvxy, label="cvxy")
#plt.scatter(t, cvxy, label="cvxy")

plt.scatter(t, cvxvx, label="cvxvx")
#plt.scatter(t, cvyvy, label="cvyvy")
plt.scatter(t, cvxvy, label="cvxvy")
#plt.scatter(t, cvxvy, label="cvxvy")


plt.xlim([dt/2, t[-1] * 2])
plt.xscale('log')
#plt.yscale('log')

plt.legend()
plt.show()
