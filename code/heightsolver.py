#----- Author: Nick Konz -----#
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import astropy as ap
import scipy as sp
import matplotlib.pyplot as plt

rC = 8.1 #kpc
v = 220 #km/s
l = np.radians(83)
delta = np.radians(5)
c = 3e5

fN = 1420.41e6 #normal
fR = 1420.5e6 #wavelengths for R, G, B
fG = 1420.71e6
fB = 1420.85e6

def w(f):
    return (3e8)/f

wN = w(fN)
wR = w(fR)
wG = w(fG)
wB = w(fB)

pm = +1;

def wDiff(w):
    return w - wN

vLOSg = c * wDiff(wG) / wN
vLOSb = c*wDiff(wB) / wN

def height(d):
    vLOSg = c * wDiff(wG) / wN
    b = np.arccos(vLOSg/v)
    a = np.pi/2 - b
    rg = rC*np.sin(l)/np.sin(a)
    vLOSb = c*wDiff(wB) / wN
    xi = np.arccos(vLOSb / v)
    et = np.pi/2 - xi
    lm = np.pi - a
    rb = rg*np.sin(lm)/np.sin(et)
    t = np.pi - lm - et
    dgb = np.sqrt(rg**2 + rb**2 - 2*rg*rb*np.cos(t))
    g = np.pi - l - a
    drg = rc*np.sin(g)/np.sin(a)
    return (drg + dgb - drg + np.sqrt(dgb**2 - 4*(drg + dgb)*drg*(np.tan(d)**2))) / 2*np.tan(d)

def circ(x,y):
    return -4 - x + x**2 + y**2

def h(dRG, dGB):
    return (dRG + dGB - dRG + pm*np.sqrt(dGB**2 - 4*(dRG + dGB)*dRG*(np.tan(delta)**2))) / 2*np.tan(delta)

def beta(dRG, dGB):
    return np.arccos(vLOSg/(v*np.cos(np.arctan(h(dRG, dGB)/dRG))))

def alpha(dRG, dGB):
    return np.pi/2 - beta(dRG, dGB)

def gamma(dRG, dGB):
    return np.pi - 1 - alpha(dRG, dGB)

def xi(dRG, dGB):
    return np.arccos(vLOSb/(v*np.cos(np.arctan(h(dRG, dGB)/(dRG+dGB)))))

def eta(dRG, dGB):
    return np.pi/2 - xi(dRG, dGB)

def lambd(dRG, dGB):
    return np.pi - alpha(dRG, dGB)

def rB(dRG, dGB):
    return rC*np.sin(lambd(dRG, dGB))/np.sin(eta(dRG, dGB))

def rG(dRG, dGB):
    return rC*np.sin(l)/np.sin(alpha(dRG, dGB))

def vtheta(dRG, dGB):
    return np.pi - lambd(dRG, dGB) - eta(dRG, dGB)

def f(dRG, dGB): # = 0
    return np.sqrt(rC**2 + rG(dRG, dGB)**2 - 2*rC*rG(dRG, dGB)*np.cos(gamma(dRG, dGB))) - dRG

def g(dRG, dGB): # = 0
    return np.sqrt(rG(dRG, dGB)**2 + rB(dRG, dGB)**2 - 2*rB(dRG, dGB)*rG(dRG, dGB)*np.cos(vtheta(dRG, dGB))) - dGB


def main():

    #parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    #parser.add_argument("N", type=str, help="number of matches generated")
    
    #args = parser.parse_args()
    #N = args.N


    
    #x = np.linspace(0, 0.14, 1000)

    x = np.linspace(0,50,1000)
    y = np.linspace(0,50,1000)
    X, Y = np.meshgrid(x,y)

    F = f(X,Y)
    G = g(X,Y)

    Z = h(X,Y)

    plt.figure(1)

    plt.imshow(Z, extent=[0, 50, 0, 50], origin='lower',
           cmap='viridis', alpha=0.8)
    plt.colorbar(label="height above Milky Way [kPC]")
    plt.axis(aspect='image');


    #for i, val in enumerate(x.tolist()):
    #    print(val, height(val))
    plt.contour(X,Y,F,[0], label="Condition 1")
    plt.contour(X,Y,G,[0], cmap=plt.cm.hot, label="Condition 2")
    plt.title("Fit for Earth's Height Above a Flat Milky Way")
    plt.ylabel('$d_{rg}$ (distance from Earth to green arm [kPC])')
    plt.xlabel('$d_{gb}$ (distance from green arm to blue arm [kPC])')
    plt.legend()
    plt.yscale('log')
    #plt.ylim(5, 0)

    #print(h(28.12, 19.59))

    plt.show()


    
main()