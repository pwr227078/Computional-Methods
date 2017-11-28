#!/usr/bin/python

from math import sin, cos, sqrt, pi
from ivp2d import rk4_2d
from matplotlib import pyplot


def dftheta(t, theta, omega):
    return omega


def dfomega(t, theta, omega):
    return -(g / l) * sin(theta)


def T1(eps):
    return (pi * sqrt(2 * l / eps / g))


g = 9.81
frEtop = .5
tE = 10.
np = 1000

l = 1.
m = 1.
E0 = frEtop * m * 2. * l * g
t0 = 0.
theta0 = 0.
omega0 = sqrt(2. * E0 / m) / l

h = (tE - t0) / np
eps = []
T = []
Tc = []


t, theta, omega = rk4_2d(t0, tE, np, theta0, omega0, dftheta, dfomega)
Etot = [0.5 * m * l ** 2 * om ** 2 +
        m * g * l * (1. - cos(thet)) for om, thet in zip(omega, theta)]

i = 0
while i < np:
    eps.append(frEtop)
    if(frEtop < 1.0):
        T.append(2 * pi / omega[i])
    else:
        T.append(T1(frEtop))
    frEtop += 0.01
    Tc.append(T[i] / (2 * pi * sqrt(l / g)) - 1)
    i += 1


pyplot.figure(1)
pyplot.plot(eps, T)
pyplot.title('Badanie przebiegu T(eps)')
pyplot.xlabel('eps')
pyplot.ylabel('T')

pyplot.figure(2)
pyplot.plot(eps, Tc)
pyplot.title('(T/T0)-1')
pyplot.xlabel('eps')
pyplot.ylabel('T')
pyplot.show()
