#!/usr/bin/python

# M. P. Polak & P. Scharoch 2017

from math import sin,cos,sqrt
from ivp2d import rk4_2d
from matplotlib import pyplot


def dftheta(t, theta, omega):
    return omega


def dfomega(t, theta, omega):
    return -(g / l) * sin(theta)


g = 9.81
frEtop = 0.1
tE = 10.
np = 1000

l = 1.
m = 1.
E0 = frEtop * m * 2. * l * g
t0 = 0.
theta0 = 0.
omega0 = sqrt(2. * E0 / m) / l

h = (tE - t0) / np

t, theta, omega = rk4_2d(t0, tE, np, theta0, omega0, dftheta, dfomega)
Etot = [0.5*m*l*l*om**2+m*g*l*(1.-cos(thet)) for om, thet in zip(omega, theta)]

pyplot.figure(1)
pyplot.plot(t, theta)
pyplot.title('For frEtop = %.2f' % frEtop)
pyplot.xlabel('t')
pyplot.ylabel(('theta'))


print(Etot)
pyplot.figure(2)
pyplot.plot(t, Etot)
pyplot.title('For frEtop = %.2f' % frEtop)
pyplot.xlabel('t')
pyplot.ylabel(('Etot'))



frEtop = 1.01
E0 = frEtop * m * 2. * l * g
t0 = 0.
theta0 = 0.
omega0 = sqrt(2. * E0 / m) / l


h = (tE - t0) / np
t, theta, omega = rk4_2d(t0, tE, np, theta0, omega0, dftheta, dfomega)
Etot = [0.5*m*l*l*om**2+m*g*l*(1.-cos(thet)) for om, thet in zip(omega, theta)]
print(Etot)
pyplot.figure(3)
pyplot.plot(t, Etot)
pyplot.title('For frEtop = %.2f' % frEtop)
pyplot.xlabel('t')
pyplot.ylabel(('Etot'))



frEtop = 10.
E0 = frEtop * m * 2. * l * g
t0 = 0.
theta0 = 0.
omega0 = sqrt(2. * E0 / m) / l

h = (tE - t0) / np
t, theta, omega = rk4_2d(t0, tE, np, theta0, omega0, dftheta, dfomega)
Etot = [0.5*m*l*l*om**2+m*g*l*(1.-cos(thet)) for om, thet in zip(omega, theta)]
print(Etot)
pyplot.figure(4)
pyplot.plot(t, Etot)
pyplot.title('For frEtop = %.2f' % frEtop)
pyplot.xlabel('t')
pyplot.ylabel('Etot')
pyplot.show()
