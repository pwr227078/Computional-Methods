#!/usr/bin/python

import matplotlib.pyplot as pyplot


def euler(a, b, n, y0, fun):
    h = (b - a) / n
    yp = [y0]
    xp = [a]
    for i in range(0, n):
        yp.append(yp[i]+h*fun(xp[i],yp[i]))
        xp.append(xp[i]+h)
    return xp, yp


def rk4(a, b, n, y0, fun):
    h = (b - a) / n
    yp = [y0]
    xp = [a]
    for i in range(0, n):
        k1 = h * fun(xp[i], yp[i])
        k2 = h * fun(xp[i] + h / 2., yp[i] + k1 / 2.)
        k3 = h * fun(xp[i] + h / 2., yp[i] + k2 / 2.)
        k4 = h * fun(xp[i] + h, yp[i] + k3)
        yp.append(yp[i] + (k1 + 2. * k2 + 2. * k3 + k4) / 6.)
        xp.append(xp[i]+h)
    return xp, yp


def ab3(a, b, n, y0, fun):
    h = (b - a) / n
    yp = [y0]
    xp = [a]
    for i in range(0,2):
        k1 = h * fun(xp[i], yp[i])
        k2 = h * fun(xp[i] + h / 2., yp[i] + k1 / 2.)
        k3 = h * fun(xp[i] + h / 2., yp[i] + k2 / 2.)
        k4 = h * fun(xp[i] + h, yp[i] + k3)
        yp.append(yp[i] + (k1 + 2. * k2 + 2. * k3 + k4) / 6.)
        xp.append(xp[i]+h)
    for i in range(2, n):
        yp.append(yp[i] + h / 12. * (23. * fun(xp[i], yp[i]) - 16. * fun(xp[i - 1], yp[i - 1]) +
                                     5. * fun(xp[i - 2], yp[i - 2])))
        xp.append(xp[i]+h)
    return xp, yp


def func(x, y):
    return -y


e1, e2 = euler(0, 1, 1, 1, func)
r1, r2 = rk4(0, 1, 1, 1, func)
a1, a2 = ab3(0, 1, 1, 1, func)
pyplot.plot(e1, e2, 'r', r1, r2, 'b', a1, a2, 'g')
pyplot.legend(['Euler', 'Runge - Kutta', 'Adams - Bashforth'])
pyplot.show()

file = open('euler_1.txt', 'w')
file.write(str(euler(0, 10, 100, 1, func)))
file.close()

file = open('rk4_1.txt', 'w')
file.write(str(rk4(0, 10, 100, 1, func)))
file.close()

file = open('ab3_1.txt', 'w')
file.write(str(ab3(0, 10, 100, 1, func)))
file.close()
