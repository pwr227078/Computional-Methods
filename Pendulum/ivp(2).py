#!/usr/bin/python

from math import log10
import matplotlib.pyplot as pyplot


def euler(a, b, n, y0, fun):
    h = (b - a) / n
    yp = [y0]
    xp = [a]
    for i in range(0, n):
        yp.append(yp[i] + h * fun(xp[i], yp[i]))
        xp.append(xp[i] + h)
    return -log10(h), yp[len(yp)-1]


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
        xp.append(xp[i] + h)
    return yp[len(yp) - 1]


def ab3(a, b, n, y0, fun):
    h = (b - a) / n
    yp = [y0]
    xp = [a]
    for i in range(0, 2):
        k1 = h * fun(xp[i], yp[i])
        k2 = h * fun(xp[i] + h / 2., yp[i] + k1 / 2.)
        k3 = h * fun(xp[i] + h / 2., yp[i] + k2 / 2.)
        k4 = h * fun(xp[i] + h, yp[i] + k3)
        yp.append(yp[i] + (k1 + 2. * k2 + 2. * k3 + k4) / 6.)
        xp.append(xp[i] + h)
    for i in range(2, n):
        yp.append(yp[i] + h / 12. * (23. * fun(xp[i], yp[i]) - 16. * fun(xp[i - 1], yp[i - 1])
                                     + 5. * fun(xp[i - 2], yp[i - 2])))
        xp.append(xp[i] + h)
    return yp[len(yp) - 1]


def func(x, y):
    return -y


x = []
y1 = []
y2 = []
y3 = []
i = 1
while i <= 10 ** 7:
    tempX, tempY = euler(0, 1, i, 1, func)
    x.append(tempX)
    y1.append(tempY)
    y2.append(rk4(0, 1, i, 1, func))
    y3.append(ab3(0, 1, i, 1, func))
    i *= 10
pyplot.plot(x, y1, 'r.', x, y2, 'b.', x, y3, 'y.')
pyplot.xlabel('-log10(h)')
pyplot.ylabel('y')
pyplot.title('Badanie zbieżności algorytmów')
pyplot.legend(['Euler', 'Runge - Kutta 4 punktowy', 'Adams - Bashforth 3 punktowy'])
pyplot.show()
