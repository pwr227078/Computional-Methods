from math import cos, sin, pi, sqrt, asin
from quadrat import simpson
import matplotlib.pyplot as pyplot


def real_amp(x):
    r = sqrt(d * d + (y - x) ** 2)
    return ampl * cos(k * r) / sqrt(r)


def imag_amp(x):
    r = (d * d + (y - x) ** 2) ** (1. / 2.)
    return ampl * sin(k * r) / sqrt(r)


ndiv = 1000
ampl = 1.
yl = -50.
yr = 50.
ny = 100
hy = (yr - yl) / ny
yarr = []
iarr = []
lmbd = 2.
k = 2. * pi / lmbd
pyplot.figure('Exercise 2 - intensity(a) chart', figsize=(9.5, 7))
pyplot.title('Intensity chart for different $a$ values\nParameters: λ = %.2f, d = 50' % lmbd)
pyplot.xlabel("Angle from the system's axis of symmetry")
pyplot.ylabel('Intensity')
d = 50.
a = 7.
for iy in range(0, ny - 1):
    y = yl + iy * hy
    i = simpson(-a / 2., a / 2., ndiv, real_amp) ** 2 + simpson(-a / 2, a / 2, ndiv, imag_amp) ** 2
    yarr.append(180 / pi * asin(y / sqrt(y ** 2 + d ** 2)))
    iarr.append(i)
pyplot.plot(yarr, iarr, 'r')
iarr.clear()

a = 10.
for iy in range(0, ny - 1):
    y = yl + iy * hy
    i = simpson(-a / 2., a / 2., ndiv, real_amp) ** 2 + simpson(-a / 2, a / 2, ndiv, imag_amp) ** 2
    iarr.append(i)
pyplot.plot(yarr, iarr, 'k')
iarr.clear()

a = 12.
for iy in range(0, ny - 1):
    y = yl + iy * hy
    i = simpson(-a / 2., a / 2., ndiv, real_amp) ** 2 + simpson(-a / 2, a / 2, ndiv, imag_amp) ** 2
    iarr.append(i)
pyplot.plot(yarr, iarr, 'b')
pyplot.legend(['a = 7', 'a = 10', 'a = 12'])

pyplot.figure('Exercise 2 - intensity(d) chart', figsize=(9.5, 7))
pyplot.title('Intensity chart for different $d$ values\nParameters: λ = %.2f, a = 12' % lmbd)
pyplot.xlabel("Angle from the system's axis of symmetry")
pyplot.ylabel('Intensity')
pyplot.plot(yarr, iarr, 'b')
yarr.clear()
iarr.clear()

d = 100.
for iy in range(0, ny - 1):
    y = yl + iy * hy
    i = simpson(-a / 2., a / 2., ndiv, real_amp) ** 2 + simpson(-a / 2, a / 2, ndiv, imag_amp) ** 2
    yarr.append(180 / pi * asin(y / sqrt(y ** 2 + d ** 2)))
    iarr.append(i)
pyplot.plot(yarr, iarr, 'y')
yarr.clear()
iarr.clear()

d = 150.
for iy in range(0, ny - 1):
    y = yl + iy * hy
    i = simpson(-a / 2., a / 2., ndiv, real_amp) ** 2 + simpson(-a / 2, a / 2, ndiv, imag_amp) ** 2
    yarr.append(180 / pi * asin(y / sqrt(y ** 2 + d ** 2)))
    iarr.append(i)
pyplot.plot(yarr, iarr, 'g')
iarr.clear()

pyplot.legend(['d = 50', 'd = 100', 'd = 150'])
pyplot.show()
