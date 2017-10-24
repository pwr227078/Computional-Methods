from math import cos, sin, pi, sqrt, asin
from quadrat import simpson
import matplotlib.pyplot as pyplot


def real_amp(x):
    r = sqrt(d * d + (y - x) ** 2)
    return ampl * cos(k * r) / sqrt(r)


def imag_amp(x):
    r = (d * d + (y - x) ** 2) ** (1. / 2.)
    return ampl * sin(k * r) / sqrt(r)


ndiv = 100
ampl = 1.
lmbd = 1.
k = 2. * pi / lmbd
d = 50.
a = 2.
yl = -50.
yr = 50.
ny = 100

hy = (yr - yl) / ny

yarr = []
iarr = []

for iy in range(0, ny - 1):
    y = yl + iy * hy
    i = simpson(-a / 2., a / 2., ndiv, real_amp) ** 2 + simpson(-a / 2, a / 2, ndiv, imag_amp) ** 2
    yarr.append(180 / pi * asin(y / sqrt(y ** 2 + d ** 2)))
    iarr.append(i)

pyplot.figure('Exercise 1 - intensity chart', figsize=(9.5, 7))
pyplot.title('Intensity chart\nParameters: λ = %.2f, a = %.2f, d = %.2f' % (lmbd, a, d))
pyplot.xlabel("Angle from the system's axis of symmetry")
pyplot.ylabel('Intensity')
pyplot.plot(yarr, iarr)
pyplot.show()
