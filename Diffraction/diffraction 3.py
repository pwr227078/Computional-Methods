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
lmbd = 3.5
k = 2. * pi / lmbd
d = 30.
yl = -200.
yr = 200.
ny = 400
hy = (yr - yl) / ny
numberOfSlits = 2
slitWidth = 10.
distanceBetweenSlits = 8
yarr = []
iarr = []
xrMAX = 0.5 * (numberOfSlits * slitWidth + (numberOfSlits - 1) * distanceBetweenSlits)

for iy in range(0, ny - 1):
    y = yl + iy * hy
    intensity = 0
    xl = -0.5 * (numberOfSlits * slitWidth + (numberOfSlits - 1) * distanceBetweenSlits)
    xr = xl + slitWidth
    while xr <= xrMAX:
        intensity += simpson(xl, xr, ndiv, real_amp) ** 2 + simpson(xl, xr, ndiv, imag_amp) ** 2
        xl += slitWidth + distanceBetweenSlits
        xr = xl + slitWidth
    yarr.append(180 / pi * asin(y / sqrt(y ** 2 + d ** 2)))
    iarr.append(intensity)

pyplot.figure('Multi-slit chart', figsize=(9.5, 7))
pyplot.title('Intensity of wave\nParameters: λ = %.2f, a = %.2f, number of slits= %.2f, distance between slits = %.2f'
             % (lmbd, slitWidth, numberOfSlits, distanceBetweenSlits))
pyplot.xlabel("Angle from the system's axis of symmetry")
pyplot.ylabel('Intensity')
pyplot.plot(yarr, iarr)
pyplot.show()
