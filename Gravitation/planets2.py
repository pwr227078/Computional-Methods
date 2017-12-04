import matplotlib.pyplot as pyplot


def F(x1, y1, forces):
    r1q = (x1 * x1 + y1 * y1) ** (1.5)
    forces[0] = -G * M * m1 * x1 / r1q
    forces[1] = -G * M * m1 * y1 / r1q
    return forces


def EL(x1, y1, v1x, v1y, enemom):
    r1 = (x1 * x1 + y1 * y1) ** (1. / 2.)
    v1s = v1x * v1x + v1y * v1y
    enemom[0] = m1 * v1s / 2. - G * m1 * M / r1
    enemom[1] = (x1 * v1y * m1 - y1 * v1x * m1)
    return enemom


G = 1.
M = 100.

m1 = 1.
x1 = [0.0] * 3
y1 = [0.0] * 3
#x1[0] = 1. #elipsa
#y1[0] = 0. #elipsa
#v1x = 3. #elipsa
#v1y = 7. #elipsa
x1[0] = -5. #hiperbola/parabola
y1[0] = 15. #hiperbola/parabola
v1x = 4. #hiperbola/parabola
v1y = -6. #hiperbola/parabola

#tE = 40. #elipsa
tE = 5. #parabola/parabola
Np = 100000

forces = [0.0] * 2
enemom = [0.0] * 2

dt = (tE - 0.0) / Np
dt2 = dt * dt

tarr = []
x1arr = []
y1arr = []
x2arr = []
y2arr = []
enearr = []
momarr = []
F(x1[0], y1[0], forces)
EL(x1[0], y1[0], v1x, v1y, enemom)

x1[1] = x1[0] + v1x * dt + 0.5 * (forces[0] / m1) * dt2
y1[1] = y1[0] + v1y * dt + 0.5 * (forces[1] / m1) * dt2

tarr.append(0.)
x1arr.append(x1[1])
y1arr.append(y1[1])
enearr.append(enemom[0])
momarr.append(enemom[1])

for n in range(1, Np - 1):
    t = 0.0 + n * dt

    F(x1[1], y1[1], forces)

    x1[2] = 2.0 * x1[1] - x1[0] + dt2 * forces[0] / m1
    y1[2] = 2.0 * y1[1] - y1[0] + dt2 * forces[1] / m1

    v1x = (x1[2] - x1[0]) / 2.0 / dt
    v1y = (y1[2] - y1[0]) / 2.0 / dt

    EL(x1[2], y1[2], v1x, v1y, enemom)

    tarr.append(t)
    x1arr.append(x1[1])
    y1arr.append(y1[1])

    enearr.append(enemom[0])
    momarr.append(enemom[1])
    x1[0] = x1[1]
    y1[0] = y1[1]

    x1[1] = x1[2]
    y1[1] = y1[2]

pyplot.plot(x1arr, y1arr, [0], [0], 'r.')
pyplot.xlabel('X')
pyplot.ylabel('Y')
pyplot.title('Trajektoria planety w ruchu wokół gwiazdy')
pyplot.legend(['Trajektoria planety', 'Gwiazda'])
pyplot.show()
