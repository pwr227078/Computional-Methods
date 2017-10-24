from math import sin, pi, log10
import matplotlib.pyplot as pyplot


def simpson(xl, xr, ndiv, fun):
    h = (xr - xl) / (ndiv * 2.)
    integral = 0.
    for i in range(0, ndiv):
        integral = integral + 2. * h *\
                                (fun(xl + i * 2. * h) + 4. *
                                    fun(xl + i * 2. * h + h) + fun(xl + i * 2. * h + 2 * h)) / 6.
    return integral


xl = 0
xr = pi
func = sin
ndiv = 6
arguments = []
integralValues = []
for j in range(1, 100):
    arguments.append(log10(j))
    integralValues.append(simpson(xl, xr, j, func))
pyplot.plot(arguments, integralValues, '.')
pyplot.title('Convergence of the Simpson algorithm\nwith respect to the grid parameter h')
pyplot.xlabel('log10(h)')
pyplot.ylabel('Integral value')
pyplot.show()

#file = open('numerical procedures 4 - output.txt', 'w')
#for i in range(0, len(integralValues)):
#    file.write('%f\t%f\n' % (i + 1, integralValues[i]))
#file.close()
