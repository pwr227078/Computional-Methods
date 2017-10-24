from math import sin
from math import pi


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
integralValue = simpson(xl, xr, 100, func)

file = open('numerical procedures 3 - output.txt', 'w')
file.write('Integral of %s from %f to %f equals: %f' % (func.__name__, xl, xr, integralValue))
file.close()
