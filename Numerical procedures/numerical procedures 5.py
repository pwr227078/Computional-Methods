from math import sqrt


def simpson(xl, xr, ndiv, fun):
    h = (xr - xl) / (ndiv * 2.)
    integral = 0.
    for i in range(0, ndiv):
        integral = integral + 2. * h *\
                                (fun(xl + i * 2. * h) + 4. *
                                    fun(xl + i * 2. * h + h) + fun(xl + i * 2. * h + 2 * h)) / 6.
    return integral


def subintegralFunction(x):
    return sqrt(1 / (1 - x ** 2))


r = 1
arcLength = simpson(0., 0.999999999, 10000000, subintegralFunction)
pi = arcLength * 2
print(pi)
