import matplotlib.pyplot as pyplot


def euler(a, b, n, fun1, y0, fun2, v0):
    h = (b - a) / n
    yp1 = [y0]
    v = [v0]
    xp = [a]
    for i in range(0, n):
        yp1.append(yp1[i] + h * fun2(xp[i], yp1[i], v[i]))
        v.append(v[i] + h * fun1(xp[i], yp1[i], v[i]))
        xp.append(xp[i]+h)
    return xp, yp1, v

def der1(x, y, v):
    return -k * y

def der2(x, y, v):
    return v


k = 1
x, y, v = euler(0, 20, 10000, der1, 0, der2, 1)
pyplot.plot(x, y)
pyplot.xlabel('x')
pyplot.ylabel('y')
pyplot.show()
