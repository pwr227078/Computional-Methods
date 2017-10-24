def fp3(x, h, fun):
    return (fun(x + h) - fun(x - h)) / (2. * h)


def fpp3(x, h, fun):
    return (fun(x + h) + fun(x - h) - 2. * fun(x)) / h ** 2


def polynomial(x):
    return x ** 3 + 3 * x ** 2 - 22


firstDerivative = fp3(1.5, 0.1, polynomial)
secondDerivative = fpp3(1., 0.1, polynomial)

file = open('numerical procedures 1 - output.txt', 'w')
file.write('First derivative value: %f \n' % firstDerivative)
file.write('Second derivative value: %f' % secondDerivative)
file.close()

