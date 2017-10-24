from math import log10, sin
import matplotlib.pyplot as pyplot


def fp3(x, h, fun):
    return (fun(x + h) - fun(x - h)) / (2. * h)


def fpp3(x, h, fun):
    return (fun(x + h) + -fun(x - h) - 2. * fun(x)) / h ** 2


arguments = []
firstDerivative = []
secondDerivative = []

x0 = 1
h = 1.
step = 10
while h > 10**(-30):
    arguments.append(-log10(h))
    firstDerivative.append(fp3(x0, h, sin))
    secondDerivative.append(fpp3(x0, h, sin))
    h /= step

pyplot.plot(arguments, firstDerivative, 'b.', arguments, secondDerivative, 'r.')
pyplot.title('Chart of first and second derivatives\nwith respect to the grid parameter h for function sin(%.2f)' % x0)
pyplot.xlabel('-log10(h)')
pyplot.legend(["sin'(%.2f)" % x0, "sin''(%.2f)" % x0])
pyplot.show()

file = open('numerical procedures 2 - output.txt', 'w')
for i in range(0, len(arguments)):
    file.write("%f\t%f\t%f\n" % (arguments[i], firstDerivative[i], secondDerivative[i]))
file.close()
