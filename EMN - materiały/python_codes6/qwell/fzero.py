#!/usr/bin/python

# M. P. Polak & P. Scharoch 2017

def bisek(xl, xr, tolx, fun):
#	i=0								#po co to?
	if (fun(xr)*fun(xl)>0.):
		return 'no zero within the specified range'
	while (xr-xl)/2.>tolx:
		x=(xr+xl)/2.
#		i=i+1							#po co to?
		if fun(xr)*fun(x)<0.:
           		xl=x
		else:
			xr=x
	return (xr+xl)/2.
