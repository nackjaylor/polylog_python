'''
Python Polylogs

A set of functions for calculating polylogarithms in Python.

Author: Jack Naylor, ACFR, 2023

Relations taken from: https://www.reed.edu/physics/faculty/crandall/papers/Polylog.pdf

'''

import numpy as np
from scipy.special import zeta, factorial, bernoulli, binom

# https://www.reed.edu/physics/faculty/crandall/papers/Polylog.pdf
def f_n_0(n, z, L):

    part = 1.0
    sum_val = 0.0
    k = 1
    zk = z
    part = zk/(k**n)
    while k - 1 < L and part > 1e-16:
        part = zk/(k**n)

        sum_val += part
        zk *= z
        k += 1
    return sum_val

def h_sum(q):

    if q == 0:
        return 0.0
    else:
        part_sum = 0.0
        for i in range(1, int(q)+1):
            part_sum += 1/i
        return part_sum

def f_n_1(n, z, L):

    part = 1.0
    sum_val = 0.0
    m = 0
    logmz = 1
    while m < L and abs(part) > 1e-16:

        if (n-m) != 1:
            part = zeta(n - m)/factorial(m)*logmz
            sum_val += part

        logmz *= np.log(z)
        m += 1
    
    return sum_val +  np.log(z)**(n-1)/factorial(n-1)*(h_sum(n-1)-np.log(-np.log(z)))


def f_n_n1(n, z, L):

    part_sum = 1.0
    sum_val = 0.0
    i = 0
    logiz = 1
    while i < L and abs(part_sum) > 1e-16:
        part_sum = bernoulli(i-n+1)[-1]/(factorial(i)*(i-n+1))*logiz
        sum_val += part_sum
        logiz *= np.log(z)
        i+=1
    
    return factorial(-n)*(-np.log(z))**(n-1)-sum_val


def bernoulli_poly(n, x):
    coeffs = bernoulli(n)
    sum_val = 0.0
    for k in range(0, n+1):
        sum_val += binom(n, k) * coeffs(n-k) * x**k
    return sum_val

def phi_func(z):
    if np.real(z) < 0.0 or z > 1.0: return 1
    else: return 0

def g_n(n,z):

    sec_term = 0.0
    if phi_func(z):
        sec_term = -2*np.pi*1j*phi_func(z)*np.log(z)**(n-1)/factorial(n-1)
    return -(2*np.pi*1j)**n/factorial(n)*bernoulli_poly(n, np.log(z)/(2*np.pi*1j)) + sec_term


def polylog(n, z, prec = 16):
    if np.imag(z) == 0:
        return polylog_real(n, z, prec = prec)
    sum_limit = np.ceil(prec*np.log2(10))

    if z == 1:
        return zeta(n)
        
    elif z == -1:
        return -(1-2**(1-n))*zeta(n)
    elif n == 0:
        return z/(1-z+np.eps)
    elif n == 1:
        return -np.log(1-z+np.eps)
    elif n == -1:
        return z/(1-z)**2
    
    elif np.abs(z) <= 1/2:
        return f_n_0(n, z, sum_limit)
    elif np.abs(z) >= 2:
        return g_n(n,z) - (-1)**n*f_n_0(n, 1/z, sum_limit)
    
    if np.sign(n) == 1:
        return f_n_1(n, z, sum_limit)
    else:
        return f_n_n1(n, z, sum_limit)

def polylog_real(n, z, prec = 16):

    sum_limit = np.ceil(prec*np.log2(10))

    if z == 1:
        return zeta(n)
        
    elif z == -1:
        return -(1-2**(1-n))*zeta(n)
    
    elif np.abs(z) < 1/4:
        return f_n_0(n, z, sum_limit)
    
    elif z < 0:
        return 2**(1-n)*polylog_real(n, z**2)-polylog_real(n, -z)
    
    return f_n_1(n, z, sum_limit)


def polylog_rec(n, z, prec = 10):

    if np.abs(z) < 2: return polylog(n, z, prec)
    return 2**(n - 1)*(polylog_rec(n, np.sqrt(z))+polylog_rec(n, -np.sqrt(z)))



def main():

    print(np.pi**2/12 - 1/2*np.log(2)**2)
    print(polylog_rec(2, 1/2))
    print(f_n_0(2, 1/2, 1000))


if __name__ == "__main__":
    main()
