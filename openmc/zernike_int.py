# Assumes python3
from mpmath import mp
import math

def radial_indefinite(a, n, m):
    # Evaluates the radial indefinite integral at a for n, m
    m = abs(m)

    # Checks if it is a valid point to evaluate at, or if a == 0
    if (n-m) % 2 > 0 or m > n or a == 0:
        return 0
    
    # Due to some weird issues with mpmath, if simultaneously n and m
    # are zero, the value of hyp3f2 is wrong.  Good news is that both
    # n and m are zero, then the solution is trivial.

    hyp = 1
    if m != 0 or n != 0:
        hyp = mp.hyp3f2(-1-n/2, -m/2-n/2, m/2-n/2, -n, -n/2, a**(-2))

    binom = mp.binomial(n, 1/2*(-m+n))

    return float(a**(2+n)/(2+n) * binom * hyp)

def polar_indefinite(a, m):
    result = 0
    absm = abs(m)
    if m > 0:
        result = math.sin(a*absm)/absm
    elif m < 0:
        result = -math.cos(a*absm)/absm
    else:
        result = a
    return result

def integrate_wedge(n, m, a, b, c, d):
    # Integrates Z_n^m over rho from a to b, and over phi from c to d.
    # c, d assumed in radians
    rho_a = radial_indefinite(a,n,m)
    rho_b = radial_indefinite(b,n,m)
    phi_c = polar_indefinite(c,m)
    phi_d = polar_indefinite(d,m)

    return (rho_b - rho_a)*(phi_d - phi_c)

def test():
    import time
    import math

    t1 = time.time()
    answer = [integrate_wedge(31, 1, 0.2, 0.7, 0, math.pi/6) for i in range(100000)]
    t2 = time.time()

    print(t2-t1)

