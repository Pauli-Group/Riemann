"""Verify the finite cusp-shift reconstruction of m, using no zeta or zero evaluator.

Requires Python 3 and python-flint 0.9.0. All numerical comparisons use Arb balls.
The cusp channel uses Euler-accelerated alternating integer sums. The independent
target uses the positive Gaussian theta pulse and its first logarithmic moment.
"""
from functools import lru_cache
from math import comb

from flint import acb, arb, ctx

ctx.prec = 165
pi = arb.pi()
I = acb(0, 1)
M = 200
X = arb(5)/2
HARMONICS = 6
LOG2 = arb(2).log()
weights = []
tail = 2**M
for j in range(M):
    tail -= comb(M,j)
    weights.append(((-1)**j * arb(tail)/2**M, arb(j+1).log()))


def disk_error(value, error):
    return value + acb(arb(0,error.upper()),arb(0,error.upper()))


@lru_cache(maxsize=None)
def arithmetic_log_derivative(k, energy):
    """zeta'/zeta(w+k), from eta and eta' at Re w=3/2."""
    s = acb(arb(3)/2+k, -energy)
    val, der = acb(0), acb(0)
    for weight, logarithm in weights:
        term = weight*(-s*logarithm).exp()
        val += term
        der -= logarithm*term
    # Euler's kth difference is the Mellin integral of e^-u (1-e^-u)^k.
    # Summing its omitted geometric factors gives 2^-M; differentiation
    # introduces log u and -digamma(s). On u<1 integrate u^(sigma-1)|log u|.
    # On u>=1 use log u<=u. This also accounts for all arithmetic-series error.
    sigma = s.real
    prefactor = arb(1)/(2**M * abs(s.gamma()))
    val_error = prefactor*sigma.gamma()
    der_error = prefactor*(sigma**(-2)+(sigma+1).gamma()) + abs(s.digamma())*val_error
    val = disk_error(val,val_error)
    der = disk_error(der,der_error)
    q = ((1-s)*LOG2).exp()
    return der/val - LOG2*q/(1-q)


def cusp_log_slope(k, energy):
    # d/dw log S0(w/2), with the chain factor 1/2 already included.
    s = acb(arb(3)/2+k,-energy)
    return (1/(s-1)-1/(s+1)+(s/2).digamma()/2-((s+1)/2).digamma()/2
            +arithmetic_log_derivative(k,energy)-arithmetic_log_derivative(k+1,energy))


def endpoint_background(s):
    return 1/s + 1/(s-1) + ((s/2).digamma()-pi.log())/2


def tail_bound(sigma):
    # |zeta'/zeta| <= sum_{n>=2} log n / n^sigma. Sum the first term
    # plus the integral from 2 to infinity (the summand is decreasing here).
    return (2**(-sigma)*LOG2
            + 2**(1-sigma)*(LOG2/(sigma-1)+1/(sigma-1)**2))


def theta_response(energy):
    z = acb(energy,1)
    def integrand(x, derivative):
        a = pi*(2*x).exp()
        theta = acb(0)
        for n in range(1,HARMONICS+1):
            w = a*n*n
            theta += (4*w*w-6*w)*(-w).exp()
        factor = x*(z*x).sin() if derivative else (z*x).cos()
        return (x/2).exp()*theta*factor

    # Bound both omitted half-line integrals using (1+x)e^(3x/2)<=e^(5x/2).
    # With a=pi*n^2, b=a*exp(2x0), each tail from x0 is bounded by
    # 2 a^(-5/4) b^(9/4)e^(-b)/(1-9/(4b)). For n>=7 consecutive
    # first-term bounds have ratio <1/2; for retained harmonics use n<=6.
    a0 = pi*(HARMONICS+1)**2
    remainder = 8*a0*(-a0).exp() + 4*HARMONICS*pi*HARMONICS**2*(
        arb(9)*X/2-pi*(2*X).exp()).exp()
    vals = [disk_error(acb.integral(lambda x, analytic: integrand(x, derivative),
                                    0,X,abs_tol=arb("1e-38"),rel_tol=arb("1e-38")),
                       remainder) for derivative in (False,True)]
    return vals[1]/vals[0]


def main():
    for energy in (10,20,30):
        direct = theta_response(energy)
        if not direct.is_finite():
            raise ArithmeticError("theta integral did not converge")
        print(f"E={energy}  direct m =",direct)
        partial = acb(0)
        for k in range(32):
            partial += cusp_log_slope(k,energy)
            n = k+1
            if n in (8,16,32):
                s = acb(arb(3)/2+n,-energy)
                approx = I*(partial+endpoint_background(s))
                error, bound = abs(approx-direct),tail_bound(s.real)
                if not (error < bound):
                    raise ArithmeticError("cusp reconstruction exceeded its sharp-cutoff bound")
                corrected = approx + I*arithmetic_log_derivative(n,energy)
                if not abs(corrected-direct) < arb("1e-27"):
                    raise ArithmeticError("independent theta/alternating-series cross-check failed")
                print(f"  N={n:2d} error {float(error):.9g}  bound {float(bound):.9g}")


if __name__ == "__main__":
    main()
