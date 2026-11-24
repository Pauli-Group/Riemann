"""Reproduce the three-resonance cusp pulse with python-flint 0.9.0 (Arb balls)."""

from math import prod

from flint import acb, arb, ctx


def remaining_energy(coefficients, exponents, time):
    """Integrate |sum a_j exp(b_j x)|^2 over x <= -time."""
    total = acb(0)
    for a, b in zip(coefficients, exponents):
        for c, d in zip(coefficients, exponents):
            rate = b + d.conjugate()
            total += a * c.conjugate() * (-time * rate).exp() / rate
    assert total.imag.contains(0)
    return total.real


if __name__ == "__main__":
    ctx.prec = 180
    zeros = [acb.zeta_zero(index) for index in (222, 223, 224)]
    assert all(zero.real == arb(1) / 2 for zero in zeros)
    heights = [zero.imag for zero in zeros]
    exponents = [zero / 2 for zero in zeros]
    coefficients = [
        acb(0, height).exp() / prod(height - other for k, other in enumerate(heights) if k != j)
        for j, height in enumerate(heights)
    ]

    # A double zero of the exit amplitude halfway through the observation.
    midpoint_terms = [a * (-2 * b).exp() for a, b in zip(coefficients, exponents)]
    assert sum(midpoint_terms).contains(0)
    assert sum(b * term for b, term in zip(exponents, midpoint_terms)).contains(0)

    retained = remaining_energy(coefficients, exponents, 4) / remaining_energy(coefficients, exponents, 0)
    assert arb("0.85750985011732128") < retained < arb("0.85750985011732130")
    print("Three-mode retained fraction:", retained)
    print("Single-mode retained fraction:", arb(-2).exp())
