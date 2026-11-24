# Hilbert–Pólya approach to the Riemann theorem

Pierre-Luc Dallaire-Demers · Pauli Group · November 2026

[Read the paper (PDF)](hilbert-polya-riemann-quantum-physicists.pdf) · [LaTeX source](hilbert-polya-riemann-quantum-physicists.tex)

How can prime numbers specify a quantum Hamiltonian whose energies are the heights of the nontrivial zeros of the Riemann zeta function? Taking the critical-line conclusion of the Riemann hypothesis as input, the paper constructs a self-adjoint Hamiltonian from the thermal fluctuations of a gas with one mode per prime. Those fluctuations fix the couplings of a nearest-neighbor spin chain. Measuring a frequency at the chain’s end gives a zero height after a change of variables, and the strength of that spectral line gives its multiplicity.

For an open wave on the modular surface, the real part of a zeta zero sets the decay rate of a resonance and its imaginary part sets the oscillation frequency. The paper uses a wave from Alpöge’s theta lift to prepare a scattering experiment on that surface. The reflected wave recovers the spin chain’s boundary response. After compensating for the common decay, a weighted trace over the resonances reproduces the smoothed prime count. A separate measurement on the shear phase circle gives a transmission-line realization with a specified electrical termination.

## Build the paper

From the repository root, run the following command with Tectonic. The vector figures are included in [figures/manuscript](figures/manuscript).

```sh
tectonic --reruns 2 hilbert-polya-riemann-quantum-physicists.tex
```

The [cusp-response check](checks/check_cusp_to_prime_response.py) reproduces the three damped-frequency comparisons in the modular scattering section. Run it with Python 3 and `python-flint==0.9.0`.

The [outgoing pulse calculation](checks/check_cusp_retention.py) uses the same dependencies to reproduce the energy retained by a superposition of three modular resonances and by each resonance in isolation.
