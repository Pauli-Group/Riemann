# Hilbert–Pólya approach to the Riemann theorem

Pierre-Luc Dallaire-Demers · Pauli Group · November 2026

[Read the paper (PDF)](hilbert-polya-riemann-quantum-physicists.pdf) · [LaTeX source](hilbert-polya-riemann-quantum-physicists.tex)

How can prime numbers specify a quantum Hamiltonian whose energies are the heights of the nontrivial zeros of the Riemann zeta function? Taking the critical-line conclusion of the Riemann hypothesis as input, the paper constructs a self-adjoint Hamiltonian from the thermal fluctuations of a gas with one mode per prime. Those fluctuations fix the couplings of a nearest-neighbor spin chain. Measuring a frequency at the chain’s end gives a zero height after a change of variables, and the strength of that spectral line gives its multiplicity.

When the reconstructed modes evolve at their signed zeta energies, the Ramsey interference signal changes slope at the logarithms of prime powers. The paper also studies a wave in Alpöge’s theta lift. The lift supplies a measurable phase. Counting each positive harmonic once gives the boundary response of a lossless transmission line. A separate static matching condition selects the short circuit whose standing waves have the zeta frequencies.

## Build the paper

From the repository root, run the following command with Tectonic. The eight vector figures are included in [figures/manuscript](figures/manuscript).

```sh
tectonic --reruns 2 hilbert-polya-riemann-quantum-physicists.tex
```
