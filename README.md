# Class-groups-of-Kummer-extensions
Sage programs used for the article "Class groups of Kummer extensions via cup products in Galois cohomology"

The program [big_modulus_seven.sage](big_modulus_seven.sage) efficiently computes the dimensions h^1_Sigma(F_p(-i)) in the case p = 7, for large sets of prime N = 1 (mod 7). This program allowed us to look at the distribution of these dimensions for N <= 100,000,000.

The program [even_invariant_polynomials.sage](even_invariant_polynomials.sage) was used to compute Kummer generators of certain extensions of Q(mu_p), unramified away from p. This program works for general (regular) p, although we only use it for p = 5 and 7.
