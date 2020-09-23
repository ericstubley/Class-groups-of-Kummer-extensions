""" given odd (regular) prime p, and even i, we want the criterion for
when H^1_Sigma(F_p(-i)) /= 0. This is when the H^1_p(F_p(1+i))
extension is split at N, so we want to find a Kummer generator for
that extension and see if its a p-th power at N. The Kummer
generator is the thing in the chi^-i eigenspace of Kummer-ified
Z[zeta_p, 1/p] """

def even_invariant_polynomial(p, i):
    """ returns the minimal polynomial of a Kummer generator of the
    H^1_p(F_p(1+i)) extension when i is odd """
    # set up the fields and relevant group of p-units
    F = GF(p)
    K.<z> = NumberField(cyclotomic_polynomial(p))
    frak_p = K.factor(p)[0][0]
    A = K.S_unit_group(proof=False, S=[frak_p])

    # select a generator of the Galois group of Q(mu_p)/Q which
    # acts on the primitive root of unity z by sending it to its rth
    # power, where r is a primitive root in F_p
    r = primitive_root(p)
    f = 0
    for j in range(len(K.automorphisms())):
        f = K.automorphisms()[j]
        if f(z) == z**r:
            break

    # set up a basis for the unit group, and find the matrix of f in
    # that basis
    gens = A.gens()

    ls = []
    for x in A.gens():
        ls.append(A(f(x)).exponents())

    M = matrix(F, (p-1)/2 + 1, ls).transpose()
    
    # select the eigenvector of M which has eigenvalue chi^-i(r)
    eigen_want = F(r)**(-i)
    vec = []
    for x in M.eigenvectors_right():
        if x[0] == eigen_want:
            vec = x[1][0]
            break

    # translate the eigenvector for M into an element of the
    # cyclotomic field, and return its minimal polynomial
    x = 1
    for j in range(len(vec)):
        x *= K(gens[j])**int(vec[j])

    return x.minpoly()

if __name__ == "__main__":
    p = 7
    for i in range(1, (p-1)/2):
        print even_invariant_polynomial(p, 2*i)
