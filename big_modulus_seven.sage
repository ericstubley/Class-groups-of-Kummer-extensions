""" Computes the dimensions h^1_Sigma(F_p(-i)) when p = 7,
i = 1, 2, 3, 4. In order to facilitate computing these numbers
for large N (we've computed N <= 100,000,000) this program 
treats the factorials needed to compute the M_i as numbers mod
th product of all N considered rather than in characteristic 0.
The main method is seven_test, which takes the upper bound on
N and produces two text files as output: one with the data of the
dimensions h^1_Sigma(F_p(-i)), one with timings of how long
the program takes to compute every 100,000 primes. Note that this
program does not actually compute ranks of class groups; the data
file output has a column for the rank, which is 0 unless the rank
is determined by the dimensions h^1_Sigma(F_p(-i)). """


p = 7

import time

def m1_m3_calc(N, fac_a, fac_b, fac_c):
    """ Given a prime N = 1 mod 7 and computes the dimensions
    h^1_Sigma(F_p(-1)) and h^1_Sigma(F_p(-3)). We check this by
    checking whether or not S_1 and S_3 are 7th powers mod N. We
    start with the factorials as numbers mod the large modulus
    computed at the beginning of the program, and reduce them mod N
    before proceeding with the computation. Note that this don't
    actually compute the S_i, but instead computes m1 and m3 which
    are 7th powers iff the S_1 and S_3 are; these mi are used to
    reduce the number of operations necessary. """
    red_fac_a = mod(fac_a, N) 
    red_fac_b = mod(fac_b, N) 
    red_fac_c = mod(fac_c, N) 

    red_fac_c3 = red_fac_c^3
    red_fac_abc3 = red_fac_a*red_fac_b*red_fac_c3

    m1 = red_fac_abc3*red_fac_b
    m3 = red_fac_abc3*red_fac_c3

    M = (N-1).divide_knowing_divisible_by(p)
    inv1 = 1 if (m1^M == 1) else 0
    inv3 = 1 if (m3^M == 1) else 0

    return inv1, inv3


def m2_calc(N):
    """ Test if the roots of f(x) = x^3 + 41x^2 + 54x + 1 are pth
    powers mod N. Something is a pth power mod N if and only if its
    (N-1)/pth power is 1, so we construct the matrix A with charpoly
    f and test whether or not 1 is an eigenvalue of A^((N-1)/p).
    Returns 1 if the roots are pth powers mod N, 0 else. """
    A = Matrix([[0,0,-1],[1,0,-54],[0,1,-41]]).change_ring(Integers(N))
    M = (N-1).divide_knowing_divisible_by(p)
    B = A^M
    charpoly_at_1 = hard_charpoly_substitute(B)
    if charpoly_at_1 == 0:
        return 1
    else:
        return 0

def m4_calc(N):
    """ Test if the roots of f(x) = x^3 - 25x^2 + 31x + 1 are pth
    powers mod N. Something is a pth power mod N if and only if its
    (N-1)/pth power is 1, so we construct the matrix A with charpoly
    f and test whether or not 1 is an eigenvalue of A^((N-1)/p).
    Returns 1 if the roots are pth powers mod N, 0 else. """
    A = Matrix([[0,0,-1],[1,0,-31],[0,1,25]]).change_ring(Integers(N))
    M = (N-1).divide_knowing_divisible_by(p)
    B = A^M
    charpoly_at_1 = hard_charpoly_substitute(B)
    if charpoly_at_1 == 0:
        return 1
    else:
        return 0

def hard_charpoly_substitute(M):
    """ Hardcoded function which evaluates the characteristic
    polynomial of a 3x3 matrix at 1. This is hardcoded rather than
    using Sage's charpoly functions to get around some very strange
    bugs with that functionality that were causing crashes. """
    t1 = (1 - M[0][0])*((1 - M[1][1])*(1 - M[2][2]) - M[1][2]*M[2][1])
    t2 = (-M[1][0])*(-M[0][1]*(1 - M[2][2]) - M[0][2]*M[2][1])
    t3 = (-M[2][0])*(M[0][1]*M[1][2] + M[0][2]*(1 - M[1][1]))

    return t1 - t2 + t3


def big_modulus(limit):
    """ Function which computes the product all of primes N = 1 
    mod 7 which are less than limit. For large values of limit this
    can take a while, so the function prints some messages
    informing the user of where the computation is at. """
    print "computing the big modulus"
    ls = []
    count = 1
    for N in primes(limit):
        if p.divides(N-1):
            ls.append(N)
        if N > count*100000:
            print "done up to {}".format(count*100000)
            count += 1

    print "computing the big product"
    bm = prod(N for N in ls)
    return bm

def seven_test(limit):
    """ This is the main function. After setting up the large 
    modulus used in computation, a loop runs over all primes N = 1
    mod 7. For each the dimensions h^1_Sigma(F_p(-i)) are calculated
    and the data is written to a file. A file is also created with
    information on how long the computation is taking, updating every
    100,000 primes N. """
    # set up the big modulus and factorials and stuff
    bm = big_modulus(limit)
    R = Integers(bm)

    fac_a_bm, fac_b_bm, fac_c_bm = R(1), R(1), R(1)
    count = 1
    prev_N = 1

    start_time = time.time()
    prev_time = start_time

    with open("./seven_data.txt", 'w') as f_data, open("./seven_data_times.txt", 'w') as f_times:
        for N in primes(limit):
            if p.divides(N-1):
                M = (N-1).divide_knowing_divisible_by(p)    

                # update the factorials
                prev_M = (prev_N-1).divide_knowing_divisible_by(p)
                fac_a_bm *= prod(k for k in [prev_M+1..M])
                fac_b_bm *= prod(k for k in [2*prev_M+1..2*M])
                fac_c_bm *= prod(k for k in [3*prev_M+1..3*M])


                # actual computations
                # note that this program does not compute ranks of
                # class groups; the data file output
                h1_1, h1_2, h1_3, h1_4 = 0, 0, 0, 0
                rank = 0

                h1_1, h1_3 = m1_m3_calc(N, fac_a_bm, fac_b_bm, fac_c_bm)

                if h1_1 == 1:
                    h1_4 = m4_calc(N)
                if h1_3 == 1:
                    h1_2 = m2_calc(N)

                dim_str = "{}{}{}{}".format(h1_1, h1_2, h1_3, h1_4)

                if dim_str == "0000":
                    rank = 1
                elif (dim_str == "1000") or (dim_str == "0010"):
                    rank = 2

                # data output
                out_str = "N = {}, dims = {}, rank = {}".format(N, dim_str, rank)
                f_data.write(out_str + '\n')

                # timing output
                # print time info every 100000
                if N > count*100000:
                    curr_time = time.time()
                    print_time = round(curr_time - prev_time, 3)
                     
                    time_str = "{} in {} seconds"
                    time_str = time_str.format(count*100000, print_time)
                    
                    print time_str
                    f_times.write(time_str + '\n')
                    count += 1
                    prev_time = curr_time

                # update the prev_N
                prev_N = N

        # one final time update before closing files
        curr_time = round(time.time() - prev_time, 3) 
        time_str = "{} in {} seconds"
        time_str = time_str.format(count*100000, curr_time)
        
        print time_str
        f_times.write(time_str + '\n')

if __name__ == "__main__":
    seven_test(100000000)
