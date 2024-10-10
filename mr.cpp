/*
 * The Miller-Rabin primality test
 *
 * Written by Charles Liu 5/5/2017
 *
 * NOTE:
 *	The code uses modular exponentiation by repeated squaring, the running
 *time of
 *	this algorithm is O(k log^3 n)
 *
 * 	The code has big room for improvements, but it does work as advertised.
 */

#include "mr.h"


using namespace std;
using namespace boost::multiprecision;

/*
 * random number engine
 */
static default_random_engine e(time(0));

/*
 * Fast calculation of `a^x mod n´ by using right-to-left binary modular
 * exponentiation.
 *
 * See http://en.wikipedia.org/wiki/Modular_exponentiation
 */
uint1024_t pow_mod(uint1024_t a, uint1024_t x, uint1024_t n) {
    /*
     * Note that this code is sensitive to overflowing for testing of large
     * prime numbers.  The `a*r´ and `a*a´ operations can overflow.  One easy
     * way of solving this is to use 128-bit precision for calculating a*b % n,
     * since the mod operator should always get us back to 64bits again.
     *
     * You can either use GCC's built-in __int128_t or use
     *
     * typedef unsigned int uint128_t __attribute__((mode(TI)));
     *
     * to create a 128-bit datatype.
     */
    uint1024_t r = 1;
    while (x) {
        if ((x & 1) == 1)
            r = a * r % n;
        a = a * a % n;
        x >>= 1;
    }
    return r;
}

/*
 * The Miller-Rabin probabilistic primality test.
 *
 * Returns true if ``n´´ is PROBABLY prime, false if it's composite.
 * The parameter ``k´´ is the accuracy.
 *
 * The running time should be somewhere around O(k log^3 n).
 *
 */

uint1024_t rrand(uint1024_t a, uint1024_t b) {
    boost::random::mt19937 rng(clock());
    boost::random::uniform_int_distribution<uint1024_t> dist(a, b);
    return dist(rng);
}

bool is_prime(uint1024_t n) {
    static std::vector<uint1024_t> mainPrimes = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997 };

    if (find(mainPrimes.begin(), mainPrimes.end(), n) != mainPrimes.end()) return true;
    for (auto &i : mainPrimes) if (n % i == 0) return false;

    // Must have ODD n greater than THREE
    if (n == 2 || n == 3)
        return true;
    if (n <= 1 || !(n & 1))
        return false;

    // Write n-1 as d*2^s by factoring powers of 2 from n-1
    int s = 0;
    uint1024_t d = n - 1;
    for (; !(d & 1); ++s, d >>= 1)
        ; // loop

    // Here, we CANNOT set it static
    // uniform_int_distribution<uint1024_t> u(2, n - 2);

    bool flg = false;
    for (int i = 0; i < int(log2(cpp_bin_float_100(n))); ++i) { // log2(uint256) << uint256
        uint1024_t a = rrand(2, n - 2);
        uint1024_t reminder = pow_mod(a, d, n);
        if (reminder == 1 | reminder == (n - 1))
            continue;
        for (int r = 1; r <= s - 1; ++r) {
            flg = false;
            reminder = pow_mod(reminder, 2, n);
            if (reminder == 1)
                return false; // n is not a prime, and a is a witness
            if (reminder == n - 1)
                // goto NEXT_WITNESS;
                flg = true;
                break;
        }
        if(!flg){
            return false;
        }
    // NEXT_WITNESS:
        // continue;
    }
    return true; // n is *probably* prime
}
