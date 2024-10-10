#ifndef _MILLER_RABIN_H_
#define _MILLER_RABIN_H_

#include <iostream>
#include <cstdlib>
#include <random>
#include <vector>
#include <ctime>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/random.hpp>

using namespace boost::multiprecision;

/*
 *
 * Параметр k не нужен, так как рекомендуется брать k = log2(m),
 * где m проверяемое число
 *
 */
// static const uint32_t DEFAULT_ACCURACY = 5;
/*
 * Fast calculation of `a^x mod n´ by using right-to-left binary modular
 * exponentiation.
 *
 * See http://en.wikipedia.org/wiki/Modular_exponentiation
 */
uint1024_t pow_mod(uint1024_t a, uint1024_t x, uint1024_t n);
/*
 * The Miller-Rabin probabilistic primality test.
 *
 * Returns true if ``n´´ is PROBABLY prime, false if it's composite.
 * The parameter ``k´´ is the accuracy.
 *
 * The running time should be somewhere around O(k log^3 n).
 *
 */
uint1024_t rrand(uint1024_t a, uint1024_t b);
/*
 * range_rand
 *
 * Генератор больших числе в заданном диапазоне.
 * Границы справа и слева включительны.
 *
 * Для ограничения диапазона числа по количеству бит нужно задать
 * правый край числа следующим образом: pow(2, n) - 1,
 * где n - это количество бит.
 *
*/
bool is_prime(uint1024_t n);
#endif // _MILLER_RABIN_H_
