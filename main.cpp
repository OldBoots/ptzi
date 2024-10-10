#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/miller_rabin.hpp>
#include <boost/random.hpp>

#include <ctime>
#include <iostream>
#include <string>

#include "mr.h"

using namespace boost::multiprecision;
using namespace boost::random;


uint1024_t GenPrime(const size_t neededBits = 1024)
{
    static mt11213b base_gen(time(nullptr));
    static independent_bits_engine<mt11213b, 1024, uint1024_t> gen(base_gen);
    uint1024_t v = (gen() << (1024 - neededBits)) >> (1024 - neededBits);
    while(!is_prime(v)) {
        v--;
    }

    return v;
}

struct Pqg {
    uint1024_t p;
    uint1024_t q;
    uint1024_t g;
};

Pqg GenPQG()
{
    Pqg result;
    while (true) {
        result.q = GenPrime(256);
        result.p = result.q * 2 + 1;
        if (is_prime(result.p)) {
            break;
        }
    }

    uint1024_t g = 2;
    while (pow_mod(g, result.q, result.p) == 1) {
        g++;
    }

    result.g = g;

    return result;
}

int main()
{
    for(unsigned i = 0; i < 10; ++i)
    {
        auto pqg = GenPQG();
        std::cout << pqg.p << " " << pqg.q << " " << pqg.g << std::endl;
    }
}
