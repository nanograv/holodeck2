/**
 *
 */

// #include <cstdlib> // For malloc and free
// #include <cmath>
// #include <iostream>
#include <memory>   // for std::unique_ptr
#include <chrono>

#include <boost/math/distributions/poisson.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

constexpr int NUM_RUNS = 10'000'000;

using RNGType          = boost::random::mt19937;
using DistPoissonType  = boost::random::poisson_distribution<>;
using VGPoissonType    = boost::random::variate_generator<RNGType&, DistPoissonType>;
using DistNormalType   = boost::random::normal_distribution<>;
using VGNormalType     = boost::random::variate_generator<RNGType&, DistNormalType>;

RNGType rng(0);


void check_timing(double mean) {
    // printf("Timing distributions with %d runs...\n", NUM_RUNS);

    int value;
    double dur;
    std::chrono::steady_clock::time_point beg, end;

    // ---- Poisson Distribution

    // printf("Poisson Distribution\n");

    DistPoissonType dist_poisson(mean);
    VGPoissonType draw_poisson(rng, dist_poisson);
    beg = std::chrono::steady_clock::now();
    for (int i = 0; i < NUM_RUNS; ++i) {
        value = draw_poisson();
    }
    end = std::chrono::steady_clock::now();
    dur = std::chrono::duration<double>(end - beg).count();

    printf("Poisson: %f seconds\n", dur);
    // printf("\t%.4e seconds/iteration\n", dur / NUM_RUNS);


    // ---- Normal Distribution

    // printf("Normal Distribution\n");

    DistNormalType dist_normal(mean, sqrt(mean));
    VGNormalType draw_normal(rng, dist_normal);

    beg = std::chrono::steady_clock::now();
    for (int i = 0; i < NUM_RUNS; ++i) {
        value = draw_normal();
    }
    end = std::chrono::steady_clock::now();
    dur = std::chrono::duration<double>(end - beg).count();

    printf("Normal: %f seconds\n", dur);
    // printf("\t%.4e seconds/iteration\n", dur / NUM_RUNS);

}


int main(int argc, char *argv[]) {

    printf("Timing distributions with %d runs...\n", NUM_RUNS);

    for (double mean : {1.0E-12, 1.0E-6, 1.0E-3, 1.0E0, 1.0E3, 1.0E6, 1.0E12}) {
        printf("\n -- Testing with mean = %.2e\n", mean);
        check_timing(mean);
    }

    printf("\nDone with tests.\n");
    return 0;
}

