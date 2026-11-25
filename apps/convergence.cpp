// Convergence study: prices one call by Monte Carlo at 1e3 ... 1e7 paths and
// compares each estimate to the exact Black-Scholes price.
// Prints CSV to stdout, so save it with:  ./build/convergence > convergence.csv
#include "black_scholes.hpp"
#include "monte_carlo.hpp"

#include <cmath>
#include <cstdio>

using namespace optionpricer;

int main() {
    const OptionParams p{100.0, 100.0, 0.05, 0.20, 1.0};  // S, K, r, sigma, T
    const std::uint64_t seed = 42;  // fixed so the CSV is the same on every run
    const double bs = bs_price(OptionType::Call, p);

    // av_ columns: same number of paths, with antithetic variates.
    std::printf("n_paths,mc_price,std_error,bs_price,abs_error,"
                "av_price,av_std_error,av_abs_error\n");
    for (std::size_t n = 1000; n <= 10000000; n *= 10) {
        const MCResult mc = mc_price(OptionType::Call, p, n, seed);
        const MCResult av = mc_price(OptionType::Call, p, n, seed, true);
        std::printf("%zu,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                    n, mc.price, mc.std_error, bs, std::fabs(mc.price - bs),
                    av.price, av.std_error, std::fabs(av.price - bs));
    }
    return 0;
}
