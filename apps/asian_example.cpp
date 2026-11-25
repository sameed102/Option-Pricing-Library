// Prices an arithmetic-average Asian call by Monte Carlo for several numbers of
// monitoring dates, next to the European call with the same inputs.
#include "black_scholes.hpp"
#include "monte_carlo.hpp"

#include <cstdio>
#include <initializer_list>  // for looping over a list of step counts

using namespace optionpricer;

int main() {
    const OptionParams p{100.0, 100.0, 0.05, 0.20, 1.0};  // S, K, r, sigma, T
    const std::size_t n = 100000;
    const std::uint64_t seed = 42;

    std::printf("S=%.2f K=%.2f r=%.2f sigma=%.2f T=%.2f  paths=%zu seed=%llu\n",
                p.S, p.K, p.r, p.sigma, p.T, n, static_cast<unsigned long long>(seed));
    const MCResult euro = mc_price(OptionType::Call, p, n, seed);
    std::printf("European call: MC %.6f (SE %.6f), Black-Scholes %.6f\n",
                euro.price, euro.std_error, bs_price(OptionType::Call, p));

    // 1 date = the average is just S_T, so it must equal the European MC price.
    // 12 = monthly, 52 = weekly, 252 = daily (trading days).
    std::printf("%-8s %12s %10s\n", "n_steps", "asian_call", "std_err");
    for (std::size_t steps : {1, 12, 52, 252}) {
        const MCResult asian = mc_asian_price(OptionType::Call, p, n, steps, seed);
        std::printf("%-8zu %12.6f %10.6f\n", steps, asian.price, asian.std_error);
    }
    return 0;
}
