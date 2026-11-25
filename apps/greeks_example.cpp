// Compares bump-and-reprice Monte Carlo Greeks with the analytic Black-Scholes Greeks.
#include "black_scholes.hpp"
#include "bump_greeks.hpp"

#include <cmath>
#include <cstdio>
#include <initializer_list>  // for looping over {Call, Put}

using namespace optionpricer;

int main() {
    const OptionParams p{100.0, 100.0, 0.05, 0.20, 1.0};  // S, K, r, sigma, T
    const std::size_t n = 1000000;
    const std::uint64_t seed = 42;

    std::printf("S=%.2f K=%.2f r=%.2f sigma=%.2f T=%.2f  paths=%zu seed=%llu\n",
                p.S, p.K, p.r, p.sigma, p.T, n, static_cast<unsigned long long>(seed));
    std::printf("%-6s %-6s %12s %12s %12s\n", "type", "greek", "bump_mc", "analytic", "abs_diff");

    for (OptionType type : {OptionType::Call, OptionType::Put}) {
        const Greeks g = mc_greeks(type, p, n, seed);
        const char* name = type == OptionType::Call ? "call" : "put";
        const char* labels[3] = {"delta", "gamma", "vega"};
        const double bumped[3] = {g.delta, g.gamma, g.vega};
        const double analytic[3] = {bs_delta(type, p), bs_gamma(p), bs_vega(p)};
        for (int i = 0; i < 3; ++i)
            std::printf("%-6s %-6s %12.6f %12.6f %12.6f\n", name, labels[i],
                        bumped[i], analytic[i], std::fabs(bumped[i] - analytic[i]));
    }
    return 0;
}
