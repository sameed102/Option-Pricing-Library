// Compares three ways of getting Greeks from Monte Carlo with the analytic Black-Scholes
// Greeks: bump-and-reprice with common random numbers, bump-and-reprice with independent
// seeds, and pathwise derivatives.
#include "black_scholes.hpp"
#include "bump_greeks.hpp"
#include "pathwise_greeks.hpp"

#include <cstdio>
#include <initializer_list>  // for looping over {Call, Put}

using namespace optionpricer;

int main() {
    const OptionParams p{100.0, 100.0, 0.05, 0.20, 1.0};  // S, K, r, sigma, T
    const std::size_t n = 1000000;
    const std::uint64_t seed = 42;

    std::printf("S=%.2f K=%.2f r=%.2f sigma=%.2f T=%.2f  paths=%zu seed=%llu\n",
                p.S, p.K, p.r, p.sigma, p.T, n, static_cast<unsigned long long>(seed));
    std::printf("%-5s %-6s %11s %11s %12s %11s %9s\n",
                "type", "greek", "analytic", "bump_crn", "bump_no_crn", "pathwise", "pw_se");

    for (OptionType type : {OptionType::Call, OptionType::Put}) {
        const char* name = type == OptionType::Call ? "call" : "put";
        const Greeks crn = mc_greeks(type, p, n, seed);
        const Greeks indep = mc_greeks(type, p, n, seed, false);
        const PathwiseGreeks pw = mc_pathwise_greeks(type, p, n, seed);

        std::printf("%-5s %-6s %11.6f %11.6f %12.6f %11.6f %9.6f\n", name, "delta",
                    bs_delta(type, p), crn.delta, indep.delta, pw.delta, pw.delta_se);
        // No pathwise gamma: the payoff's slope jumps at the strike (see pathwise_greeks.hpp).
        std::printf("%-5s %-6s %11.6f %11.6f %12.6f %11s %9s\n", name, "gamma",
                    bs_gamma(p), crn.gamma, indep.gamma, "n/a", "n/a");
        std::printf("%-5s %-6s %11.6f %11.6f %12.6f %11.6f %9.6f\n", name, "vega",
                    bs_vega(p), crn.vega, indep.vega, pw.vega, pw.vega_se);
    }
    return 0;
}
