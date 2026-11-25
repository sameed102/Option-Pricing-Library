// Prices one call and one put by Monte Carlo and compares to Black-Scholes.
#include "black_scholes.hpp"
#include "monte_carlo.hpp"

#include <cstdio>

using namespace optionpricer;

int main() {
    const OptionParams p{100.0, 100.0, 0.05, 0.20, 1.0};  // S, K, r, sigma, T
    const std::size_t n = 100000;
    const std::uint64_t seed = 42;

    std::printf("S=%.2f K=%.2f r=%.2f sigma=%.2f T=%.2f  paths=%zu seed=%llu\n",
                p.S, p.K, p.r, p.sigma, p.T, n, static_cast<unsigned long long>(seed));
    std::printf("%-6s %10s %10s %10s %10s\n", "type", "mc_price", "std_err", "bs_price", "z_score");

    for (OptionType type : {OptionType::Call, OptionType::Put}) {
        const MCResult mc = mc_price(type, p, n, seed);
        const double bs = bs_price(type, p);
        // How many standard errors the estimate is from the exact answer.
        const double z = (mc.price - bs) / mc.std_error;
        std::printf("%-6s %10.6f %10.6f %10.6f %10.3f\n",
                    type == OptionType::Call ? "call" : "put", mc.price, mc.std_error, bs, z);
    }

    // Same seed => identical result. This is what makes runs reproducible.
    const MCResult a = mc_price(OptionType::Call, p, n, seed);
    const MCResult b = mc_price(OptionType::Call, p, n, seed);
    std::printf("same seed twice gives same price: %s\n", a.price == b.price ? "yes" : "no");
    return 0;
}
