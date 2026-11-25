// Prints Black-Scholes prices and Greeks for one at-the-money option.
#include "black_scholes.hpp"

#include <cstdio>

using namespace optionpricer;

int main() {
    const OptionParams p{100.0, 100.0, 0.05, 0.20, 1.0};  // S, K, r, sigma, T

    std::printf("S=%.2f K=%.2f r=%.2f sigma=%.2f T=%.2f\n", p.S, p.K, p.r, p.sigma, p.T);
    std::printf("%-6s %10s %10s %10s %10s\n", "type", "price", "delta", "gamma", "vega");
    std::printf("%-6s %10.6f %10.6f %10.6f %10.6f\n", "call",
                bs_price(OptionType::Call, p), bs_delta(OptionType::Call, p),
                bs_gamma(p), bs_vega(p));
    std::printf("%-6s %10.6f %10.6f %10.6f %10.6f\n", "put",
                bs_price(OptionType::Put, p), bs_delta(OptionType::Put, p),
                bs_gamma(p), bs_vega(p));

    // Round trip: feed the call price back in and recover sigma.
    const double iv = implied_vol(OptionType::Call, bs_price(OptionType::Call, p),
                                  p.S, p.K, p.r, p.T);
    std::printf("implied vol from the call price: %.10f\n", iv);
    return 0;
}
