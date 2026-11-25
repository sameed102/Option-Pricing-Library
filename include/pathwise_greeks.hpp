// Pathwise Greeks: differentiate each simulated discounted payoff with respect to the
// input, then average. No bump, so no finite-difference bias.
#pragma once

#include "black_scholes.hpp"  // OptionType, OptionParams

#include <cstddef>
#include <cstdint>

namespace optionpricer {

// Each Greek is a Monte Carlo average, so it comes with a standard error.
// There is no pathwise gamma: the payoff's first derivative jumps at the strike,
// so its second derivative is zero on every path and the average would be 0.
struct PathwiseGreeks {
    double delta, delta_se;
    double vega, vega_se;  // per 1.00 change in sigma, like bs_vega
};

PathwiseGreeks mc_pathwise_greeks(OptionType type, const OptionParams& p,
                                  std::size_t n_paths, std::uint64_t seed);

}  // namespace optionpricer
