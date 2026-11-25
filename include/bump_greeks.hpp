// Greeks by bump-and-reprice: Monte Carlo prices at bumped inputs,
// combined with central finite differences.
#pragma once

#include "black_scholes.hpp"  // OptionType, OptionParams

#include <cstddef>
#include <cstdint>

namespace optionpricer {

struct Greeks {
    double delta;  // dV/dS
    double gamma;  // d2V/dS2
    double vega;   // dV/dsigma, per 1.00 change in sigma (same units as bs_vega)
};

// With common_random_numbers = true (the default), every repricing uses the same seed,
// so all of them see the same random numbers. Set it to false to give each repricing
// its own seed (seed, seed+1, ...). That exists only to show how much noisier it is.
Greeks mc_greeks(OptionType type, const OptionParams& p, std::size_t n_paths,
                 std::uint64_t seed, bool common_random_numbers = true);

}  // namespace optionpricer
