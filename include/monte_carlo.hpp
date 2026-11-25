// Monte Carlo pricer for European options under geometric Brownian motion.
#pragma once

#include "black_scholes.hpp"  // OptionType, OptionParams

#include <cstddef>
#include <cstdint>
#include <optional>

namespace optionpricer {

// A Monte Carlo estimate is a sample mean, so we report it with its standard error.
struct MCResult {
    double price;      // mean of the discounted payoffs
    double std_error;  // sample standard deviation / sqrt(n_paths)
};

// Simulates n_paths terminal prices and averages the discounted payoffs.
// Pass a seed to get the same answer every run; leave it empty for a random seed.
MCResult mc_price(OptionType type, const OptionParams& p, std::size_t n_paths,
                  std::optional<std::uint64_t> seed = std::nullopt);

}  // namespace optionpricer
