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
// antithetic = true pairs each draw Z with -Z: n_paths/2 draws, n_paths paths in total
// (n_paths must then be even), so both modes cost the same number of payoffs.
// control_variate = true corrects the estimate using the discounted stock price,
// whose expectation is known exactly (it is S).
MCResult mc_price(OptionType type, const OptionParams& p, std::size_t n_paths,
                  std::optional<std::uint64_t> seed = std::nullopt,
                  bool antithetic = false, bool control_variate = false);

// Arithmetic-average Asian option: the payoff uses the average of the stock price on
// n_steps equally spaced dates (T/n_steps, 2T/n_steps, ..., T) instead of S_T.
MCResult mc_asian_price(OptionType type, const OptionParams& p, std::size_t n_paths,
                        std::size_t n_steps,
                        std::optional<std::uint64_t> seed = std::nullopt);

}  // namespace optionpricer
