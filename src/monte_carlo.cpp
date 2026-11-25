#include "monte_carlo.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

namespace optionpricer {

namespace {

// Payoff at expiry for a given terminal price.
double payoff(OptionType type, double ST, double K) {
    if (type == OptionType::Call)
        return std::max(ST - K, 0.0);
    return std::max(K - ST, 0.0);
}

}  // namespace

MCResult mc_price(OptionType type, const OptionParams& p, std::size_t n_paths,
                  std::optional<std::uint64_t> seed, bool antithetic) {
    if (antithetic && n_paths % 2 != 0)
        throw std::invalid_argument("antithetic needs an even number of paths");
    // Each sample is one path, or one (Z, -Z) pair when antithetic.
    const std::size_t n_samples = antithetic ? n_paths / 2 : n_paths;
    if (n_samples < 2)
        throw std::invalid_argument("need at least 2 samples for a standard error");
    if (p.S <= 0 || p.K <= 0 || p.sigma <= 0 || p.T <= 0)
        throw std::invalid_argument("S, K, sigma and T must be positive");

    // Mersenne Twister: fixed seed if given, otherwise seeded from the OS.
    std::mt19937_64 rng(seed ? *seed : std::random_device{}());
    std::normal_distribution<double> Z(0.0, 1.0);

    // Under the risk-neutral measure, GBM has the exact solution
    //   S_T = S * exp((r - sigma^2/2) T + sigma sqrt(T) Z),  Z ~ N(0,1),
    // so we can jump straight to expiry in one step (no time grid needed).
    const double drift = (p.r - 0.5 * p.sigma * p.sigma) * p.T;
    const double vol = p.sigma * std::sqrt(p.T);
    const double discount = std::exp(-p.r * p.T);

    // Keep running sums of X and X^2 to get the mean and variance in one pass.
    double sum = 0.0, sum_sq = 0.0;
    for (std::size_t i = 0; i < n_samples; ++i) {
        const double z = Z(rng);
        double X = discount * payoff(type, p.S * std::exp(drift + vol * z), p.K);
        if (antithetic) {
            // -z is also N(0,1). Averaging the pair cancels much of the noise, because
            // when one path ends high the other ends low.
            const double X_anti = discount * payoff(type, p.S * std::exp(drift - vol * z), p.K);
            X = 0.5 * (X + X_anti);
        }
        sum += X;
        sum_sq += X * X;
    }

    // The pair averages are i.i.d. (the two paths in a pair are not independent),
    // so the standard error uses n_samples, not n_paths.
    const double n = static_cast<double>(n_samples);
    const double mean = sum / n;
    // Unbiased sample variance: (sum X^2 - n mean^2) / (n - 1)
    const double var = (sum_sq - n * mean * mean) / (n - 1.0);
    return {mean, std::sqrt(var / n)};
}

}  // namespace optionpricer
