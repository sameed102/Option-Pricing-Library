#include "monte_carlo.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

namespace optionpricer {

namespace {

// Payoff at expiry for a given terminal (or average) price.
double payoff(OptionType type, double ST, double K) {
    if (type == OptionType::Call)
        return std::max(ST - K, 0.0);
    return std::max(K - ST, 0.0);
}

void check_params(const OptionParams& p) {
    if (p.S <= 0 || p.K <= 0 || p.sigma <= 0 || p.T <= 0)
        throw std::invalid_argument("S, K, sigma and T must be positive");
}

// Mersenne Twister: fixed seed if given, otherwise seeded from the OS.
std::mt19937_64 make_rng(std::optional<std::uint64_t> seed) {
    return std::mt19937_64(seed ? *seed : std::random_device{}());
}

// Mean and standard error of n i.i.d. samples, from the running sums of X and X^2.
MCResult mean_and_se(double sum, double sum_sq, double n) {
    const double mean = sum / n;
    // Unbiased sample variance: (sum X^2 - n mean^2) / (n - 1)
    const double var = (sum_sq - n * mean * mean) / (n - 1.0);
    return {mean, std::sqrt(var / n)};
}

}  // namespace

MCResult mc_price(OptionType type, const OptionParams& p, std::size_t n_paths,
                  std::optional<std::uint64_t> seed, bool antithetic,
                  bool control_variate) {
    if (antithetic && n_paths % 2 != 0)
        throw std::invalid_argument("antithetic needs an even number of paths");
    // Each sample is one path, or one (Z, -Z) pair when antithetic.
    const std::size_t n_samples = antithetic ? n_paths / 2 : n_paths;
    if (n_samples < 2)
        throw std::invalid_argument("need at least 2 samples for a standard error");
    check_params(p);

    std::mt19937_64 rng = make_rng(seed);
    std::normal_distribution<double> Z(0.0, 1.0);

    // Under the risk-neutral measure, GBM has the exact solution
    //   S_T = S * exp((r - sigma^2/2) T + sigma sqrt(T) Z),  Z ~ N(0,1),
    // so we can jump straight to expiry in one step (no time grid needed).
    const double drift = (p.r - 0.5 * p.sigma * p.sigma) * p.T;
    const double vol = p.sigma * std::sqrt(p.T);
    const double discount = std::exp(-p.r * p.T);

    // X = discounted payoff (what we want), C = discounted S_T (the control).
    // Running sums give means, variances and the covariance in one pass.
    double sum_x = 0.0, sum_xx = 0.0, sum_c = 0.0, sum_cc = 0.0, sum_xc = 0.0;
    for (std::size_t i = 0; i < n_samples; ++i) {
        const double z = Z(rng);
        double ST = p.S * std::exp(drift + vol * z);
        double X = discount * payoff(type, ST, p.K);
        double C = discount * ST;
        if (antithetic) {
            // -z is also N(0,1). Averaging the pair cancels much of the noise, because
            // when one path ends high the other ends low.
            const double ST_anti = p.S * std::exp(drift - vol * z);
            X = 0.5 * (X + discount * payoff(type, ST_anti, p.K));
            C = 0.5 * (C + discount * ST_anti);
        }
        sum_x += X;
        sum_xx += X * X;
        sum_c += C;
        sum_cc += C * C;
        sum_xc += X * C;
    }

    // The samples (paths, or pair averages) are i.i.d.; the two paths in a pair are
    // not, so the standard error uses n_samples, not n_paths.
    const double n = static_cast<double>(n_samples);
    if (!control_variate)
        return mean_and_se(sum_x, sum_xx, n);

    // Control variate: the discounted stock is a martingale under the risk-neutral
    // measure, so E[C] = S exactly. If this sample's mean(C) came out above S, mean(X)
    // probably did too (X and C are correlated), so we correct it:
    //   price = mean(X) - beta * (mean(C) - S),  beta = Cov(X, C) / Var(C).
    // This beta is the slope of a regression of X on C, and it minimises the variance.
    const double mean_x = sum_x / n, mean_c = sum_c / n;
    const double var_x = (sum_xx - n * mean_x * mean_x) / (n - 1.0);
    const double var_c = (sum_cc - n * mean_c * mean_c) / (n - 1.0);
    const double cov_xc = (sum_xc - n * mean_x * mean_c) / (n - 1.0);
    const double beta = cov_xc / var_c;
    // Var(X - beta C) = Var(X) - Cov(X,C)^2 / Var(C) = Var(X) (1 - rho^2).
    // shortcut: beta is estimated from the same sample, which adds an O(1/n) bias;
    // negligible at these path counts, use a separate pilot run if it ever matters.
    const double var_adj = var_x - cov_xc * cov_xc / var_c;
    return {mean_x - beta * (mean_c - p.S), std::sqrt(var_adj / n)};
}

MCResult mc_asian_price(OptionType type, const OptionParams& p, std::size_t n_paths,
                        std::size_t n_steps, std::optional<std::uint64_t> seed) {
    if (n_paths < 2)
        throw std::invalid_argument("need at least 2 paths for a standard error");
    if (n_steps < 1)
        throw std::invalid_argument("need at least 1 time step");
    check_params(p);

    std::mt19937_64 rng = make_rng(seed);
    std::normal_distribution<double> Z(0.0, 1.0);

    // The payoff depends on the whole path, so we step through time. Each step uses
    // the same exact GBM solution as above, over dt instead of T.
    const double dt = p.T / static_cast<double>(n_steps);
    const double drift = (p.r - 0.5 * p.sigma * p.sigma) * dt;
    const double vol = p.sigma * std::sqrt(dt);
    const double discount = std::exp(-p.r * p.T);

    double sum = 0.0, sum_sq = 0.0;
    for (std::size_t i = 0; i < n_paths; ++i) {
        double S = p.S, total = 0.0;
        for (std::size_t j = 0; j < n_steps; ++j) {
            S *= std::exp(drift + vol * Z(rng));  // price at the next monitoring date
            total += S;
        }
        const double X = discount * payoff(type, total / static_cast<double>(n_steps), p.K);
        sum += X;
        sum_sq += X * X;
    }
    return mean_and_se(sum, sum_sq, static_cast<double>(n_paths));
}

}  // namespace optionpricer
