#include "pathwise_greeks.hpp"

#include <cmath>
#include <random>
#include <stdexcept>

namespace optionpricer {

PathwiseGreeks mc_pathwise_greeks(OptionType type, const OptionParams& p,
                                  std::size_t n_paths, std::uint64_t seed) {
    if (n_paths < 2)
        throw std::invalid_argument("need at least 2 paths for a standard error");
    if (p.S <= 0 || p.K <= 0 || p.sigma <= 0 || p.T <= 0)
        throw std::invalid_argument("S, K, sigma and T must be positive");

    // Same generator and seed as mc_price, so the paths are the same ones it prices.
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> Z(0.0, 1.0);

    const double sqrtT = std::sqrt(p.T);
    const double drift = (p.r - 0.5 * p.sigma * p.sigma) * p.T;
    const double vol = p.sigma * sqrtT;
    const double discount = std::exp(-p.r * p.T);
    const double sign = (type == OptionType::Call) ? 1.0 : -1.0;  // d(payoff)/d(S_T) when in the money

    double sum_d = 0.0, sum_dd = 0.0, sum_v = 0.0, sum_vv = 0.0;
    for (std::size_t i = 0; i < n_paths; ++i) {
        const double z = Z(rng);
        const double ST = p.S * std::exp(drift + vol * z);
        const bool in_the_money = (type == OptionType::Call) ? ST > p.K : ST < p.K;

        // Chain rule: dX/dtheta = e^{-rT} * d(payoff)/d(S_T) * d(S_T)/dtheta, where
        // d(payoff)/d(S_T) is +1 (call) or -1 (put) in the money and 0 otherwise, and
        //   d(S_T)/dS     = S_T / S
        //   d(S_T)/dsigma = S_T (sqrt(T) Z - sigma T)
        double d = 0.0, v = 0.0;
        if (in_the_money) {
            d = discount * sign * ST / p.S;
            v = discount * sign * ST * (sqrtT * z - p.sigma * p.T);
        }
        sum_d += d;
        sum_dd += d * d;
        sum_v += v;
        sum_vv += v * v;
    }

    // Mean and standard error, as in mc_price.
    const double n = static_cast<double>(n_paths);
    auto mean_se = [n](double sum, double sum_sq, double& mean, double& se) {
        mean = sum / n;
        se = std::sqrt((sum_sq - n * mean * mean) / (n - 1.0) / n);
    };
    PathwiseGreeks g{};
    mean_se(sum_d, sum_dd, g.delta, g.delta_se);
    mean_se(sum_v, sum_vv, g.vega, g.vega_se);
    return g;
}

}  // namespace optionpricer
