#include "black_scholes.hpp"

#include <cmath>
#include <stdexcept>

namespace optionpricer {

namespace {

// d1 and d2 would divide by zero if sigma or T were 0, so reject those inputs.
void check(const OptionParams& p) {
    if (p.S <= 0 || p.K <= 0 || p.sigma <= 0 || p.T <= 0)
        throw std::invalid_argument("S, K, sigma and T must be positive");
}

// d1 = [ln(S/K) + (r + sigma^2/2) T] / (sigma sqrt(T))
double d1(const OptionParams& p) {
    return (std::log(p.S / p.K) + (p.r + 0.5 * p.sigma * p.sigma) * p.T)
           / (p.sigma * std::sqrt(p.T));
}

// d2 = d1 - sigma sqrt(T)
double d2(const OptionParams& p) {
    return d1(p) - p.sigma * std::sqrt(p.T);
}

}  // namespace

// N(x) = 0.5 * erfc(-x / sqrt(2)); erfc is in <cmath> and is accurate in the tails.
double norm_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

// phi(x) = exp(-x^2/2) / sqrt(2 pi)
double norm_pdf(double x) {
    const double pi = 3.14159265358979323846;
    return std::exp(-0.5 * x * x) / std::sqrt(2.0 * pi);
}

// Call = S N(d1) - K e^{-rT} N(d2)
// Put  = K e^{-rT} N(-d2) - S N(-d1)
double bs_price(OptionType type, const OptionParams& p) {
    check(p);
    const double discK = p.K * std::exp(-p.r * p.T);  // strike discounted to today
    if (type == OptionType::Call)
        return p.S * norm_cdf(d1(p)) - discK * norm_cdf(d2(p));
    return discK * norm_cdf(-d2(p)) - p.S * norm_cdf(-d1(p));
}

// Call delta = N(d1), Put delta = N(d1) - 1
double bs_delta(OptionType type, const OptionParams& p) {
    check(p);
    if (type == OptionType::Call)
        return norm_cdf(d1(p));
    return norm_cdf(d1(p)) - 1.0;
}

// Gamma = phi(d1) / (S sigma sqrt(T))
double bs_gamma(const OptionParams& p) {
    check(p);
    return norm_pdf(d1(p)) / (p.S * p.sigma * std::sqrt(p.T));
}

// Vega = S phi(d1) sqrt(T)
double bs_vega(const OptionParams& p) {
    check(p);
    return p.S * norm_pdf(d1(p)) * std::sqrt(p.T);
}

}  // namespace optionpricer
