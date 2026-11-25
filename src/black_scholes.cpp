#include "black_scholes.hpp"

#include <algorithm>
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

double implied_vol(OptionType type, double price, double S, double K, double r, double T) {
    // No-arbitrage bounds. The BS price rises with sigma from the lower bound
    // (sigma -> 0) to the upper bound (sigma -> infinity), so a solution exists only
    // strictly between them, and it is unique.
    //   call: max(S - K e^{-rT}, 0) < C < S
    //   put:  max(K e^{-rT} - S, 0) < P < K e^{-rT}
    const double discK = K * std::exp(-r * T);
    const double lower = (type == OptionType::Call) ? std::max(S - discK, 0.0)
                                                    : std::max(discK - S, 0.0);
    const double upper = (type == OptionType::Call) ? S : discK;
    if (!(price > lower && price < upper))
        throw std::invalid_argument("price is outside the no-arbitrage bounds");

    // Newton's method on f(sigma) = bs_price(sigma) - price, whose derivative is vega:
    //   sigma_new = sigma - f(sigma) / vega(sigma)
    OptionParams p{S, K, r, 0.2, T};  // start from 20% vol
    for (int i = 0; i < 100; ++i) {
        const double diff = bs_price(type, p) - price;
        if (std::fabs(diff) < 1e-12 * price)
            return p.sigma;
        const double next = p.sigma - diff / bs_vega(p);
        // shortcut: a step that would make sigma negative halves it instead; fine for
        // normal inputs, switch to a bracketed method (e.g. Brent) for extreme ones.
        p.sigma = next > 0 ? next : 0.5 * p.sigma;
    }
    throw std::runtime_error("implied_vol: Newton's method did not converge");
}

}  // namespace optionpricer
