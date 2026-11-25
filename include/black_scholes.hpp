// Black-Scholes closed-form prices and analytic Greeks for European options.
#pragma once

namespace optionpricer {

// Call pays max(S_T - K, 0); Put pays max(K - S_T, 0).
enum class OptionType { Call, Put };

// The five Black-Scholes inputs, grouped so every function takes the same thing.
struct OptionParams {
    double S;      // spot price today
    double K;      // strike price
    double r;      // continuously compounded risk-free rate (0.05 = 5%)
    double sigma;  // annualised volatility (0.20 = 20%)
    double T;      // time to expiry in years
};

// Standard normal CDF, N(x) = P(Z <= x).
double norm_cdf(double x);

// Standard normal PDF, phi(x).
double norm_pdf(double x);

// Option price.
double bs_price(OptionType type, const OptionParams& p);

// dPrice/dS.
double bs_delta(OptionType type, const OptionParams& p);

// d2Price/dS2 (same for calls and puts).
double bs_gamma(const OptionParams& p);

// dPrice/dsigma per 1.00 change in sigma (same for calls and puts).
// Divide by 100 to get the "per 1 vol point" number traders quote.
double bs_vega(const OptionParams& p);

}  // namespace optionpricer
