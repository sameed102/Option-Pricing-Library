// Unit tests. No framework: check() prints PASS/FAIL, and main returns 1 if any
// check failed, which is how CTest knows the test failed.
#include "black_scholes.hpp"
#include "bump_greeks.hpp"
#include "monte_carlo.hpp"
#include "pathwise_greeks.hpp"

#include <cmath>
#include <cstdio>
#include <exception>
#include <initializer_list>

using namespace optionpricer;

namespace {

int failures = 0;

void check(bool ok, const char* name, double got, double want, double tol) {
    std::printf("%s  %-40s got %12.6f  want %12.6f  tol %.2e\n",
                ok ? "PASS" : "FAIL", name, got, want, tol);
    if (!ok) ++failures;
}

// Passes if calling f() throws an exception (used for invalid inputs).
template <class F>
void check_throws(F f, const char* name) {
    bool threw = false;
    try {
        f();
    } catch (const std::exception&) {
        threw = true;
    }
    std::printf("%s  %s\n", threw ? "PASS" : "FAIL", name);
    if (!threw) ++failures;
}

// A few option setups: at the money, in the money, out of the money.
const OptionParams cases[] = {
    {100.0, 100.0, 0.05, 0.20, 1.0},
    {110.0, 100.0, 0.03, 0.25, 0.5},
    {90.0, 100.0, 0.01, 0.30, 2.0},
};

// Put-call parity: C - P = S - K e^{-rT}. It holds for any model, so the
// Black-Scholes prices must satisfy it up to floating-point rounding.
void test_put_call_parity() {
    for (const OptionParams& p : cases) {
        const double lhs = bs_price(OptionType::Call, p) - bs_price(OptionType::Put, p);
        const double rhs = p.S - p.K * std::exp(-p.r * p.T);
        const double tol = 1e-10;
        check(std::fabs(lhs - rhs) < tol, "put-call parity", lhs, rhs, tol);
    }
}

// MC is random, so it can't match BS exactly. If the pricer is right, the error is
// roughly N(0, SE^2), and |error| < 4 SE fails only about 6 times in 100,000.
// Run for each mode: plain, antithetic, control variate.
void test_mc_matches_bs() {
    const std::size_t n = 100000;
    const char* names[3] = {"MC within 4 SE of BS", "MC (antithetic) within 4 SE of BS",
                            "MC (control var.) within 4 SE of BS"};
    for (const OptionParams& p : cases)
        for (OptionType type : {OptionType::Call, OptionType::Put})
            for (int mode = 0; mode < 3; ++mode) {
                const MCResult mc = mc_price(type, p, n, 42, mode == 1, mode == 2);
                const double bs = bs_price(type, p);
                const double tol = 4.0 * mc.std_error;
                check(std::fabs(mc.price - bs) < tol, names[mode], mc.price, bs, tol);
            }
}

// Antithetic pairs paths, so an odd path count is an input error, not something
// to round silently.
void test_invalid_inputs() {
    const OptionParams p = cases[0];
    check_throws([&] { mc_price(OptionType::Call, p, 1001, 42, true); },
                 "antithetic rejects an odd number of paths");
    check_throws([&] { implied_vol(OptionType::Call, 150.0, p.S, p.K, p.r, p.T); },
                 "implied vol rejects a call price above S");
}

// Implied vol must undo bs_price: price at a known sigma, invert, get sigma back.
void test_implied_vol_round_trip() {
    for (const OptionParams& c : cases)
        for (OptionType type : {OptionType::Call, OptionType::Put})
            for (double sigma : {0.1, 0.3, 0.6}) {
                OptionParams p = c;
                p.sigma = sigma;
                const double iv = implied_vol(type, bs_price(type, p), p.S, p.K, p.r, p.T);
                const double tol = 1e-8;
                check(std::fabs(iv - sigma) < tol, "implied vol round trip", iv, sigma, tol);
            }
}

// With one monitoring date the average is S_T, so the Asian pricer must reproduce
// the European MC price on the same seed. With many dates, averaging lowers the
// volatility of the underlying, so the Asian call must be cheaper than the European.
void test_asian() {
    const OptionParams p = cases[0];
    const std::size_t n = 100000;
    const double euro_mc = mc_price(OptionType::Call, p, n, 42).price;
    const double asian1 = mc_asian_price(OptionType::Call, p, n, 1, 42).price;
    check(std::fabs(asian1 - euro_mc) < 1e-12, "Asian with 1 date equals European MC",
          asian1, euro_mc, 1e-12);

    const MCResult asian = mc_asian_price(OptionType::Call, p, n, 252, 42);
    const double euro_bs = bs_price(OptionType::Call, p);
    check(asian.price + 4.0 * asian.std_error < euro_bs,
          "Asian call (+4 SE) below European call", asian.price, euro_bs, 0.0);
}

// Pathwise Greeks are unbiased averages with their own SE, so use the 4 SE rule.
void test_pathwise_greeks() {
    const std::size_t n = 1000000;
    for (const OptionParams& p : cases)
        for (OptionType type : {OptionType::Call, OptionType::Put}) {
            const PathwiseGreeks g = mc_pathwise_greeks(type, p, n, 42);
            const double want_d = bs_delta(type, p), want_v = bs_vega(p);
            check(std::fabs(g.delta - want_d) < 4.0 * g.delta_se,
                  "pathwise delta within 4 SE", g.delta, want_d, 4.0 * g.delta_se);
            check(std::fabs(g.vega - want_v) < 4.0 * g.vega_se,
                  "pathwise vega within 4 SE", g.vega, want_v, 4.0 * g.vega_se);
        }
}

// Bump Greeks carry MC noise and finite-difference bias, so compare within 1%
// of the analytic value.
void test_bump_greeks() {
    const std::size_t n = 1000000;
    for (const OptionParams& p : cases)
        for (OptionType type : {OptionType::Call, OptionType::Put}) {
            const Greeks g = mc_greeks(type, p, n, 42);
            const double want[3] = {bs_delta(type, p), bs_gamma(p), bs_vega(p)};
            const double got[3] = {g.delta, g.gamma, g.vega};
            const char* names[3] = {"bump delta within 1%", "bump gamma within 1%",
                                    "bump vega within 1%"};
            for (int i = 0; i < 3; ++i) {
                const double tol = 0.01 * std::fabs(want[i]);
                check(std::fabs(got[i] - want[i]) < tol, names[i], got[i], want[i], tol);
            }
        }
}

}  // namespace

int main() {
    test_put_call_parity();
    test_mc_matches_bs();
    test_bump_greeks();
    test_invalid_inputs();
    test_implied_vol_round_trip();
    test_asian();
    test_pathwise_greeks();
    std::printf("%d failure(s)\n", failures);
    return failures == 0 ? 0 : 1;
}
