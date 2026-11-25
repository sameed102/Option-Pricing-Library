// Unit tests. No framework: check() prints PASS/FAIL, and main returns 1 if any
// check failed, which is how CTest knows the test failed.
#include "black_scholes.hpp"
#include "bump_greeks.hpp"
#include "monte_carlo.hpp"

#include <cmath>
#include <cstdio>
#include <initializer_list>

using namespace optionpricer;

namespace {

int failures = 0;

void check(bool ok, const char* name, double got, double want, double tol) {
    std::printf("%s  %-40s got %12.6f  want %12.6f  tol %.2e\n",
                ok ? "PASS" : "FAIL", name, got, want, tol);
    if (!ok) ++failures;
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
void test_mc_matches_bs() {
    const std::size_t n = 100000;
    for (const OptionParams& p : cases)
        for (OptionType type : {OptionType::Call, OptionType::Put})
            for (bool antithetic : {false, true}) {
                const MCResult mc = mc_price(type, p, n, 42, antithetic);
                const double bs = bs_price(type, p);
                const double tol = 4.0 * mc.std_error;
                const char* name = antithetic ? "MC (antithetic) within 4 SE of BS"
                                              : "MC within 4 SE of BS";
                check(std::fabs(mc.price - bs) < tol, name, mc.price, bs, tol);
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
    std::printf("%d failure(s)\n", failures);
    return failures == 0 ? 0 : 1;
}
