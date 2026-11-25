#include "bump_greeks.hpp"

#include "monte_carlo.hpp"

namespace optionpricer {

Greeks mc_greeks(OptionType type, const OptionParams& p, std::size_t n_paths,
                 std::uint64_t seed) {
    // Bump sizes: 1% of spot, and 1 vol point. Smaller bumps cut the
    // finite-difference bias but amplify Monte Carlo noise (we divide by h or h^2).
    const double hS = 0.01 * p.S;
    const double hSigma = 0.01;

    // Reprice with the same seed each time, so the only thing that changes
    // between prices is the bumped input, not the random draws.
    auto price = [&](const OptionParams& q) { return mc_price(type, q, n_paths, seed).price; };

    OptionParams up = p, down = p;
    up.S += hS;
    down.S -= hS;
    const double V = price(p);
    const double V_up = price(up);
    const double V_down = price(down);

    OptionParams vol_up = p, vol_down = p;
    vol_up.sigma += hSigma;
    vol_down.sigma -= hSigma;
    const double V_vol_up = price(vol_up);
    const double V_vol_down = price(vol_down);

    // Central differences:
    //   delta = (V(S+h) - V(S-h)) / 2h
    //   gamma = (V(S+h) - 2V(S) + V(S-h)) / h^2
    //   vega  = (V(sigma+k) - V(sigma-k)) / 2k
    return {
        (V_up - V_down) / (2.0 * hS),
        (V_up - 2.0 * V + V_down) / (hS * hS),
        (V_vol_up - V_vol_down) / (2.0 * hSigma),
    };
}

}  // namespace optionpricer
