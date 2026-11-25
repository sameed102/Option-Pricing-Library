#include "bump_greeks.hpp"

#include "monte_carlo.hpp"

namespace optionpricer {

Greeks mc_greeks(OptionType type, const OptionParams& p, std::size_t n_paths,
                 std::uint64_t seed, bool common_random_numbers) {
    // Bump sizes: 1% of spot, and 1 vol point. Smaller bumps cut the
    // finite-difference bias but amplify Monte Carlo noise (we divide by h or h^2).
    const double hS = 0.01 * p.S;
    const double hSigma = 0.01;

    // With CRN, reprice with the same seed each time, so the only thing that changes
    // between prices is the bumped input, not the random draws. Without CRN, the
    // k-th repricing uses seed + k, so each price has its own independent noise.
    auto price = [&](const OptionParams& q, std::uint64_t k) {
        return mc_price(type, q, n_paths, common_random_numbers ? seed : seed + k).price;
    };

    OptionParams up = p, down = p;
    up.S += hS;
    down.S -= hS;
    const double V = price(p, 0);
    const double V_up = price(up, 1);
    const double V_down = price(down, 2);

    OptionParams vol_up = p, vol_down = p;
    vol_up.sigma += hSigma;
    vol_down.sigma -= hSigma;
    const double V_vol_up = price(vol_up, 3);
    const double V_vol_down = price(vol_down, 4);

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
