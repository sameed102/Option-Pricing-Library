# optionpricer

A small C++17 library for pricing European options. It has two methods so far:

1. **Black-Scholes closed form**: exact prices for calls and puts, plus analytic delta, gamma and vega.
2. **Monte Carlo**: simulates the terminal stock price under geometric Brownian motion and returns a price estimate with its standard error. Passing a seed makes a run reproducible; leaving it out uses a random seed.

It uses only the standard library (`<cmath>`, `<random>`) and builds with CMake.

## Layout

```
include/black_scholes.hpp   OptionType, OptionParams, Black-Scholes prices and Greeks
src/black_scholes.cpp
include/monte_carlo.hpp     MCResult, mc_price
src/monte_carlo.cpp
apps/bs_example.cpp         prints Black-Scholes prices and Greeks for one option
apps/mc_example.cpp         compares Monte Carlo to Black-Scholes for one option
```

## Build and run

Requires CMake 3.14 or newer and a C++17 compiler. The build defaults to Release, because Monte Carlo runs slowly in Debug.

```bash
cmake -S . -B build
```

```bash
cmake --build build
```

```bash
./build/bs_example
```

```bash
./build/mc_example
```

## Formulas

Inputs: spot $S$, strike $K$, risk-free rate $r$ (continuously compounded), volatility $\sigma$, time to expiry $T$ in years. $N$ is the standard normal CDF and $\varphi$ its density.

### Black-Scholes

$$d_1 = \frac{\ln(S/K) + (r + \tfrac{1}{2}\sigma^2)T}{\sigma\sqrt{T}}, \qquad d_2 = d_1 - \sigma\sqrt{T}$$

| Quantity | Call | Put |
|---|---|---|
| Price | $S N(d_1) - K e^{-rT} N(d_2)$ | $K e^{-rT} N(-d_2) - S N(-d_1)$ |
| Delta $\partial V/\partial S$ | $N(d_1)$ | $N(d_1) - 1$ |
| Gamma $\partial^2 V/\partial S^2$ | $\dfrac{\varphi(d_1)}{S\sigma\sqrt{T}}$ | same |
| Vega $\partial V/\partial \sigma$ | $S\varphi(d_1)\sqrt{T}$ | same |

Vega is per 1.00 change in $\sigma$. Divide by 100 to get the per-vol-point figure traders quote.

$N(x)$ is computed as $\tfrac{1}{2}\,\mathrm{erfc}(-x/\sqrt{2})$, which stays accurate in the tails.

### Monte Carlo

Under the risk-neutral measure the stock follows GBM, $dS = rS\,dt + \sigma S\,dW$. That equation has an exact solution at expiry, so each path needs only one standard normal draw $Z$:

$$S_T = S \exp\!\left((r - \tfrac{1}{2}\sigma^2)T + \sigma\sqrt{T}\,Z\right)$$

For $n$ paths with discounted payoffs $X_i = e^{-rT}\max(S_T^{(i)} - K, 0)$ (call) or $e^{-rT}\max(K - S_T^{(i)}, 0)$ (put):

$$\hat V = \frac{1}{n}\sum_{i=1}^n X_i, \qquad s^2 = \frac{\sum X_i^2 - n\hat V^2}{n-1}, \qquad \mathrm{SE} = \frac{s}{\sqrt{n}}$$

Random numbers come from `std::mt19937_64` and `std::normal_distribution`. `mc_price` takes an optional seed: with one, the same seed always gives the same price; without one, the generator is seeded from `std::random_device`.

## Results

All parameters: $S = 100$, $K = 100$, $r = 0.05$, $\sigma = 0.20$, $T = 1$.

Output of `./build/bs_example`:

```
S=100.00 K=100.00 r=0.05 sigma=0.20 T=1.00
type        price      delta      gamma       vega
call    10.450584   0.636831   0.018762  37.524035
put      5.573526  -0.363169   0.018762  37.524035
```

Output of `./build/mc_example` (100,000 paths, seed 42; running it again gives identical output):

```
S=100.00 K=100.00 r=0.05 sigma=0.20 T=1.00  paths=100000 seed=42
type     mc_price    std_err   bs_price    z_score
call    10.512032   0.046856  10.450584      1.311
put      5.549716   0.027297   5.573526     -0.872
same seed twice gives same price: yes
```

`z_score` is (MC price − BS price) / standard error. Both estimates land within 1.4 standard errors of the exact Black-Scholes price, as you would expect from sampling noise.

## Status

Done: Black-Scholes prices and Greeks; the Monte Carlo pricer.

Planned: a convergence study with CSV output, antithetic variates, bump-and-reprice Greeks with common random numbers, unit tests, and a matplotlib plot of error against path count.
