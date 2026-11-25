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
include/bump_greeks.hpp     Greeks, mc_greeks (bump-and-reprice)
src/bump_greeks.cpp
apps/convergence.cpp        Monte Carlo error at 1e3 ... 1e7 paths, as CSV
apps/greeks_example.cpp     compares bump Greeks to analytic Greeks
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

```bash
./build/convergence > convergence.csv
```

```bash
./build/greeks_example
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

**Antithetic variates** (`mc_price(..., seed, true)`): if $Z \sim N(0,1)$ then so is $-Z$. Each draw is used twice, and the two discounted payoffs are averaged:

$$Y_i = \tfrac{1}{2}\left(X(Z_i) + X(-Z_i)\right), \qquad i = 1, \dots, n/2$$

The estimate is the mean of the $n/2$ values $Y_i$, and its standard error is $s_Y/\sqrt{n/2}$. The two paths in a pair are not independent, but the pairs are, so the SE is counted in pairs. Both modes evaluate $n$ payoffs, so the comparison is fair. A call payoff rises with $Z$, so $X(Z)$ and $X(-Z)$ are negatively correlated, and $\mathrm{Var}(Y)$ is less than half of $\mathrm{Var}(X)$.

### Bump-and-reprice Greeks

`mc_greeks` reprices with Monte Carlo at bumped inputs and takes central differences, with $h = 0.01\,S$ and $k = 0.01$:

$$\Delta \approx \frac{V(S+h) - V(S-h)}{2h}, \qquad \Gamma \approx \frac{V(S+h) - 2V(S) + V(S-h)}{h^2}, \qquad \mathcal{V} \approx \frac{V(\sigma+k) - V(\sigma-k)}{2k}$$

All five prices use the same seed (**common random numbers**). Each path sees the same $Z$ in every repricing, so most of the Monte Carlo noise cancels in the differences. With independent seeds, the noise in each price would be divided by $h$ (or $h^2$ for gamma) and would swamp the result. Central differences have $O(h^2)$ bias. A smaller bump means less bias but more noise.

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

### Convergence

Output of `./build/convergence` (at-the-money call, seed 42, saved as `convergence.csv`):

```
n_paths,mc_price,std_error,bs_price,abs_error,av_price,av_std_error,av_abs_error
1000,10.778279,0.491342,10.450584,0.327695,10.635034,0.340532,0.184450
10000,10.629807,0.148182,10.450584,0.179223,10.516702,0.104093,0.066118
100000,10.512032,0.046856,10.450584,0.061449,10.452037,0.032809,0.001454
1000000,10.472343,0.014750,10.450584,0.021760,10.449871,0.010425,0.000713
10000000,10.450101,0.004655,10.450584,0.000482,10.454197,0.003291,0.003613
```

The `av_` columns use antithetic variates with the same number of paths.

| n_paths | plain SE | antithetic SE | ratio |
|---:|---:|---:|---:|
| 1,000 | 0.491342 | 0.340532 | 0.693 |
| 10,000 | 0.148182 | 0.104093 | 0.702 |
| 100,000 | 0.046856 | 0.032809 | 0.700 |
| 1,000,000 | 0.014750 | 0.010425 | 0.707 |
| 10,000,000 | 0.004655 | 0.003291 | 0.707 |

Antithetic variates cut the standard error by about 30% at every path count. A ratio of about 0.70 means the variance roughly halves ($0.70^2 \approx 0.49$), which matches doubling the number of plain paths at no extra cost.

Each tenfold increase in paths cuts the standard error by about 3.2, close to $\sqrt{10} \approx 3.16$, as $\mathrm{SE} = s/\sqrt{n}$ predicts. The absolute error is noisier, but stays below 1.5 standard errors at every path count. Each run uses seed 42, so a smaller run's paths are the first paths of every larger run, and the errors in different rows are correlated.

### Bump-and-reprice Greeks

Output of `./build/greeks_example` (1,000,000 paths per repricing, seed 42; running it again gives identical output):

```
S=100.00 K=100.00 r=0.05 sigma=0.20 T=1.00  paths=1000000 seed=42
type   greek       bump_mc     analytic     abs_diff
call   delta      0.637532     0.636831     0.000701
call   gamma      0.018825     0.018762     0.000063
call   vega      37.630871    37.524035     0.106837
put    delta     -0.362805    -0.363169     0.000364
put    gamma      0.018825     0.018762     0.000063
put    vega      37.452576    37.524035     0.071459
```

Every bump Greek is within 0.3% of its analytic value. The call and put gammas match exactly. On every path, call payoff minus put payoff is $S_T - K$, which is linear in $S$, so its second difference is zero.

## Status

Done: Black-Scholes prices and Greeks; the Monte Carlo pricer; the convergence study; antithetic variates; bump-and-reprice Greeks.

Planned: unit tests, and a matplotlib plot of error against path count.
