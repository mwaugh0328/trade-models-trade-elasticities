# Price Gap Moments Analysis

This document describes the analysis performed in [`pricegap-moments-analysis.ipynb`](pricegap-moments-analysis.ipynb), which examines the properties of the bilateral price gap data used in the trade elasticity estimation.

## Overview

The price gap moment `d_ni` is the key data input for identifying the trade elasticity θ. It measures bilateral price discrepancies across countries using ICP (International Comparison Program) price level data on traded goods categories. The notebook pools data across all three sample years (2004, 2011, 2017) and examines whether the price gaps correlate with standard gravity variables — distance, shared borders, and tariffs — as theory predicts they should.

## Data Construction

For each country pair (n, i) and year, the price gap is constructed as:

```
d_ni = max_k(log p_n^k - log p_i^k) - mean_k(log p_n^k - log p_i^k)
```

where `k` indexes traded goods categories (ICP "basic headings"). The `d_ni2` variant uses the 85th percentile instead of the max.

The number of basic headings varies by year:
- **2004**: 62 traded categories
- **2011**: 71 traded categories
- **2017**: 64 traded categories

The notebook merges these with bilateral gravity variables (distance, border, language) from `data/top30_gravity_data.csv` and tariff data from `data/tariffs-{year}.csv`, then filters out observations with extreme trade shares (≈ 0 or ≈ 1).

## Results

After filtering out observations with extreme trade shares (trade share ≈ 0 or ≈ 1), 2,604 of the original 2,610 observations remain.

### 1. Summary Statistics

#### `d_ni` (max price gap)

| Year | N    | N Basic Headings | Mean  | Median | Min   | Max   | Std   |
|------|------|------------------|-------|--------|-------|-------|-------|
| All  | 2604 | —                | 0.938 | 0.855  | 0.185 | 3.128 | 0.418 |
| 2004 |  866 | 62               | 0.924 | 0.898  | 0.252 | 2.268 | 0.321 |
| 2011 |  868 | 71               | 0.980 | 0.829  | 0.185 | 3.128 | 0.514 |
| 2017 |  870 | 64               | 0.911 | 0.835  | 0.206 | 2.258 | 0.394 |

The mean `d_ni` is close to 1 across all years, meaning the max log-price deviation within a country pair typically exceeds the mean deviation by about one standard deviation's worth. The 2011 sample has the widest range (max = 3.13) and highest dispersion (std = 0.51).

#### `d_ni2` (85th percentile price gap)

| Year | N    | N Basic Headings | Mean  | Median | Min   | Max   | Std   |
|------|------|------------------|-------|--------|-------|-------|-------|
| All  | 2604 | —                | 0.343 | 0.344  | 0.076 | 0.773 | 0.126 |
| 2004 |  866 | 62               | 0.369 | 0.380  | 0.100 | 0.765 | 0.119 |
| 2011 |  868 | 71               | 0.338 | 0.335  | 0.076 | 0.735 | 0.128 |
| 2017 |  870 | 64               | 0.323 | 0.310  | 0.081 | 0.773 | 0.125 |

The `d_ni2` measure is substantially smaller (mean ≈ 0.34 vs. 0.94) and less dispersed, as expected from using the 85th percentile instead of the maximum. There is a slight downward trend across years (0.37 → 0.34 → 0.32), possibly reflecting declining trade barriers over time.

### 2. Correlation Analysis

The notebook computes Pearson correlations (with p-values) between log price gap moments and three gravity variables. The expected signs are: positive for distance (more distant pairs → larger price gaps), negative for shared border (contiguity → smaller gaps), and positive for tariffs (higher barriers → larger gaps).

#### Correlations for `d_ni`

| Year | N    | corr(log dni, log dist) | p-value | corr(log dni, border) | p-value | corr(log dni, log(1+tariff)) | p-value |
|------|------|------------------------:|--------:|----------------------:|--------:|-----------------------------:|--------:|
| All  | 2604 |                   0.485 |  0.0000 |                -0.232 |  0.0000 |                        0.350 |  0.0000 |
| 2004 |  866 |                   0.308 |  0.0000 |                -0.123 |  0.0003 |                        0.265 |  0.0000 |
| 2011 |  868 |                   0.562 |  0.0000 |                -0.263 |  0.0000 |                        0.425 |  0.0000 |
| 2017 |  870 |                   0.558 |  0.0000 |                -0.294 |  0.0000 |                        0.430 |  0.0000 |

#### Correlations for `d_ni2`

| Year | N    | corr(log dni2, log dist) | p-value | corr(log dni2, border) | p-value | corr(log dni2, log(1+tariff)) | p-value |
|------|------|-------------------------:|--------:|-----------------------:|--------:|------------------------------:|--------:|
| All  | 2604 |                    0.537 |  0.0000 |                 -0.241 |  0.0000 |                         0.359 |  0.0000 |
| 2004 |  866 |                    0.384 |  0.0000 |                 -0.158 |  0.0000 |                         0.244 |  0.0000 |
| 2011 |  868 |                    0.654 |  0.0000 |                 -0.309 |  0.0000 |                         0.433 |  0.0000 |
| 2017 |  870 |                    0.575 |  0.0000 |                 -0.255 |  0.0000 |                         0.423 |  0.0000 |

All correlations are highly significant (p ≈ 0) and in the expected direction across every year and for both measures. Distance is the strongest correlate (0.31–0.65), followed by tariffs (0.24–0.43), with border effects more modest (-0.12 to -0.31). The `d_ni2` measure shows slightly stronger distance correlations than the `d_ni` measure, suggesting the 85th percentile may be less noisy than the max.

### 3. Regression Analysis

The notebook estimates OLS regressions of the form:

```
log(d_ni) = α_year + β₁ · border + β₂ · log(dist) + β₃ · log(1 + 0.01 · tariff) + ε_ni
```

with year fixed effects, for both `d_ni` and `d_ni2`.

#### Regression: log(d_ni)

| Variable          | Coef    | Std Err | t       | p-value | [95% CI]          |
|-------------------|--------:|--------:|--------:|--------:|-------------------|
| Intercept         | -1.555  |   0.070 | -22.34  |   0.000 | [-1.692, -1.419]  |
| year = 2011       |  0.020  |   0.018 |   1.11  |   0.267 | [-0.015,  0.056]  |
| year = 2017       | -0.019  |   0.018 |  -1.03  |   0.303 | [-0.055,  0.017]  |
| border            | -0.077  |   0.034 |  -2.28  |   0.023 | [-0.143, -0.011]  |
| log(dist)         |  0.169  |   0.009 |  18.74  |   0.000 | [ 0.151,  0.186]  |
| log(1+tariff)     |  1.755  |   0.203 |   8.66  |   0.000 | [ 1.358,  2.153]  |

R² = 0.260, N = 2,604, F = 182.6

#### Regression: log(d_ni2)

| Variable          | Coef    | Std Err | t       | p-value | [95% CI]          |
|-------------------|--------:|--------:|--------:|--------:|-------------------|
| Intercept         | -2.615  |   0.062 | -42.05  |   0.000 | [-2.737, -2.493]  |
| year = 2011       | -0.097  |   0.016 |  -5.97  |   0.000 | [-0.128, -0.065]  |
| year = 2017       | -0.138  |   0.016 |  -8.44  |   0.000 | [-0.170, -0.106]  |
| border            | -0.037  |   0.030 |  -1.24  |   0.216 | [-0.097,  0.022]  |
| log(dist)         |  0.190  |   0.008 |  23.65  |   0.000 | [ 0.174,  0.206]  |
| log(1+tariff)     |  1.199  |   0.181 |   6.62  |   0.000 | [ 0.844,  1.554]  |

R² = 0.325, N = 2,604, F = 250.4

#### Interpretation

- **Distance** is highly significant in both regressions (t ≈ 19–24). A 10% increase in bilateral distance raises the price gap by about 1.7–1.9%.
- **Tariffs** are also highly significant (t ≈ 7–9). The large coefficients (1.2–1.8) mean that tariff variation generates substantial price gap variation — this is the key identifying variation for θ.
- **Border** is marginally significant for `d_ni` (p = 0.023) and insignificant for `d_ni2` (p = 0.216), suggesting that contiguity effects are mostly absorbed by distance.
- **Year effects** are insignificant for `d_ni` but highly significant for `d_ni2`, with negative coefficients for 2011 and 2017 relative to 2004, consistent with the declining trend in `d_ni2` summary statistics.
- The `d_ni2` regression has higher R² (0.325 vs. 0.260), confirming that the 85th percentile measure has a stronger signal-to-noise ratio than the max.

## Connection to the Estimation

The price gap moments analyzed here are the same objects matched in the GMM estimation (see `model-code/estimate-trade-models.jl`):

- **Exactly identified**: The estimator sets `E[d_ni^data - d_ni^model(θ)] = 0` using the `d_ni` moment.
- **Over-identified**: The estimator additionally matches `d_ni2` and `cov(d_ni, log dist)`.

The correlation and regression results in this notebook therefore provide intuition for what the moments are capturing and whether the data has sufficient variation to identify θ. Strong correlations with trade costs suggest the moments are informative; weak correlations would raise concerns about identification power.

## How to Run

Open `pricegap-moments-analysis.ipynb` in VS Code or Jupyter and run all cells. Requires Python with `pandas`, `numpy`, `scipy`, and `statsmodels`.
