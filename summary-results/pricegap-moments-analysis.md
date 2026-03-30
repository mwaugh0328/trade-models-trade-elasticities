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

The notebook computes Pearson correlations (with p-values) between the price gap moments and three gravity variables. The expected signs are: positive for distance (more distant pairs → larger price gaps), negative for shared border (contiguity → smaller gaps), and positive for tariffs (higher barriers → larger gaps).

#### Correlations for `d_ni`

| Year | N    | corr(dni, log dist) | p-value | corr(dni, border) | p-value | corr(dni, log(1+tariff)) | p-value |
|------|------|--------------------:|--------:|------------------:|--------:|-------------------------:|--------:|
| All  | 2604 |               0.422 |  0.0000 |            -0.183 |  0.0000 |                    0.338 |  0.0000 |
| 2004 |  866 |               0.278 |  0.0000 |            -0.108 |  0.0014 |                    0.261 |  0.0000 |
| 2011 |  868 |               0.472 |  0.0000 |            -0.196 |  0.0000 |                    0.409 |  0.0000 |
| 2017 |  870 |               0.501 |  0.0000 |            -0.240 |  0.0000 |                    0.450 |  0.0000 |

#### Correlations for `d_ni2`

| Year | N    | corr(dni2, log dist) | p-value | corr(dni2, border) | p-value | corr(dni2, log(1+tariff)) | p-value |
|------|------|---------------------:|--------:|-------------------:|--------:|--------------------------:|--------:|
| All  | 2604 |                0.501 |  0.0000 |             -0.208 |  0.0000 |                     0.349 |  0.0000 |
| 2004 |  866 |                0.371 |  0.0000 |             -0.144 |  0.0000 |                     0.246 |  0.0000 |
| 2011 |  868 |                0.605 |  0.0000 |             -0.263 |  0.0000 |                     0.404 |  0.0000 |
| 2017 |  870 |                0.536 |  0.0000 |             -0.222 |  0.0000 |                     0.418 |  0.0000 |

All correlations are highly significant (p ≈ 0) and in the expected direction across every year and for both measures. Distance is the strongest correlate (0.28–0.61), followed by tariffs (0.25–0.45), with border effects more modest (-0.11 to -0.26). The `d_ni2` measure shows slightly stronger distance correlations than the `d_ni` measure, suggesting the 85th percentile may be less noisy than the max.

### 3. Regression Analysis

The notebook estimates OLS regressions of the form:

```
d_ni = α_year + β₁ · border + β₂ · log(dist) + β₃ · log(1 + 0.01 · tariff) + ε_ni
```

with year fixed effects, for both `d_ni` and `d_ni2`. Since `d_ni` is already in log-price units (i.e., `exp(d_ni)` represents the trade cost level), the dependent variable is not logged.

#### Regression: d_ni

| Variable          | Coef    | Std Err | t       | p-value | [95% CI]          |
|-------------------|--------:|--------:|--------:|--------:|-------------------|
| Intercept         | -0.235  |   0.069 |  -3.40  |   0.001 | [-0.370, -0.099]  |
| year = 2011       |  0.078  |   0.018 |   4.34  |   0.000 | [ 0.043,  0.113]  |
| year = 2017       |  0.017  |   0.018 |   0.93  |   0.354 | [-0.019,  0.052]  |
| border            | -0.032  |   0.034 |  -0.96  |   0.339 | [-0.098,  0.034]  |
| log(dist)         |  0.135  |   0.009 |  15.11  |   0.000 | [ 0.117,  0.152]  |
| log(1+tariff)     |  1.973  |   0.201 |   9.81  |   0.000 | [ 1.579,  2.367]  |

R² = 0.212, N = 2,604, F = 139.9

#### Regression: d_ni2

| Variable          | Coef    | Std Err | t       | p-value | [95% CI]          |
|-------------------|--------:|--------:|--------:|--------:|-------------------|
| Intercept         | -0.084  |   0.020 |  -4.25  |   0.000 | [-0.122, -0.045]  |
| year = 2011       | -0.027  |   0.005 |  -5.26  |   0.000 | [-0.037, -0.017]  |
| year = 2017       | -0.040  |   0.005 |  -7.76  |   0.000 | [-0.050, -0.030]  |
| border            | -0.001  |   0.010 |  -0.11  |   0.909 | [-0.020,  0.018]  |
| log(dist)         |  0.055  |   0.003 |  21.53  |   0.000 | [ 0.050,  0.060]  |
| log(1+tariff)     |  0.394  |   0.057 |   6.85  |   0.000 | [ 0.281,  0.506]  |

R² = 0.287, N = 2,604, F = 209.6

#### Interpretation

- **Distance** is highly significant in both regressions (t ≈ 15–22). A 10% increase in bilateral distance raises the price gap `d_ni` by about 0.013 (and `d_ni2` by about 0.005). Since `exp(d_ni)` represents the trade cost level, these effects translate into meaningful increases in trade costs.
- **Tariffs** are also highly significant (t ≈ 7–10). The large coefficient on `d_ni` (1.97) means that tariff variation generates substantial price gap variation — this is the key identifying variation for θ.
- **Border** is insignificant in both regressions (p = 0.34 for `d_ni`, p = 0.91 for `d_ni2`), suggesting that contiguity effects are absorbed by distance.
- **Year effects**: The 2011 dummy is significant and positive for `d_ni` (consistent with the higher mean in 2011), while both 2011 and 2017 are significantly negative for `d_ni2`, consistent with the declining trend in `d_ni2` summary statistics.
- The `d_ni2` regression has higher R² (0.287 vs. 0.212), confirming that the 85th percentile measure has a stronger signal-to-noise ratio than the max.

## Connection to the Estimation

The price gap moments analyzed here are the same objects matched in the GMM estimation (see `model-code/estimate-trade-models.jl`):

- **Exactly identified**: The estimator sets `E[d_ni^data - d_ni^model(θ)] = 0` using the `d_ni` moment.
- **Over-identified**: The estimator additionally matches `d_ni2` and `cov(d_ni, log dist)`.

The correlation and regression results in this notebook therefore provide intuition for what the moments are capturing and whether the data has sufficient variation to identify θ. Strong correlations with trade costs suggest the moments are informative; weak correlations would raise concerns about identification power.

## How to Run

Open `pricegap-moments-analysis.ipynb` in VS Code or Jupyter and run all cells. Requires Python with `pandas`, `numpy`, `scipy`, and `statsmodels`.
