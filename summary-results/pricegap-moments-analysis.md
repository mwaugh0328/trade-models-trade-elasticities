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

## Analysis Sections

### 1. Summary Statistics

The notebook reports summary statistics (mean, median, min, max, standard deviation) for both `d_ni` and `d_ni2` by year and pooled across all years. These characterize the distribution of price gaps in the data and help assess whether there are significant differences across sample periods.

### 2. Correlation Analysis

The notebook computes Pearson correlations (with p-values and 90% confidence intervals via Fisher z-transformation) between the log price gap moments and three gravity variables:

| Correlation | Expected Sign | Rationale |
|------------|--------------|-----------|
| `corr(log d_ni, log dist)` | Positive | Greater distance → higher trade costs → larger price gaps |
| `corr(log d_ni, border)` | Negative | Shared border → lower trade costs → smaller price gaps |
| `corr(log d_ni, log(1 + tariff))` | Positive | Higher tariffs → higher trade barriers → larger price gaps |

These correlations serve as a sanity check that the price gap data behaves consistently with trade theory. Significant correlations in the expected direction support using `d_ni` as a moment for estimating trade costs — if price gaps were pure noise, they would not correlate systematically with observables.

The analysis is performed for both `d_ni` (max price gap) and `d_ni2` (85th percentile), and separately by year as well as pooled.

### 3. Regression Analysis

The notebook estimates OLS regressions of the form:

```
log(d_ni) = α_year + β₁ · border + β₂ · log(dist) + β₃ · log(1 + 0.01 · tariff) + ε_ni
```

with year fixed effects, for both `d_ni` and `d_ni2`. This is a multivariate version of the correlation analysis that controls for year effects and examines the partial contribution of each gravity variable.

Key questions addressed:
- **Does distance predict price gaps?** A positive β₂ means more distant country pairs have larger price discrepancies, consistent with iceberg trade costs.
- **Do tariffs predict price gaps?** A positive β₃ provides direct evidence that policy trade barriers translate into price differences, which is the mechanism that identifies θ from tariff variation.
- **Are results stable across years?** Year fixed effects absorb level shifts, so the coefficients reflect within-year cross-sectional variation.

## Connection to the Estimation

The price gap moments analyzed here are the same objects matched in the GMM estimation (see `model-code/estimate-trade-models.jl`):

- **Exactly identified**: The estimator sets `E[d_ni^data - d_ni^model(θ)] = 0` using the `d_ni` moment.
- **Over-identified**: The estimator additionally matches `d_ni2` and `cov(d_ni, log dist)`.

The correlation and regression results in this notebook therefore provide intuition for what the moments are capturing and whether the data has sufficient variation to identify θ. Strong correlations with trade costs suggest the moments are informative; weak correlations would raise concerns about identification power.

## How to Run

Open `pricegap-moments-analysis.ipynb` in VS Code or Jupyter and run all cells. Requires Python with `pandas`, `numpy`, `scipy`, and `statsmodels`.
