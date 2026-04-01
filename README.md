# Trade Models, Trade Elasticities, and the Gains from Trade

Repository for **[Trade Models, Trade Elasticities, and the Gains from Trade](trade-elasticities.pdf)** (JIE-ISOT version).

This codebase estimates the trade elasticity (θ) — a key parameter governing the gains from trade — under four different structural trade models using international price gap data and bilateral trade flows. The key insight is that different trade models imply different mappings from θ to observable price gaps, so the same data yields different estimates depending on the assumed model.

**Key documents:**
- [model-description.md](model-description.md) — economic environment and simulation structure for each model
- [summary-results/estimation-results-summary.md](summary-results/estimation-results-summary.md) — all estimation results
- [summary-results/pricegap-moments-analysis.md](summary-results/pricegap-moments-analysis.md) — price gap data properties and correlations

## Models

| Model | Key Feature | Trade Elasticity |
|-------|-------------|-----------------|
| **Eaton-Kortum (EK)** | Ricardian comparative advantage with Fréchet productivity draws | θ (Fréchet shape) |
| **BEJK** | Bertrand competition with two productivity draws per firm | θ (Fréchet shape), σ (elasticity of substitution) |
| **Krugman** | Monopolistic competition, love of variety, no comparative advantage | θ (governs CES elasticity η = θ + 1) |
| **Melitz** | Heterogeneous firms with export selection (only productive firms export) | θ (Pareto shape), σ (CES elasticity) |

Krugman also has a "Model 2" variant (Krugman-model2) that incorporates asymmetric country size effects on the extensive margin.

## Repository Structure

```
├── build-data-code/       # Data construction pipeline
│   ├── make_price_data.jl           # Main driver: builds price gap and trade share data (2004, 2011)
│   ├── make2017_data.jl             # Builds 2017 data (uses BACI instead of Feenstra)
│   ├── functions_make_trade_share.jl # Shared functions for trade share construction
│   ├── clean-trade-data.ipynb       # Cleans BACI trade data for 2017
│   ├── clean-grossoutput-data.ipynb # Cleans UNIDO gross output data
│   └── ...                         # Raw data subdirectories
│
├── data/                  # Processed datasets used in estimation
│   ├── pricegap-df-{year}.csv       # Bilateral price gap moments (d_ni) by year
│   ├── tradeshare-df-{year}.csv     # Bilateral trade share matrices by year
│   ├── tariffs-{year}.csv           # Bilateral tariff data
│   ├── top30_gravity_data.csv       # Gravity variables (distance, border, etc.)
│   └── top30-codes-names.csv        # Country codes and names for the 30-country sample
│
├── model-code/            # Estimation and simulation code
│   ├── gravity-tools.jl             # Gravity equation estimation and trade cost construction
│   ├── trade-environment.jl         # Trade equilibrium solver (wages, trade balance)
│   ├── simulate-trade-models.jl     # Monte Carlo simulation of all four trade models
│   ├── estimate-trade-models.jl     # GMM estimation: moment computation and optimization
│   ├── estimate-ek-model.jl         # Runner script for EK estimation
│   ├── estimate-bejk-model.jl       # Runner script for BEJK estimation
│   ├── estimate-krugman-model.jl    # Runner script for Krugman estimation
│   ├── estimate-melitz-model.jl     # Runner script for Melitz estimation
│   ├── estimate-krugman-model2.jl   # Runner script for Krugman Model 2 variant
│   ├── estimate-bejk-model-sigma-loop.jl # BEJK estimation over multiple σ values
│   ├── pricegap-moments.jl          # Exploratory analysis of price gap moments
│   ├── pricegap-moments-allyears.jl # Price gap analysis pooled across years
│   └── run-estimation-all-models.jl # Master script to run all estimations
│
├── results/               # Raw estimation output (CSV files per model/method/run)
│
├── summary-results/       # Aggregated results across models and years
│   ├── combined-exact-results.csv           # Just-identified estimates (all models)
│   ├── combined-over-results.csv            # Over-identified estimates (all models)
│   ├── combined-bejk-sigma-*.csv            # BEJK sensitivity to σ
│   ├── combined-krugman-model2-*.csv        # Krugman Model 2 results
│   ├── estimation-results-summary.md        # Tabular summary of all results
│   ├── pricegap-moments-analysis.md         # Discussion of price gap data properties
│   ├── pricegap-moments-analysis.ipynb      # Notebook: correlations and regressions on price gaps
│   └── combine-results.ipynb               # Notebook that aggregates result CSVs
│
└── slides/                # Presentation materials
```

## Methodology

### Step 1: Data Construction (`build-data-code/`)

The pipeline constructs two key objects for each year (2004, 2011, 2017):

- **Price gap moments** (`d_ni`): Using ICP price level data across ~60-70 traded goods categories, the code computes bilateral price gap statistics. For each country pair (n, i), `d_ni` = max(log p_n - log p_i) - mean(log p_n - log p_i), capturing the largest price discrepancy relative to average price differences. Additional moments (`d_ni2`, `d_ni3`) use the 85th and 75th percentiles instead of the max.

- **Trade shares** (`π_ni`): Bilateral trade flows (Feenstra for 2004/2011, BACI for 2017) divided by total absorption (gross output + imports - exports) give expenditure shares.

The sample consists of the 30 largest manufacturing economies. For full details on the data pipeline, sources, and cleaning steps, see the [build-data-code README](build-data-code/README.md).

### Step 2: Gravity Estimation (`gravity-tools.jl`)

A standard gravity equation is estimated on log normalized trade shares:

```
log(π_ni / π_nn) = S_i + θ_m,n + Σ_k β_k · dist_bin_k + γ · border + ε_ni
```

This recovers:
- **Technology parameters** (`S_i`): Exporter fixed effects capturing comparative advantage
- **Trade costs** (`d_ni`): Constructed from distance bins, border effects, and asymmetric destination effects
- **Error variance** (`σ_ν`): Used in the parametric bootstrap

### Step 3: Simulation (`simulate-trade-models.jl`)

Given θ, trade costs `d`, and technology parameters `S`, each model is simulated via Monte Carlo:

- **EK**: Draw Fréchet productivities, find lowest-cost supplier for each good in each country
- **BEJK**: Draw two Fréchet productivities per firm, apply Bertrand pricing (price = min of second-best cost, monopoly markup)
- **Krugman**: All firms charge CES markup over marginal cost; no selection on goods
- **Melitz**: Draw Pareto productivities; only firms clearing the export cost cutoff sell abroad. Trade shares computed over the "common set" of goods available everywhere.

Each simulation produces trade shares and a matrix of prices, from which simulated `d_ni` moments are computed.

### Step 4: GMM Estimation (`estimate-trade-models.jl`)

The trade elasticity θ is estimated by matching data moments to simulated moments:

- **Exactly identified** ("exact"): Match `mean(d_ni^data - d_ni^model(θ)) = 0`. Single moment, single parameter. Solved by minimizing the squared moment.

- **Over-identified** ("over"): Match three moments — `d_ni` (max price gap), `d_ni2` (85th percentile price gap), and `cov(d_ni, log(distance))` — with an optimal weighting matrix. Reports a J-statistic for the overidentifying restriction test.

Confidence intervals are constructed via parametric bootstrap: re-draw gravity residuals → re-estimate gravity → re-simulate → re-estimate θ.

## Sample Use Case

To estimate the EK model:

```julia
# From the model-code/ directory:
cd("model-code")
include("estimate-ek-model.jl")
```

This will:
1. Load price gap and trade share data for years 2004, 2011, 2017
2. Run gravity regressions to recover trade costs and technology
3. Estimate θ by matching simulated EK price gap moments to data
4. Run 100 bootstrap replications for confidence intervals
5. Write results to `results/ek-estimate-{method}-{date}.csv`

To run all models at once:
```julia
cd("model-code")
include("run-estimation-all-models.jl")
```

### Finding Results

Results are saved to `results/` as individual CSV files named `{model}-estimate-{method}-{date}.csv`. The `summary-results/combine-results.ipynb` notebook aggregates these into combined files. See [`summary-results/estimation-results-summary.md`](summary-results/estimation-results-summary.md) for a tabular overview of all estimates organized by model and year. For an analysis of the price gap data properties — including correlations with gravity variables and regression diagnostics — see the [price gap moments analysis](summary-results/pricegap-moments-analysis.md).

### Key Output Columns

| Column | Description |
|--------|-------------|
| `θ` | Estimated trade elasticity |
| `J` | J-statistic (overidentification test; 0 for exactly identified) |
| `p10`, `p90` | 10th and 90th percentile of bootstrap θ distribution |
| `Jpercentile` | Fraction of bootstrap J-stats below the point estimate J |
| `model_moment`, `data_moment` | Matched moments at the estimate |
| `model` | Trade model (ek, bejk, krugman, melitz) |
| `method` | Identification strategy (exact or over) |
| `year` | Data year (2004, 2011, 2017) |

## Dependencies

Julia packages required:
- `FixedEffectModels`, `DataFrames`, `CSV` — data handling and gravity estimation
- `Distributions`, `Random`, `StatsBase` — simulation
- `MINPACK` — nonlinear equation solving
- `Optimization`, `OptimizationOptimJL`, `OptimizationPRIMA` — GMM optimization (BOBYQA)
- `Parameters` — struct unpacking
- `HypothesisTests` — correlation tests
- `Plots` — visualization
- `SpecialFunctions`, `BenchmarkTools`, `LinearAlgebra` — utilities

Python (for data cleaning notebooks):
- `pandas`, `numpy`

## References

- Simonovska, I. and Waugh, M.E. "Trade Models, Trade Elasticities, and the Gains from Trade"
- Eaton, J. and Kortum, S. (2002). "Technology, Geography, and Trade." *Econometrica*
- Bernard, A., Eaton, J., Jensen, J.B., and Kortum, S. (2003). "Plants and Productivity in International Trade." *AER*
- Melitz, M.J. (2003). "The Impact of Trade on Intra-Industry Reallocations and Aggregate Industry Productivity." *Econometrica*
- Chaney, T. (2008). "Distorted Gravity: The Intensive and Extensive Margins of International Trade." *AER*
- Waugh, M.E. (2010). "International Trade and Income Differences." *AER*
