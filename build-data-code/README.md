# Data Construction Pipeline

This directory contains all the code and raw data inputs used to build the processed datasets in `data/` that feed into the trade elasticity estimation.

## Output Files

The pipeline produces the following files in `data/`:

| File | Description |
|------|-------------|
| `pricegap-df-{year}.csv` | Bilateral price gap moments (`d_ni`, `d_ni2`, `d_ni3`) and normalized trade shares for each country pair |
| `tradeshare-df-{year}.csv` | Bilateral trade share matrix merged with gravity variables |
| `tariffs-{year}.csv` | Bilateral applied tariff rates |
| `top30_gravity_data.csv` | Gravity covariates (distance bins, border, language) for the 30-country sample |

## Driver Scripts

| Script | Years | Description |
|--------|-------|-------------|
| `make_price_data.jl` | 2004, 2011 | Main driver. Loads ICP price data and Feenstra trade flows, constructs trade shares and price gap moments, writes `pricegap-df-{year}.csv` and `tradeshare-df-{year}.csv` |
| `make2017_data.jl` | 2017 | Same as above but uses BACI trade flow data (cleaned by `clean-trade-data.ipynb`) and pre-processed 2017 ICP prices |
| `make_trade_share.jl` | 2004, 2011 | Standalone script to build only the trade share matrices (alternative to `make_price_data.jl`) |

## Shared Functions

`functions_make_trade_share.jl` contains the core functions used by all driver scripts:

- **`construct_tradematrix()`** — Aggregates raw bilateral trade flows (by ISIC code) into a country-by-country trade matrix
- **`aggregate_drop()`** — Merges country aggregates (Belgium+Netherlands, China+HK+Macau, Singapore+Malaysia) and selects the top 30 countries
- **`construct_tradeshare()`** — Converts trade flows into expenditure shares: `π_ni = X_ni / (GO_n + M_n - E_n)`
- **`construct_tradeshare_baci()`** — Same as above but for the BACI trade data format (used for 2017)
- **`adjust_price_data()`** — Processes raw ICP PPP price data: drops non-traded categories, adjusts by exchange rates, selects the 30-country sample
- **`build_dni()`** — Computes price gap moments from the price matrix: `d_ni = max(log p_n - log p_i) - mean(log p_n - log p_i)` and the 85th/75th percentile variants

## Raw Data Sources

### ICP Price Data

| Folder | Year | Format | Notes |
|--------|------|--------|-------|
| `2004-price-data/` | 2004 | `.mat` (MATLAB) | Pre-processed price matrix from earlier version of the project |
| `2011-price-data/` | 2011 | `.xlsx` | ICP basic headings with PPP price levels |
| `2017-price-gap/` | 2017 | `.csv`, `.xlsx` | ICP 2017 round data; cleaned via `ICP_micro_price_data_PPP.ipynb` |

### Trade Flow Data

| Folder | Years | Source | Notes |
|--------|-------|--------|-------|
| `feenstra-trade-data/` | 2004, 2011 | Feenstra World Trade Flows | Bilateral flows by ISIC 3-digit manufacturing codes. Note: this data source may no longer be publicly available |
| `baci-trade-data/` | 2017 | BACI (CEPII) | Cleaned and aggregated via `clean-trade-data.ipynb` |

### Gross Output Data

| Folder | Years | Source |
|--------|-------|--------|
| `UNIDO-data/` | 2004, 2011, 2017 | UNIDO INDSTAT manufacturing gross output. Cleaned via `clean-grossoutput-data.ipynb` |

### Gravity Variables

| Folder | Description |
|--------|-------------|
| `make-gravity-var/` | Constructs bilateral gravity variables (distance, border, shared language, EU/EFTA membership). Uses Stata `.do` files and a Python notebook (`make-gravity-var.ipynb`). Output: `top30_gravity_data.csv` |

## Data Cleaning Notebooks

| Notebook | Language | Purpose |
|----------|----------|---------|
| `clean-grossoutput-data.ipynb` | Python | Cleans UNIDO gross output data for 2004, 2011, 2017 |
| `clean-trade-data.ipynb` | Python | Cleans and aggregates BACI trade data for 2017 |
| `2017-price-gap/ICP_micro_price_data_PPP.ipynb` | Python | Processes raw ICP 2017 price data into traded-goods price matrix |

## Country Sample

The file `drop_30.csv` defines which countries are included in the 30-country sample. Country aggregation rules (hardcoded in `aggregate_drop()`):
- Belgium + Netherlands
- China + Hong Kong + Macau
- Singapore + Malaysia

## How to Rebuild the Data

From the repository root:

```julia
# Build 2004 and 2011 data
include("build-data-code/make_price_data.jl")

# Build 2017 data
include("build-data-code/make2017_data.jl")
```

The Python notebooks should be run first if the raw data has changed. Gravity variables are pre-built in `make-gravity-var/top30_gravity_data.csv` and copied to `data/`.



