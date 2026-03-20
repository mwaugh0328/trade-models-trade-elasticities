# Estimation Results Summary

Summary of trade elasticity (θ) estimates across models, identification strategies, and years. Confidence intervals are 80% (p10–p90) from parametric bootstrap with 100 replications.

---

## 1. Baseline Results: Exactly Identified

Single moment condition: `E[d_ni^data - d_ni^model(θ)] = 0`.

| Model | Year | θ | 80% CI (p10–p90) | Data Moment | Model Moment |
|-------|------|-----|-------------------|-------------|--------------|
| **BEJK** | 2004 | 2.82 | 2.73 – 2.93 | 0.93 | 0.93 |
| **BEJK** | 2011 | 2.73 | 2.66 – 2.85 | 0.98 | 0.98 |
| **BEJK** | 2017 | 2.87 | 2.76 – 3.00 | 0.91 | 0.91 |
| **Melitz** | 2004 | 3.50 | 3.06 – 3.90 | 0.93 | 0.93 |
| **Melitz** | 2011 | 3.30 | 3.10 – 3.83 | 0.98 | 0.98 |
| **Melitz** | 2017 | 3.32 | 2.83 – 3.73 | 0.91 | 0.91 |
| **EK** | 2004 | 4.09 | 3.92 – 4.23 | 0.93 | 0.93 |
| **EK** | 2011 | 3.94 | 3.81 – 4.08 | 0.98 | 0.98 |
| **EK** | 2017 | 4.19 | 4.05 – 4.34 | 0.91 | 0.91 |
| **Krugman** | 2004 | 5.19 | 4.96 – 5.57 | 0.93 | 0.93 |
| **Krugman** | 2011 | 4.90 | 4.73 – 5.12 | 0.98 | 0.98 |
| **Krugman** | 2017 | 5.17 | 4.94 – 5.41 | 0.91 | 0.91 |

**Key takeaway**: The ordering θ\_BEJK < θ\_Melitz < θ\_EK < θ\_Krugman is consistent across all years. Models with richer micro-structure (Bertrand pricing, selection) compress observed price gaps and require lower θ to match the data.

---

## 2. Baseline Results: Over-Identified

Three moment conditions: d\_ni (max price gap), d\_ni2 (85th pctile price gap), cov(d\_ni, log distance). Optimal weighting matrix.

| Model | Year | θ | (p10–p90) | J-stat | J pctile | Moment 1 (model/data) | Moment 2 (model/data) | Moment 3 (model/data) |
|-------|------|-----|--------|--------|----------|----------------------|----------------------|----------------------|
| **BEJK** | 2004 | 2.76 | 2.71 – 2.83 | 23.3 | 0.68 | 0.95 / 0.93 | 0.36 / 0.37 | 0.07 / 0.09 |
| **BEJK** | 2011 | 2.80 | 2.74 – 2.88 | 173.2 | 1.00 | 0.96 / 0.98 | 0.34 / 0.34 | 0.07 / 0.25 |
| **BEJK** | 2017 | 2.88 | 2.81 – 2.95 | 180.7 | 1.00 | 0.91 / 0.91 | 0.33 / 0.32 | 0.06 / 0.20 |
| **EK** | 2004 | 4.00 | 3.91 – 4.12 | 10.7 | 0.38 | 0.95 / 0.93 | 0.36 / 0.37 | 0.09 / 0.09 |
| **EK** | 2011 | 4.07 | 3.99 – 4.18 | 136.8 | 0.99 | 0.95 / 0.98 | 0.34 / 0.34 | 0.09 / 0.25 |
| **EK** | 2017 | 4.18 | 4.09 – 4.31 | 139.3 | 1.00 | 0.91 / 0.91 | 0.32 / 0.32 | 0.07 / 0.20 |
| **Krugman** | 2004 | 3.76 | 3.56 – 3.92 | 1390.5 | 1.00 | 1.28 / 0.93 | 0.25 / 0.37 | 0.23 / 0.09 |
| **Krugman** | 2011 | 3.49 | 3.34 – 3.60 | 983.9 | 1.00 | 1.38 / 0.98 | 0.25 / 0.34 | 0.23 / 0.25 |
| **Krugman** | 2017 | 3.82 | 3.64 – 3.98 | 1411.7 | 1.00 | 1.23 / 0.91 | 0.21 / 0.32 | 0.20 / 0.20 |
| **Melitz** | 2004 | 3.59 | 3.38 – 3.91 | 133.2 | 0.86 | 0.90 / 0.93 | 0.34 / 0.37 | 0.21 / 0.09 |
| **Melitz** | 2011 | 3.15 | 2.91 – 3.67 | 60.8 | 0.53 | 1.03 / 0.98 | 0.29 / 0.34 | 0.25 / 0.25 |
| **Melitz** | 2017 | 3.08 | 2.60 – 3.82 | 287.5 | 0.92 | 0.98 / 0.91 | 0.23 / 0.32 | 0.25 / 0.20 |

**Key takeaway**: Krugman shows very large J-statistics and systematic moment mismatches (model moment 1 >> data), indicating the model is rejected by the overidentifying restrictions. EK performs best in 2004 (J=10.7). BEJK and Melitz are intermediate.

---

## 3. BEJK Sensitivity to σ (Elasticity of Substitution)

### Exactly Identified

| σ | Year | θ | (p10–p90) |
|---|------|-----|--------|
| 1.5 | 2004 | 2.82 | 2.73 – 2.93 |
| 1.5 | 2011 | 2.73 | 2.66 – 2.85 |
| 1.5 | 2017 | 2.87 | 2.76 – 3.00 |
| 2.0 | 2004 | 3.07 | 2.95 – 3.17 |
| 2.0 | 2011 | 2.99 | 2.89 – 3.09 |
| 2.0 | 2017 | 3.12 | 3.00 – 3.22 |
| 2.5 | 2004 | 3.28 | 3.13 – 3.36 |
| 2.5 | 2011 | 3.19 | 3.08 – 3.30 |
| 2.5 | 2017 | 3.33 | 3.21 – 3.43 |
| 3.0 | 2004 | 3.43 | 3.28 – 3.53 |
| 3.0 | 2011 | 3.33 | 3.22 – 3.45 |
| 3.0 | 2017 | 3.49 | 3.36 – 3.61 |

### Over-Identified

| σ | Year | θ | (p10–p90) | J-stat | J pctile |
|---|------|-----|--------|--------|----------|
| 1.5 | 2004 | 2.76 | 2.71 – 2.83 | 23.3 | 0.68 |
| 1.5 | 2011 | 2.80 | 2.74 – 2.88 | 173.2 | 1.00 |
| 1.5 | 2017 | 2.88 | 2.81 – 2.95 | 180.7 | 1.00 |
| 2.0 | 2004 | 2.96 | 2.91 – 3.03 | 61.6 | 0.97 |
| 2.0 | 2011 | 2.97 | 2.92 – 3.05 | 180.1 | 1.00 |
| 2.0 | 2017 | 3.07 | 3.00 – 3.14 | 174.3 | 1.00 |
| 2.5 | 2004 | 3.13 | 3.08 – 3.20 | 86.3 | 1.00 |
| 2.5 | 2011 | 3.13 | 3.07 – 3.21 | 192.5 | 1.00 |
| 2.5 | 2017 | 3.24 | 3.18 – 3.32 | 183.9 | 1.00 |
| 3.0 | 2004 | 3.26 | 3.21 – 3.34 | 93.3 | 0.99 |
| 3.0 | 2011 | 3.25 | 3.19 – 3.34 | 195.1 | 0.99 |
| 3.0 | 2017 | 3.37 | 3.31 – 3.45 | 189.5 | 1.00 |

**Key takeaway**: θ increases monotonically with σ. As σ rises, Bertrand markups shrink, reducing the model's ability to compress price gaps, so a higher θ is needed. By σ = 3.0, BEJK estimates approach the Melitz range.

---

## 4. Krugman Model 2 Results

Incorporates asymmetric country size (number of varieties proportional to expenditure).

### Over-Identified

| Year | θ | (p10–p90) | J-stat | J pctile | Moment 1 (model/data) |
|------|-----|--------|--------|----------|----------------------|
| 2004 | 3.72 | 3.55 – 3.93 | 122.1 | 0.89 | 0.92 / 0.93 |
| 2011 | 3.39 | 3.03 – 3.57 | 45.5 | 0.42 | 1.05 / 0.98 |
| 2017 | 3.51 | 3.23 – 4.02 | 45.4 | 0.36 | 0.95 / 0.91 |

### Exactly Identified

| Year | θ | (p10–p90) |
|------|-----|--------|
| 2004 | 3.68 | 3.28 – 4.04 |
| 2011 | 3.62 | 3.31 – 4.12 |
| 2017 | 3.65 | 3.25 – 4.06 |

**Key takeaway**: Model 2 substantially lowers θ compared to baseline Krugman (from ~5.0 to ~3.5–3.7) and dramatically improves moment fit (J drops from ~1000 to ~50-120). Accounting for country size effects matters significantly for the Krugman framework.

---

## 5. Melitz Model 2 Results (Sensitivity to σ)

### Exactly Identified

| σ | Year | θ | (p10–p90) |
|---|------|-----|--------|
| 1.5 | 2004 | 1.65 | 1.51 – 1.80 |
| 1.5 | 2011 | 1.50 | 1.38 – 1.58 |
| 1.5 | 2017 | 1.50 | 1.35 – 1.78 |
| 2.0 | 2004 | 2.12 | 1.98 – 2.33 |
| 2.0 | 2011 | 2.05 | 1.90 – 2.18 |
| 2.0 | 2017 | 2.02 | 1.88 – 2.16 |
| 2.5 | 2004 | 2.47 | 2.32 – 2.73 |
| 2.5 | 2011 | 2.46 | 2.32 – 2.61 |
| 2.5 | 2017 | 2.45 | 2.28 – 2.61 |
| 3.0 | 2004 | 2.77 | 2.53 – 3.02 |
| 3.0 | 2011 | 2.79 | 2.57 – 2.98 |
| 3.0 | 2017 | 2.82 | 2.60 – 3.01 |

### Over-Identified

| σ | Year | θ | (p10–p90) | J-stat | J pctile |
|---|------|-----|--------|--------|----------|
| 2.0 | 2004 | 2.05 | 1.93 – 2.51 | 149.5 | 0.49 |
| 2.0 | 2011 | 2.03 | 1.96 – 2.50 | 12.3 | 0.10 |
| 2.0 | 2017 | 1.99 | 1.92 – 9.20 | 48.0 | 0.15 |
| 2.5 | 2004 | 2.55 | 2.40 – 2.98 | 156.1 | 0.71 |
| 2.5 | 2011 | 2.50 | 2.33 – 2.76 | 40.2 | 0.33 |
| 2.5 | 2017 | 2.52 | 2.34 – 2.85 | 26.8 | 0.29 |
| 3.0 | 2004 | 2.84 | 2.71 – 3.19 | 145.4 | 0.77 |
| 3.0 | 2011 | 2.93 | 2.75 – 3.24 | 38.2 | 0.31 |
| 3.0 | 2017 | 2.92 | 2.72 – 3.18 | 44.7 | 0.49 |

**Key takeaway**: Melitz Model 2 produces the lowest θ estimates of any specification, especially at low σ. Like BEJK, θ increases with σ. The over-identified J-statistics are moderate, suggesting reasonable fit, particularly for 2011 and 2017 at σ ≥ 2.0.

---

## Summary Across Models (Exactly Identified, Baseline σ)

| Model | 2004 | 2011 | 2017 | Avg |
|-------|------|------|------|-----|
| BEJK (σ=1.5) | 2.82 | 2.73 | 2.87 | 2.81 |
| Melitz | 3.50 | 3.30 | 3.32 | 3.37 |
| EK | 4.09 | 3.94 | 4.19 | 4.07 |
| Krugman | 5.19 | 4.90 | 5.17 | 5.09 |
