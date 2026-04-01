# Model Descriptions

This document describes the economic environment and simulation structure for each trade model estimated in this project. All models share a common empirical framework — a gravity equation estimated on bilateral trade data for 30 countries — but differ in the micro-foundations that generate trade patterns and prices.

---

## Common Environment

### Data and Parameters

All models take as inputs:

- **$S_i$** (technology): Country-specific productivity parameters recovered from the gravity regression's exporter fixed effects. Formally, $S_i = \exp(\hat{\alpha}_i^{\mathrm{exporter}})$ where $\hat{\alpha}_i$ is the estimated exporter fixed effect.
- **$d_{ni}$** (trade costs): Bilateral iceberg trade costs constructed from the gravity regression. For country $n$ importing from $i$:

$$d_{ni} = \exp\left(-\frac{1}{\theta}\left[\hat{\delta}(\mathrm{dist}_{ni}) + \hat{\gamma} \cdot \mathrm{border}_{ni} + \hat{\theta}_{m,n}\right]\right)$$

where $\hat{\delta}(\cdot)$ are distance-bin coefficients, $\hat{\gamma}$ is the border effect, and $\hat{\theta}_{m,n}$ is an importer-specific asymmetric trade cost (following Eaton–Kortum's normalization). Trade costs are floored at 1.

- **$\theta$** (trade elasticity): The parameter being estimated. It governs the dispersion of productivity draws (Fréchet shape in EK/BEJK, Pareto shape in Melitz) or the demand elasticity (Krugman).
- **$\sigma$** (elasticity of substitution): Relevant for BEJK and Melitz models. Baseline value is 1.5; sensitivity analysis explores $\sigma \in \lbrace 1.5, 2.0, 2.5, 3.0 \rbrace$.

### Gravity Equation

The gravity regression takes the trade-share form from Eaton and Kortum (2002):

$$\log\left(\frac{\pi_{ni}}{\pi_{nn}}\right) = \alpha_i^{\mathrm{exp}} - \alpha_n^{\mathrm{imp}} + \delta(\mathrm{dist}_{ni}) + \gamma \cdot \mathrm{border}_{ni} + \varepsilon_{ni}$$

with six distance bins (0–375, 375–750, 750–1500, 1500–3000, 3000–6000, 6000+ miles). The exporter fixed effects yield $S_i$, the importer fixed effects yield asymmetric trade costs $\theta_{m,n}$, and the bilateral variables determine $d_{ni}$.

### Simulation Approach

Each model simulates a large number of goods ( $N = 100,\!000$ – $150,\!000$ ) and computes:
1. **Trade shares** $\pi_{ni}$: fraction of country $n$'s expenditure on goods from country $i$.
2. **Prices**: good-level prices observed in each country, used to construct the price gap moments $d_{ni}$ and $d_{ni2}$.

The estimation matches simulated moments (price gaps) to data moments via GMM.

---

## 1. Eaton–Kortum (EK)

### Economic Environment

- **Market structure**: Perfect competition. Each good is produced by a single country — whichever has the lowest delivered cost.
- **Productivity**: Country $i$'s unit cost for good $g$ is drawn from a Fréchet distribution: $\Pr(z_{ig} \le z) = \exp(-S_i \cdot z^{-\theta})$. Higher $S_i$ means country $i$ is more likely to draw low costs.
- **Trade**: Country $n$ buys good $g$ from whichever source $i$ offers the lowest CIF price $p_{nig} = d_{ni} \cdot c_{ig}$. There is no markup — prices equal marginal cost.
- **Price aggregation**: The CES price index uses $p^{1-\sigma}$ weighting, but since pricing is competitive, $\sigma$ only affects the price index, not the trade pattern.

### Simulation

1. Draw $N$ independent uniform random variables per country; transform via $c = (\log u / (-S_i))^{1/\theta}$ to get Fréchet-distributed marginal costs.
2. For each good, find the minimum CIF price across all exporters: $p_{ng} = \min_{i} \lbrace d_{ni} \cdot c_{ig} \rbrace$.
3. Assign trade to the lowest-cost supplier; aggregate using CES weights $p^{1-\sigma}$.
4. All goods exist in all countries → no common-set selection needed.

### Key Properties

- Every good is supplied by exactly one country (winner-take-all).
- All goods are available everywhere, so the "common set" for price gap computation is the full set.
- The trade elasticity is $\theta$ — a higher $\theta$ means less productivity dispersion, more equal unit costs, and smaller price gaps.

---

## 2. Bernard–Eaton–Jensen–Kortum (BEJK)

### Economic Environment

- **Market structure**: Bertrand competition. In each good, the lowest-cost producer sets price to just undercut the second-lowest competitor, subject to a CES markup ceiling.
- **Productivity**: Same Fréchet draws as EK. Each country draws *two* independent productivities per good: a first-best and a second-best. The second-best draw is conditional on being worse than the first: $c_2 = (\log u_2 / (-S_i) + c_1^{\theta})^{1/\theta}$.
- **Pricing rule**: The winning firm charges:

$$p = \min\left(d_{ni^{*}} \cdot c_{2,i^{*}},\ \min_{j \neq i^{*}} d_{n,j} \cdot c_{1,j},\ \frac{\sigma}{\sigma-1} \cdot d_{ni^{*}} \cdot c_{1,i^{*}}\right)$$

where $i^*$ is the lowest-cost international supplier. The price is the minimum of: (a) the second-best domestic draw (times trade cost), (b) the best foreign competitor's delivered cost, or (c) the monopoly markup.

- **Role of $\sigma$**: The CES markup $\sigma/(\sigma-1)$ caps the price. With low $\sigma$ (e.g., 1.5), markups are large (3×), giving the winner room to charge above marginal cost. As $\sigma \to \infty$, Bertrand pricing collapses to competitive pricing and BEJK → EK.

### Simulation

1. Draw first-best marginal costs from Fréchet (same as EK, same seed).
2. Draw second-best marginal costs conditional on being worse than first-best.
3. For each good, find the two lowest international delivered costs.
4. Price = min(second-best cost, CES markup × lowest cost).
5. Aggregate trade shares with CES weights.

### Key Properties

- Variable markups compress observed prices relative to marginal costs. For a given $\theta$, price gaps are smaller than under EK because the Bertrand mechanism limits how far prices can diverge.
- To match the same data moments, BEJK requires a lower $\theta$ than EK.
- Increasing $\sigma$ shrinks the markup ceiling, reducing the compression, and BEJK estimates approach EK.

---

## 3. Krugman

### Economic Environment

- **Market structure**: Monopolistic competition with CES demand. Each firm produces a unique variety and charges a constant markup over marginal cost.
- **Productivity**: Unlike EK/BEJK, there are no productivity draws across goods within a country. All firms in country $i$ have the same unit cost $c_i = (\log(0.5)/(-S_i))^{1/\theta}$ (a degenerate draw at the median).
- **Number of varieties**: Each country produces $N$ varieties (symmetric across countries). All varieties from all countries are consumed everywhere.
- **Demand elasticity**: $\eta = \theta + 1$. The markup is $\eta/(\eta-1) = (\theta+1)/\theta$.
- **Pricing**: The delivered price of any variety from $i$ in market $n$ is $p_{ni} = d_{ni} \cdot \eta/(\eta-1) \cdot c_i$.

### Simulation

1. All firms in country $i$ produce at the same cost (no heterogeneity within country).
2. Every variety is available in every market at price = trade cost × markup × unit cost.
3. Trade shares are CES aggregates: $\pi_{ni} \propto N \cdot (d_{ni} \cdot c_i)^{1-\eta}$.
4. Prices are sampled uniformly across the $N_{\mathrm{cntry}} \times N_{\mathrm{goods}}$ varieties.

### Key Properties

- No within-country heterogeneity → all price variation comes from *between*-country cost differences amplified by trade costs.
- No extensive margin: all varieties are traded everywhere.
- The demand elasticity $\eta = \theta + 1$ directly governs trade flows, so $\theta$ plays a different role than in EK/BEJK.
- Produces the highest $\theta$ estimates because the model has the least internal mechanism to generate price variation — all of it must come from $\theta$.

---

## 4. Melitz

### Economic Environment

- **Market structure**: Monopolistic competition with heterogeneous firms and export selection. Firms must pay a fixed cost to export; only the most productive survive.
- **Productivity**: Firm productivity $\varphi$ follows a Pareto distribution with shape parameter $\theta$: $\Pr(\varphi \ge x) = (x/\varphi^{*})^{-\theta}$, conditional on entry. Higher $\theta$ means less dispersion (more firms near the cutoff).
- **Entry and selection**: A firm from country $i$ with marginal cost $c$ can sell in market $n$ only if its delivered price $d_{ni} \cdot p_i(c)$ clears market $n$'s cost cutoff $\bar{c}_n$:

$$d_{ni} \cdot \frac{\sigma}{\sigma-1} \cdot c \leq \frac{\sigma}{\sigma-1} \cdot \bar{c}_n$$

Under the Simonovska–Waugh (SW) / Chaney (2008) assumption that fixed export costs are proportional to origin country size ($f_{ij} = f \cdot L_i$), the cutoff simplifies:

$$\bar{c}_n = \Phi_n^{-1/\theta}, \qquad \Phi_n = \sum_k S_k \cdot d_{nk}^{-\theta}$$

- **Domestic share**: $\pi_{nn} = S_n / \Phi_n$ — the fraction of varieties in $n$ produced domestically.
- **Markup**: CES markup $\sigma/(\sigma-1)$, same for all firms (no variable markups).

### Simulation

1. For each country $j$, determine the number of active domestic varieties: proportional to $\pi_{jj}$.
2. Draw uniform random variables, scale by $\pi_{jj}$, and transform: $p = \frac{\sigma}{\sigma-1} \cdot (u \cdot \pi_{jj} / S_j)^{1/\theta}$ to get Pareto-distributed prices.
3. For each good, check whether foreign firms clear the destination's cost cutoff: $d_{ni} \cdot p_{ig} \leq \frac{\sigma}{\sigma-1} \cdot \bar{c}_n$.
4. Aggregate trade shares using CES weights $p^{1-\sigma}$.
5. Identify the **common set**: goods sold in *all* countries (only the most productive firms export everywhere). Sample from this common set for price gap computation.

### Key Properties

- **Export selection** compresses the observed distribution of prices. Only the most productive (lowest-price) firms export to all markets, so the common set used for price gap moments is truncated — right tail prices are missing.
- $\sigma$ affects both the markup level and the selection cutoff.

---

## 5. Krugman Model 2 — Asymmetric Country Size

### Difference from Baseline Krugman

Baseline Krugman assumes each country produces the same number of varieties $N$. Krugman Model 2 relaxes this: **the number of varieties produced by country $i$ is proportional to its expenditure $E_i$**.

This changes two things:

1. **Trade shares**: Expenditure-weighted CES aggregation. The contribution of exporter $i$ to importer $n$'s price index is weighted by $E_i$:

$$\pi_{ni} \propto E_i \cdot (d_{ni} \cdot c_i)^{1-\eta}$$

Larger countries produce more varieties and account for a larger share of the import bundle.

2. **Price sampling**: When sampling prices for moment computation, varieties are drawn with probability proportional to $E_i$ (the exporter's expenditure), rather than uniformly. This reflects the fact that a random good is more likely to originate from a large economy.

### Effect on Estimates

Accounting for asymmetric country size substantially lowers $\theta$ estimates relative to baseline Krugman (from $\approx 5.0$ to $\approx 3.5$) and dramatically improves model fit (J-statistics drop from $\sim 1000$ to $\sim 50$–$120$). The baseline Krugman's symmetric-variety assumption forces $\theta$ to do too much work — it must explain size-driven variation that is better captured by $E_i$.

---

## Summary: Model Ordering

The models form a consistent ordering in $\theta$ estimates:

$$\theta_{\mathrm{BEJK}} < \theta_{\mathrm{Melitz}} < \theta_{\mathrm{EK}} < \theta_{\mathrm{Krugman}}$$

This ordering reflects the richness of each model's micro-structure:
- **BEJK** has Bertrand variable markups that compress prices → needs the lowest $\theta$.
- **Melitz** has export selection that truncates the price distribution → needs lower $\theta$ than EK.
- **EK** is perfectly competitive with no selection → moderate $\theta$.
- **Krugman** has no within-country heterogeneity and no selection → needs the highest $\theta$ to generate sufficient price variation.

The Krugman Model 2 variant incorporates country size via expenditure $E_i$, producing lower $\theta$ than baseline Krugman because expenditure-driven variation substitutes for productivity dispersion in explaining the data.
