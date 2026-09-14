# Replication of Kim & Burgess (2001) — Oil Mixture/Water Interfacial Tension

This repository reproduces the quantitative plots in:

**H. Kim and D. J. Burgess (2001), "Prediction of Interfacial Tension between Oil Mixtures and Water", Journal of Colloid and Interface Science, 241, 509–513. DOI: 10.1006/jcis.2001.7655.**

## Figures included

- **Figure 2:** Interfacial tension vs time for mineral oil/n-butyl acetate at several volume fractions.
- **Figure 3:** Five oil-mixture systems, with experimental points, nonlinear-regression curves, and curves calculated from Eq. (6).
- **Figure 4:** alpha vs difference in pure-oil interfacial tension.
- **Figure 5:** Mineral oil + ethyl caprate.
- **Figure 6:** Mineral oil + octanol.

Figure 1 is an experimental apparatus schematic, not a quantitative graph, so it is not recreated as a numerical plot.

## Important data note

The article prints experimental measurements as graph markers but does **not** provide the underlying numerical raw dataset. Therefore:

1. The **model curves** are generated from the equations and parameter values printed in the article.
2. The experimental point CSVs are **approximate digitizations** of the published figures and should not be treated as the authors' raw data.
3. Table 1 and Table 2 values are transcribed from the article.

## Main equations

### Eq. (5)

```text
gamma = (gamma1 - gamma2) * exp(-alpha * V^0.7) + gamma2
```

where:
- `gamma1`, `gamma2` = pure oil/water interfacial tensions
- `alpha` = exponential coefficient
- `V` = volume fraction of the lower-IFT oil

### Eq. (6)

```text
alpha = 0.178 * Delta_gamma - 0.0132
```

## Run

```bash
pip install -r requirements.txt
python src/replicate_all_figures.py
```

The PNG and SVG files are written to `figures/`.

## Data files

- `table1_pure_oils.csv`
- `fig2_time_series_digitized_approx.csv`
- `fig3_a_mineral_nBA.csv`
- `fig3_b_benzene_benzyl_alcohol.csv`
- `fig3_c_mineral_ipm.csv`
- `fig3_d_mineral_benzene.csv`
- `fig3_e_ipm_nBA.csv`
- `fig4_alpha_delta_gamma.csv`
- `fig5_mineral_ethyl_caprate.csv`
- `fig6_mineral_octanol.csv`

## Reproducibility

The script uses only NumPy, pandas and Matplotlib. No proprietary regression software is required for the model curves because the paper reports the fitted alpha values in Table 2 and Figure 3.

For a stricter replication study, replace the approximate digitized CSV values with data obtained using a calibrated plot-digitization workflow such as WebPlotDigitizer, then rerun the script.
