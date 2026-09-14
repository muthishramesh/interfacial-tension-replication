"""
Replication of all quantitative figures in:
Kim, H. & Burgess, D. J. (2001),
"Prediction of Interfacial Tension between Oil Mixtures and Water",
Journal of Colloid and Interface Science 241, 509-513.
doi:10.1006/jcis.2001.7655

Figures reproduced:
Fig. 2 - IFT vs time for mineral oil/n-butyl acetate
Fig. 3 - five oil-mixture systems (experimental points + model curves)
Fig. 4 - alpha vs delta-gamma
Fig. 5 - mineral oil/ethyl caprate
Fig. 6 - mineral oil/octanol

IMPORTANT:
The article provides plotted experimental points, but not the underlying
numeric raw datasets. Therefore the CSV experimental point files are
approximate digitizations from the published raster figures. The analytical
model curves use the equations and parameter values printed in the paper.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)

# -----------------------------
# Model equations from the paper
# -----------------------------
def ift_model(V, gamma1, gamma2, alpha):
    """Eq. (5): gamma = (gamma1-gamma2)*exp(-alpha*V^0.7)+gamma2."""
    V = np.asarray(V, dtype=float)
    return (gamma1 - gamma2) * np.exp(-alpha * V**0.7) + gamma2

def alpha_from_delta_gamma(delta_gamma):
    """Eq. (6): alpha = 0.178*Delta_gamma - 0.0132."""
    return 0.178 * np.asarray(delta_gamma) - 0.0132


PURE = {
    "Mineral oil": 49.0,
    "Benzene": 32.7,
    "Isopropyl myristate": 24.2,
    "n-Butyl acetate": 13.0,
    "Benzyl alcohol": 3.3,
}

SYSTEMS = {
    "a": {
        "title": "(a) Mineral oil + n-butyl acetate",
        "file": "fig3_a_mineral_nBA.csv",
        "gamma1": 49.0, "gamma2": 10.3,
        "alpha_exp": 6.8, "alpha_calc": 6.9,
        "xlabel": "Volume Fraction of n-Butylacetate in Mineral Oil",
    },
    "b": {
        "title": "(b) Benzene + benzyl alcohol",
        "file": "fig3_b_benzene_benzyl_alcohol.csv",
        "gamma1": 29.5, "gamma2": 3.3,
        "alpha_exp": 5.6, "alpha_calc": 5.2,
        "xlabel": "Volume Fraction of Benzyl Alcohol in Benzene",
    },
    "c": {
        "title": "(c) Mineral oil + isopropyl myristate",
        "file": "fig3_c_mineral_ipm.csv",
        "gamma1": 24.8, "gamma2": 24.2,
        "alpha_exp": 4.0, "alpha_calc": 4.4,
        "xlabel": "Volume Fraction of Isopropylmyristate in Mineral Oil",
    },
    "d": {
        "title": "(d) Mineral oil + benzene",
        "file": "fig3_d_mineral_benzene.csv",
        "gamma1": 17.5, "gamma2": 32.7,
        "alpha_exp": 3.3, "alpha_calc": 3.1,
        "xlabel": "Volume Fraction of Benzene in Mineral Oil",
    },
    "e": {
        "title": "(e) Isopropyl myristate + n-butyl acetate",
        "file": "fig3_e_ipm_nBA.csv",
        "gamma1": 16.8, "gamma2": 13.0,
        "alpha_exp": 2.9, "alpha_calc": 3.0,
        "xlabel": "Volume Fraction of n-Butylacetate in Isopropylmyristate",
    },
}

def style_axes(ax):
    ax.tick_params(direction="in", top=True, right=True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)

# -----------------------------
# Fig. 2
# -----------------------------
def plot_fig2():
    df = pd.read_csv(DATA / "fig2_time_series_digitized_approx.csv")
    fig, ax = plt.subplots(figsize=(7, 4.8))
    for vf, g in df.groupby("volume_fraction"):
        ax.plot(
            g["time_min"], g["interfacial_tension_mN_m"],
            marker="o", markersize=3.5, linewidth=1.0,
            label=f"Volume Fraction = {vf:g}"
        )
    ax.set_xlim(0, 40)
    ax.set_ylim(0, 60)
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Interfacial tension (mN/m)")
    ax.set_title("Fig. 2 replication — approximate digitization")
    ax.legend(fontsize=7, loc="center right", frameon=True)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(OUT / "fig2_replicated.png", dpi=300)
    fig.savefig(OUT / "fig2_replicated.svg")
    plt.close(fig)

# -----------------------------
# Fig. 3
# -----------------------------
def plot_fig3():
    fig, axes = plt.subplots(3, 2, figsize=(10, 12))
    axes = axes.ravel()

    for i, (key, s) in enumerate(SYSTEMS.items()):
        ax = axes[i]
        df = pd.read_csv(DATA / s["file"])
        V = np.linspace(0, 1, 500)

        y_exp = ift_model(V, s["gamma1"], s["gamma2"], s["alpha_exp"])
        y_calc = ift_model(V, s["gamma1"], s["gamma2"], s["alpha_calc"])

        ax.plot(V, y_exp, linewidth=1.4, label=rf"solid: $\alpha={s['alpha_exp']}$")
        ax.plot(V, y_calc, "--", linewidth=1.1, label=rf"dashed: $\alpha={s['alpha_calc']}$")
        ax.scatter(
            df["volume_fraction"], df["experimental_gamma"],
            s=24, facecolors="black", edgecolors="black",
            zorder=3, label="digitized experimental points"
        )

        ax.set_xlim(0, 1)
        ax.set_ylim(0, 50)
        ax.set_xlabel(s["xlabel"], fontsize=8)
        ax.set_ylabel("Interfacial Tension (mN/m)")
        ax.set_title(s["title"], loc="left", fontsize=10, fontweight="bold")
        ax.legend(fontsize=7, frameon=False)
        style_axes(ax)

    axes[-1].axis("off")
    fig.suptitle(
        "Fig. 3 replication — model curves from Eq. (5) and approximate digitization",
        y=0.995, fontsize=12
    )
    fig.tight_layout()
    fig.savefig(OUT / "fig3_replicated.png", dpi=300)
    fig.savefig(OUT / "fig3_replicated.svg")
    plt.close(fig)

# -----------------------------
# Fig. 4
# -----------------------------
def plot_fig4():
    df = pd.read_csv(DATA / "fig4_alpha_delta_gamma.csv")

    # Eq. (6), used for the straight regression/prediction line.
    x = np.linspace(15, 40, 300)
    y = alpha_from_delta_gamma(x)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(
        df["delta_gamma_mN_m"], df["alpha_experimental"],
        s=35, facecolors="black", edgecolors="black",
        label="Experimental α"
    )
    ax.plot(x, y, linewidth=1.2, label=r"Eq. (6): $\alpha=0.178\Delta\gamma-0.0132$")

    ax.set_xlim(15, 40)
    ax.set_ylim(2, 8)
    ax.set_xlabel("Delta Interfacial Tension (mN/m)")
    ax.set_ylabel(r"$\alpha$")
    ax.set_title(r"Fig. 4 replication — $\alpha$ vs. $\Delta\gamma$")
    ax.legend(frameon=False)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(OUT / "fig4_replicated.png", dpi=300)
    fig.savefig(OUT / "fig4_replicated.svg")
    plt.close(fig)

# -----------------------------
# Fig. 5
# -----------------------------
def plot_fig5():
    df = pd.read_csv(DATA / "fig5_mineral_ethyl_caprate.csv")
    V = np.linspace(0, 1, 500)

    # Paper: pure ethyl caprate = 23.6 mN/m and alpha = 4.7.
    y = ift_model(V, 49.0, 23.6, 4.7)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(V, y, linewidth=1.3, label=r"Paper prediction: $\gamma=26.4e^{-4.7V^{0.7}}+23.6$")
    ax.scatter(
        df["volume_fraction"], df["experimental_gamma"],
        s=45, facecolors="white", edgecolors="black",
        label="Digitized experimental points", zorder=3
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 50)
    ax.set_xlabel("Volume Fraction of Ethyl Caprate in Mineral Oil")
    ax.set_ylabel("Interfacial Tension (mN/m)")
    ax.set_title("Fig. 5 replication — mineral oil + ethyl caprate")
    ax.legend(frameon=False, fontsize=8)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(OUT / "fig5_replicated.png", dpi=300)
    fig.savefig(OUT / "fig5_replicated.svg")
    plt.close(fig)

# -----------------------------
# Fig. 6
# -----------------------------
def plot_fig6():
    df = pd.read_csv(DATA / "fig6_mineral_octanol.csv")
    V = np.linspace(0, 1, 500)

    # Octanol pure IFT is read approximately from the plotted endpoint.
    # Eq. (6) gives alpha from Delta-gamma = 49 - 4.0.
    gamma_octanol = 4.0
    alpha = float(alpha_from_delta_gamma(49.0 - gamma_octanol))
    y = ift_model(V, 49.0, gamma_octanol, alpha)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(V, y, linewidth=1.3, label=rf"Eq. (5)+(6), $\alpha={alpha:.2f}$")
    ax.scatter(
        df["volume_fraction"], df["experimental_gamma"],
        s=45, facecolors="white", edgecolors="black",
        label="Digitized experimental points", zorder=3
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 50)
    ax.set_xlabel("Volume Fraction of Octanol in Mineral Oil")
    ax.set_ylabel("Interfacial Tension (mN/m)")
    ax.set_title("Fig. 6 replication — mineral oil + octanol")
    ax.legend(frameon=False, fontsize=8)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(OUT / "fig6_replicated.png", dpi=300)
    fig.savefig(OUT / "fig6_replicated.svg")
    plt.close(fig)

if __name__ == "__main__":
    plot_fig2()
    plot_fig3()
    plot_fig4()
    plot_fig5()
    plot_fig6()
    print(f"Figures written to: {OUT}")
