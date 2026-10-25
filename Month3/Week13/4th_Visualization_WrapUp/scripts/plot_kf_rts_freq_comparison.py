"""
plot_kf_rts_freq_comparison.py

Month 3 Milestone — KF vs RTS, Single-Frequency vs Dual-Frequency
Produces a 2x2 grid:
    Top-left:     Single-freq   3D error vs time    (KF vs RTS)
    Top-right:    Dual-freq     3D error vs time    (KF vs RTS)
    Bottom-left:  Single-freq   Position 1-sigma vs time (KF vs RTS)
    Bottom-right: Dual-freq     Position 1-sigma vs time (KF vs RTS)

Both input CSVs are expected to have the same column structure produced by
kf_spp_processor.cpp (Week 12 / Week 13 Day 3):
    t_gps, kf_lat, kf_lon, kf_h, rts_lat, rts_lon, rts_h, kf_sigma, rts_sigma

Edit SINGLE_FREQ_CSV / DUAL_FREQ_CSV below to point at your actual files:
  - single-freq: Month3/Week12/2nd_RTS_Implementation/data/trajectory_rts.csv
  - dual-freq:   Month3/Week13/3rd_DualFreq_KF_RTS/data/trajectory_rts.csv
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ─────────────────────────────────────────────────────────────
# CONFIG — point these at your real files
# ─────────────────────────────────────────────────────────────
SINGLE_FREQ_CSV = "../../../Week12/2nd_RTS_Implementation/data/trajectory_rts_crctd.csv"
DUAL_FREQ_CSV   = "../data/trajectory_rts_crctd.csv"

#TRUE_LAT, TRUE_LON, TRUE_H = 33.578381, -117.064265, 472.0
TRUE_LAT = 33.578243000    # deg
TRUE_LON = -117.064599000  # deg
TRUE_H   = 470.061         # m, ellipsoidal
R_EARTH = 6378137.0
LAT_RAD = np.radians(TRUE_LAT)


def err3d(lat, lon, h):
    """3D position error (m) vs the known BILL00USA truth coordinate."""
    e = np.radians(lon - TRUE_LON) * R_EARTH * np.cos(LAT_RAD)
    n = np.radians(lat - TRUE_LAT) * R_EARTH
    dz = h - TRUE_H
    return np.sqrt(e**2 + n**2 + dz**2)


def load_and_process(csv_path):
    """Load a trajectory_rts.csv and compute KF/RTS error + hours column."""
    df = pd.read_csv(csv_path)
    df['t_hours']  = df['t_gps'] / 3600.0
    df['kf_err']   = err3d(df['kf_lat'],  df['kf_lon'],  df['kf_h'])
    df['rts_err']  = err3d(df['rts_lat'], df['rts_lon'], df['rts_h'])
    return df


# ─────────────────────────────────────────────────────────────
# LOAD DATA
# ─────────────────────────────────────────────────────────────
sf = load_and_process(SINGLE_FREQ_CSV)
df_ = load_and_process(DUAL_FREQ_CSV)   # 'df' is a common pandas alias; renamed to avoid clash

print(f"Single-freq: KF mean={sf['kf_err'].mean():.2f}m  RTS mean={sf['rts_err'].mean():.2f}m")
print(f"Dual-freq:   KF mean={df_['kf_err'].mean():.2f}m  RTS mean={df_['rts_err'].mean():.2f}m")

# ─────────────────────────────────────────────────────────────
# FIGURE — 2x2 grid
# ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(15, 10), facecolor='white')

KF_COLOR  = '#2E75B6'
RTS_COLOR = '#1E7145'

# ---------- TOP ROW: 3D error vs time ----------
ax = axes[0, 0]
ax.plot(sf['t_hours'], sf['kf_err'],  color=KF_COLOR,  lw=1.0, alpha=0.85, label='KF')
ax.plot(sf['t_hours'], sf['rts_err'], color=RTS_COLOR, lw=1.6, label='RTS')
ax.set_title(f"Single-Frequency — 3D Error vs Time\n"
             f"KF mean={sf['kf_err'].mean():.1f}m   RTS mean={sf['rts_err'].mean():.1f}m",
             fontsize=11.5, fontweight='bold', color='#1B3A6B')
ax.set_xlabel('Time (hours)'); ax.set_ylabel('3D error (m)')
ax.legend(fontsize=9.5); ax.grid(alpha=0.3); ax.set_xlim(0, 24)

ax = axes[0, 1]
ax.plot(df_['t_hours'], df_['kf_err'],  color=KF_COLOR,  lw=1.0, alpha=0.85, label='KF')
ax.plot(df_['t_hours'], df_['rts_err'], color=RTS_COLOR, lw=1.6, label='RTS')
ax.set_title(f"Dual-Frequency — 3D Error vs Time\n"
             f"KF mean={df_['kf_err'].mean():.1f}m   RTS mean={df_['rts_err'].mean():.1f}m",
             fontsize=11.5, fontweight='bold', color='#1B3A6B')
ax.set_xlabel('Time (hours)'); ax.set_ylabel('3D error (m)')
ax.legend(fontsize=9.5); ax.grid(alpha=0.3); ax.set_xlim(0, 24)

# Match y-axis scale across the top row for a fair visual comparison
top_ymax = max(sf['kf_err'].max(), df_['kf_err'].max()) * 1.05
axes[0, 0].set_ylim(0, top_ymax)
axes[0, 1].set_ylim(0, top_ymax)

# ---------- BOTTOM ROW: Position 1-sigma (uncertainty) vs time ----------
ax = axes[1, 0]
ax.plot(sf['t_hours'], sf['kf_sigma'],  color=KF_COLOR,  lw=1.2, label='KF 1-sigma')
ax.plot(sf['t_hours'], sf['rts_sigma'], color=RTS_COLOR, lw=1.6, label='RTS 1-sigma')
ax.set_title("Single-Frequency — Position Uncertainty vs Time",
             fontsize=11.5, fontweight='bold', color='#1B3A6B')
ax.set_xlabel('Time (hours)'); ax.set_ylabel('Position 1-sigma (m)')
ax.set_yscale('log')
ax.legend(fontsize=9.5); ax.grid(alpha=0.3, which='both'); ax.set_xlim(0, 24)

ax = axes[1, 1]
ax.plot(df_['t_hours'], df_['kf_sigma'],  color=KF_COLOR,  lw=1.2, label='KF 1-sigma')
ax.plot(df_['t_hours'], df_['rts_sigma'], color=RTS_COLOR, lw=1.6, label='RTS 1-sigma')
ax.set_title("Dual-Frequency — Position Uncertainty vs Time",
             fontsize=11.5, fontweight='bold', color='#1B3A6B')
ax.set_xlabel('Time (hours)'); ax.set_ylabel('Position 1-sigma (m)')
ax.set_yscale('log')
ax.legend(fontsize=9.5); ax.grid(alpha=0.3, which='both'); ax.set_xlim(0, 24)

fig.suptitle('BILL00USA — KF vs RTS: Single-Frequency vs Dual-Frequency (24 hours)',
             fontsize=15, fontweight='bold', color='#1B3A6B', y=1.01)
plt.tight_layout()
plt.savefig('kf_rts_freq_comparison_crctd.png', dpi=200, bbox_inches='tight', facecolor='white')
print("Saved kf_rts_freq_comparison_crctd.png")
