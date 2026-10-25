"""
plot_final_sp3clk_comparison_crctd.py

CORRECTED copy of plot_final_sp3clk_comparison.py (Month 3 / Week 14 review, 2026-10-24).
The original stays untouched. Changes are marked [CRCTD]:
  1. Input: the corrected Day 4 output, trajectory_dualfreq_SP3CLK_RTS_crctd.csv.
  2. TRUTH: the old (33.578381, -117.064265, 472.0) came from an early, buggy SPP run and is
     34.6 m away from BILL's IGS position (found in the Month 2 Week 9 review). Truth is now
     computed from the IGS ECEF coordinates, the same numbers the C++ code uses.
  3. Errors: exact ECEF -> local East/North/Up rotation at the truth point (was a
     spherical-Earth approximation, ~0.4% scale error; small, but no reason to keep it).
  4. Benchmark: 13.84 m was the reception-time-bug pipeline. The corrected references are
     Week 13 (dual-freq, broadcast, KF+RTS) = 0.797 m and Week 12 (single-freq) = 2.502 m.
  5. Scatter: the RTS result is ONE point (static station, Q = 0), so it is drawn as a
     marker; as a line it was hidden under the KF track.
  6. Labels: kf_sigma / rts_sigma are sqrt(P_xx) = the ECEF X sigma, so the axis says so.
  7. Prints the RTS offset in E/N/U, to cross-check the C++ summary line.

Produces a 4-panel figure:
    Top-left:     Horizontal scatter — KF trajectory vs RTS-smoothed vs Truth
    Top-right:    3D error vs time — KF vs RTS
    Bottom-left:  Error distribution histogram — KF vs RTS
    Bottom-right: Position uncertainty vs time — KF vs RTS (log scale)

Run from a subfolder of Month3/Week14/4th_SP3_KF_RTS/ (CSV path is ../data/...), as before.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ─────────────────────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────────────────────
CSV_PATH = "../data/trajectory_dualfreq_SP3CLK_RTS_crctd.csv"   # [CRCTD] 1

# [CRCTD] 2. BILL00USA IGS position (ECEF, m) — same values as TRUE_X/Y/Z in the C++ code
TRUE_X, TRUE_Y, TRUE_Z = -2420420.501, -4737131.551, 3507827.665

KF_COLOR  = "#2E75B6"
RTS_COLOR = "#1E7145"
TRUTH_COLOR = "gold"

# [CRCTD] 4. Corrected references (static KF + RTS, transmit-time-corrected SPP)
BENCH_DUAL_BRDC = 0.797    # Week 13: dual-frequency, broadcast orbits/clocks
BENCH_SINGLE    = 2.502    # Week 12: single-frequency + Klobuchar


# ─────────────────────────────────────────────────────────────
# [CRCTD] 3. Exact WGS84 geometry
# ─────────────────────────────────────────────────────────────
A_WGS = 6378137.0
F_WGS = 1.0 / 298.257223563
E2    = F_WGS * (2.0 - F_WGS)


def lla_to_ecef(lat_deg, lon_deg, h):
    """Geodetic (deg, deg, m) -> ECEF (m). Works on numpy arrays / pandas Series."""
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    N = A_WGS / np.sqrt(1.0 - E2 * np.sin(lat) ** 2)
    return ((N + h) * np.cos(lat) * np.cos(lon),
            (N + h) * np.cos(lat) * np.sin(lon),
            (N * (1.0 - E2) + h) * np.sin(lat))


def ecef_to_lla(x, y, z):
    """ECEF (m) -> geodetic (deg, deg, m), iterative."""
    p = np.hypot(x, y)
    lat = np.arctan2(z, p * (1.0 - E2))
    for _ in range(10):
        N = A_WGS / np.sqrt(1.0 - E2 * np.sin(lat) ** 2)
        h = p / np.cos(lat) - N
        lat = np.arctan2(z, p * (1.0 - E2 * N / (N + h)))
    return np.degrees(lat), np.degrees(np.arctan2(y, x)), h


TRUE_LAT, TRUE_LON, TRUE_H = ecef_to_lla(TRUE_X, TRUE_Y, TRUE_Z)   # 33.578243, -117.064599, 470.061
_phi, _lam = np.radians(TRUE_LAT), np.radians(TRUE_LON)


def err_components(lat, lon, h):
    """East / North / Up error (m) vs BILL00USA truth: ECEF difference rotated into local ENU."""
    x, y, z = lla_to_ecef(lat, lon, h)
    dx, dy, dz = x - TRUE_X, y - TRUE_Y, z - TRUE_Z
    east  = -np.sin(_lam) * dx + np.cos(_lam) * dy
    north = -np.sin(_phi) * np.cos(_lam) * dx - np.sin(_phi) * np.sin(_lam) * dy + np.cos(_phi) * dz
    up    =  np.cos(_phi) * np.cos(_lam) * dx + np.cos(_phi) * np.sin(_lam) * dy + np.sin(_phi) * dz
    return east, north, up


# ─────────────────────────────────────────────────────────────
# LOAD + VERIFY COLUMNS + PROCESS
# ─────────────────────────────────────────────────────────────
df = pd.read_csv(CSV_PATH)

# Safety check: confirm the actual column names before assuming they match
print("Columns found in CSV:", df.columns.tolist())
print(f"Truth (from IGS ECEF): lat {TRUE_LAT:.9f}  lon {TRUE_LON:.9f}  h {TRUE_H:.4f} m")
print()

required_cols = ['kf_lat', 'kf_lon', 'kf_h', 'rts_lat', 'rts_lon', 'rts_h']
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise SystemExit(
        f"ERROR: expected columns {missing} not found in CSV.\n"
        f"Update the column names below to match what's actually in your file."
    )

df['t_hours'] = df['t_gps'] / 3600.0   # 0-24 h only because BILL's file is a Sunday

kf_east, kf_north, kf_up    = err_components(df['kf_lat'],  df['kf_lon'],  df['kf_h'])
rts_east, rts_north, rts_up = err_components(df['rts_lat'], df['rts_lon'], df['rts_h'])

df['kf_3d']  = np.sqrt(kf_east**2  + kf_north**2  + kf_up**2)
df['rts_3d'] = np.sqrt(rts_east**2 + rts_north**2 + rts_up**2)

kf_err  = df['kf_3d']
rts_err = df['rts_3d']
kf_mean, kf_rms   = kf_err.mean(),  np.sqrt((kf_err**2).mean())
rts_mean, rts_rms = rts_err.mean(), np.sqrt((rts_err**2).mean())

print(f"SP3/CLK + KF:        mean={kf_mean:.3f}m   RMS={kf_rms:.3f}m")
print(f"SP3/CLK + KF + RTS:  mean={rts_mean:.3f}m   RMS={rts_rms:.3f}m")
# [CRCTD] 7. RTS is flat for a static station: report its E/N/U offset (compare with the C++ line)
print(f"RTS offset (E, N, U): ({rts_east.mean():+.3f}, {rts_north.mean():+.3f}, {rts_up.mean():+.3f}) m"
      f"   horizontal {np.hypot(rts_east.mean(), rts_north.mean()):.3f} m")
print(f"[Reference] Week 13 dual-freq broadcast + KF + RTS: {BENCH_DUAL_BRDC:.3f} m")
print(f"[Reference] Week 12 single-freq + KF + RTS:         {BENCH_SINGLE:.3f} m")
print()

# ─────────────────────────────────────────────────────────────
# FIGURE — 4 panels
# ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(15, 11), facecolor='white')

# ---------- TOP-LEFT: Horizontal scatter ----------
ax = axes[0, 0]
ax.plot(kf_east, kf_north, color=KF_COLOR, lw=0.8, alpha=0.6,
        label=f'KF (mean 3D={kf_mean:.2f}m)')
# [CRCTD] 5. RTS = one point for a static station -> marker, drawn on top of the KF track
ax.scatter([rts_east.mean()], [rts_north.mean()], s=120, c=RTS_COLOR, marker='o',
           edgecolors='black', linewidths=1.2, zorder=9,
           label=f'RTS (3D={rts_mean:.2f}m)')
ax.scatter([0], [0], s=400, c=TRUTH_COLOR, marker='*',
           edgecolors='black', linewidths=1.8, zorder=10, label='Truth (IGS)')
ax.axhline(0, color='#404040', lw=0.5)
ax.axvline(0, color='#404040', lw=0.5)
ax.set_xlabel('East error (m)', fontsize=11)
ax.set_ylabel('North error (m)', fontsize=11)
ax.set_title('Horizontal Scatter — SP3/CLK + KF vs RTS',
              fontsize=12.5, fontweight='bold', color='#1B3A6B')
ax.legend(fontsize=10, loc='best')
ax.grid(alpha=0.3)
ax.set_aspect('equal')

# ---------- TOP-RIGHT: 3D error vs time ----------
ax = axes[0, 1]
ax.plot(df['t_hours'], df['kf_3d'], color=KF_COLOR, lw=0.8, alpha=0.6, label='KF')
ax.plot(df['t_hours'], df['rts_3d'], color=RTS_COLOR, lw=1.8, label='RTS')
ax.axhline(BENCH_DUAL_BRDC, color='#C55A11', lw=1.5, linestyle='--', alpha=0.8,      # [CRCTD] 4
           label=f'Wk13 dual-freq broadcast RTS ({BENCH_DUAL_BRDC:.2f}m)')
ax.set_xlabel('Time (hours)', fontsize=11)
ax.set_ylabel('3D error (m)', fontsize=11)
ax.set_title('3D Error vs Time', fontsize=12.5, fontweight='bold', color='#1B3A6B')
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
ax.set_xlim(0, 24)
ax.set_ylim(bottom=0)

# ---------- BOTTOM-LEFT: Error histogram ----------
ax = axes[1, 0]
bins = np.linspace(0, max(df['kf_3d'].max(), df['rts_3d'].max()) * 1.05, 45)
ax.hist(df['kf_3d'], bins=bins, color=KF_COLOR, alpha=0.5, edgecolor='black',
        label=f'KF (RMS={kf_rms:.2f}m)')
ax.hist(df['rts_3d'], bins=bins, color=RTS_COLOR, alpha=0.6, edgecolor='black',
        label=f'RTS (RMS={rts_rms:.2f}m)')
ax.axvline(BENCH_DUAL_BRDC, color='#C55A11', lw=1.5, linestyle='--', alpha=0.8,       # [CRCTD] 4
           label='Wk13 broadcast RTS')
ax.set_xlabel('3D error (m)', fontsize=11)
ax.set_ylabel('Count', fontsize=11)
ax.set_title('Error Distribution', fontsize=12.5, fontweight='bold', color='#1B3A6B')
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

# ---------- BOTTOM-RIGHT: Uncertainty vs time ----------
ax = axes[1, 1]
if 'kf_sigma' in df.columns and 'rts_sigma' in df.columns:
    ax.plot(df['t_hours'], df['kf_sigma'], color=KF_COLOR, lw=1.3, label='KF 1-sigma (X)')
    ax.plot(df['t_hours'], df['rts_sigma'], color=RTS_COLOR, lw=1.8, label='RTS 1-sigma (X)')
    ax.set_yscale('log')
    ax.legend(fontsize=10)
else:
    ax.text(0.5, 0.5, 'sigma columns not found in CSV\n(check kf_sigma / rts_sigma column names)',
            ha='center', va='center', transform=ax.transAxes, fontsize=10, color='#999')
ax.set_xlabel('Time (hours)', fontsize=11)
ax.set_ylabel('ECEF X 1-sigma (m)', fontsize=11)       # [CRCTD] 6. sqrt(P_xx), not a 3D sigma
ax.set_title('Position Uncertainty vs Time', fontsize=12.5, fontweight='bold', color='#1B3A6B')
ax.grid(alpha=0.3, which='both')
ax.set_xlim(0, 24)

fig.suptitle('Month 3 Finale — BILL00USA: Dual-Freq + SP3/CLK + KF/RTS (24 hours, corrected)',
              fontsize=15.5, fontweight='bold', color='#1B3A6B', y=1.01)
plt.tight_layout()
plt.savefig('final_sp3clk_comparison_crctd.png', dpi=200, bbox_inches='tight', facecolor='white')
print("\nSaved final_sp3clk_comparison_crctd.png")
