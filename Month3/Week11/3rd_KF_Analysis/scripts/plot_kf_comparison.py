# plot_kf_comparison.py
# Compare raw SPP vs Kalman-filtered on real BILL00USA data
 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
 
# Load raw SPP trajectory (Week 9)
raw = pd.read_csv('../data/trajectory_crctd.csv')
raw = raw[raw['converged'] == 1].reset_index(drop=True)
 
# Load KF-filtered trajectory (Week 11)
kf = pd.read_csv('../data/trajectory_kf_crctd.csv')
 
print(f'Raw epochs: {len(raw)}')
print(f'KF epochs:  {len(kf)}')
 
# BILL00USA truth (from RINEX obs header APPROX POSITION XYZ)
TRUE_LAT = 33.578243000    # deg
TRUE_LON = -117.064599000  # deg
TRUE_H   = 470.061         # m, ellipsoidal
# Convert lat/lon differences to East/North meters
R_EARTH = 6378137.0
lat_rad = np.radians(TRUE_LAT)
 
# Raw errors
raw['east_err']  = np.radians(raw['lon_deg'] - TRUE_LON) * R_EARTH * np.cos(lat_rad)
raw['north_err'] = np.radians(raw['lat_deg'] - TRUE_LAT) * R_EARTH
raw['h_err']     = raw['height_m'] - TRUE_H
raw['3d_err']    = np.sqrt(raw['east_err']**2 + raw['north_err']**2 + raw['h_err']**2)
raw['t_hours']   = raw['t_gps'] / 3600.0
 
# KF errors
kf['east_err']   = np.radians(kf['kf_lon'] - TRUE_LON) * R_EARTH * np.cos(lat_rad)
kf['north_err']  = np.radians(kf['kf_lat'] - TRUE_LAT) * R_EARTH
kf['h_err']      = kf['kf_h'] - TRUE_H
kf['3d_err']     = np.sqrt(kf['east_err']**2 + kf['north_err']**2 + kf['h_err']**2)
kf['t_hours']    = kf['t_gps'] / 3600.0
 
# Summary statistics
print(f"Raw 3D error: mean={raw['3d_err'].mean():.1f}m, RMS={np.sqrt((raw['3d_err']**2).mean()):.1f}m")
print(f"KF  3D error: mean={kf['3d_err'].mean():.1f}m, RMS={np.sqrt((kf['3d_err']**2).mean()):.1f}m")

# 2x2 grid of comparison plots
fig, axes = plt.subplots(2, 2, figsize=(14, 10), facecolor='white')
# ── PLOT 1: Horizontal scatter — the story in one image ──
ax = axes[0, 0]
ax.scatter(raw['east_err'], raw['north_err'],
           s=6, c='#C55A11', alpha=0.3, edgecolors='none',
           label=f'Raw SPP ({len(raw)} fixes)')
ax.plot(kf['east_err'], kf['north_err'],
        color='#2E75B6', lw=1.5, alpha=0.9, label='KF trajectory')
ax.scatter([0], [0], s=350, c='#1E7145', marker='*',
           edgecolors='black', linewidths=2, zorder=10, label='Truth')
ax.axhline(0, color='#404040', lw=0.5)
ax.axvline(0, color='#404040', lw=0.5)
ax.set_xlabel('East error (m)')
ax.set_ylabel('North error (m)')
ax.set_title('Horizontal Scatter — Raw vs Filtered',
             fontsize=12, fontweight='bold', color='#1B3A6B')
ax.grid(True, alpha=0.3); ax.legend(loc='upper right', fontsize=9)
ax.set_aspect('equal')
# ── PLOT 2: 3D error time series ──
ax = axes[0, 1]
ax.plot(raw['t_hours'], raw['3d_err'],
        color='#C55A11', lw=0.4, alpha=0.5, label='Raw SPP error')
ax.plot(kf['t_hours'], kf['3d_err'],
        color='#2E75B6', lw=1.5, label='KF error')
ax.set_xlabel('Time (hours)')
ax.set_ylabel('3D error (m)')
ax.set_title('3D Position Error vs Time',
             fontsize=12, fontweight='bold', color='#1B3A6B')
ax.grid(True, alpha=0.3); ax.legend(loc='upper right', fontsize=9)
ax.set_xlim(0, 24)
# ── PLOT 3: Error distribution comparison ──
ax = axes[1, 0]
ax.hist(raw['3d_err'], bins=40, color='#C55A11', alpha=0.6,
        edgecolor='black',
        label=f"Raw (mean={raw['3d_err'].mean():.1f}m)")
ax.hist(kf['3d_err'], bins=40, color='#2E75B6', alpha=0.6,
        edgecolor='black',
        label=f"KF (mean={kf['3d_err'].mean():.1f}m)")
ax.set_xlabel('3D error (m)')
ax.set_ylabel('Count')
ax.set_title('Error Distribution Comparison',
             fontsize=12, fontweight='bold', color='#1B3A6B')
ax.grid(True, alpha=0.3); ax.legend(loc='upper right', fontsize=9)
# ── PLOT 4: KF uncertainty over time ──
ax = axes[1, 1]
ax.plot(kf['t_hours'], kf['sqrt_P00'], color='#5B2C8D', lw=2,
        label='KF position 1-sigma (x)')
ax.axhline(raw['sigma0_m'].mean(), color='#C55A11', lw=1.5, linestyle='--',
           alpha=0.7,
           label=f"Raw sigma-0 mean ({raw['sigma0_m'].mean():.1f}m)")
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Position uncertainty (m)')
ax.set_title('Filter Confidence Grows Over Time',
             fontsize=12, fontweight='bold', color='#1B3A6B')
ax.set_yscale('log')
ax.grid(True, alpha=0.3, which='both')
ax.legend(loc='upper right', fontsize=9)
ax.set_xlim(0, 24)
fig.suptitle('BILL00USA — Raw SPP vs Kalman-Filtered (24 hours)',
             fontsize=14, fontweight='bold', color='#1B3A6B', y=1.00)
plt.tight_layout()
plt.savefig('../data/kf_comparison_crctd.png', dpi=200,
            bbox_inches='tight', facecolor='white')
plt.savefig('../data/kf_comparison.pdf', bbox_inches='tight')
plt.close()
print('Saved kf_comparison_crctd.png and .pdf')

