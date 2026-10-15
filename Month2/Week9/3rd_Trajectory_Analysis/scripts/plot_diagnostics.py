# ========================================================================#
#                                   Header + Load                         #
# ========================================================================#
# plot_diagnostics.py
# Engineering diagnostic plots for gnss-positioning-engine
 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
 
# ── LOAD DATA ──
df = pd.read_csv('../data/trajectory_crctd.csv')
df = df[df['converged'] == 1].reset_index(drop=True)
 
# Truth for 3D error computation
# BILL00USA truth, converted from IGS XYZ (-2420420.501, -4737131.551, 3507827.665), WGS84
TRUE_LAT = 33.578243000    # deg
TRUE_LON = -117.064599000  # deg
TRUE_H   = 470.061         # m, ellipsoidal
R_EARTH  = 6378137.0
 
lat_rad = np.radians(TRUE_LAT)
df['east_err']  = np.radians(df['lon_deg'] - TRUE_LON) * R_EARTH * np.cos(lat_rad)
df['north_err'] = np.radians(df['lat_deg'] - TRUE_LAT) * R_EARTH
df['h_err']     = df['height_m'] - TRUE_H
df['3d_err']    = np.sqrt(df['east_err']**2 + df['north_err']**2 + df['h_err']**2)
df['t_hours']   = df['t_gps'] / 3600.0
 
print(f'Loaded {len(df)} converged epochs')

# ========================================================================#
#                         Setup 2x2 Figure                                #
# ========================================================================#
fig, axes = plt.subplots(2, 2, figsize=(14, 10), facecolor='white')

# ── PLOT 1: PDOP histogram ──
ax = axes[0, 0]
ax.hist(df['pdop'], bins=40, color='#5B2C8D', alpha=0.7, edgecolor='black')
ax.axvline(df['pdop'].mean(), color='#C55A11', lw=2, linestyle='--',
           label=f"Mean PDOP = {df['pdop'].mean():.2f}")
ax.axvline(3.0, color='#1E7145', lw=1.5, linestyle=':',
           label='Good/Marginal threshold')
ax.set_xlabel('PDOP')
ax.set_ylabel('Count')
ax.set_title('PDOP Distribution',
             fontsize=12, fontweight='bold', color='#1B3A6B')
ax.grid(True, alpha=0.3)
ax.legend(loc='upper right', fontsize=9)

# ── PLOT 2: PDOP time series ──
ax = axes[0, 1]
ax.plot(df['t_hours'], df['pdop'], color='#5B2C8D', lw=0.5, alpha=0.7)
ax.axhline(3.0, color='#1E7145', lw=1.5, linestyle=':', label='Good (< 3)')
ax.axhline(6.0, color='#C55A11', lw=1.5, linestyle=':', label='Marginal (< 6)')
ax.set_xlabel('Time (hours)')
ax.set_ylabel('PDOP')
ax.set_title('PDOP vs Time',
             fontsize=12, fontweight='bold', color='#1B3A6B')
ax.grid(True, alpha=0.3)
ax.legend(loc='upper right', fontsize=9)
ax.set_xlim(0, 24)

# ── PLOT 3: Correlation check ──
ax = axes[1, 0]
ax.scatter(df['pdop'], df['sigma0_m'], s=6, c='#2E75B6',
           alpha=0.4, edgecolors='none')
# Fit a line to show the correlation
z = np.polyfit(df['pdop'], df['sigma0_m'], 1)
xfit = np.linspace(df['pdop'].min(), df['pdop'].max(), 100)
ax.plot(xfit, z[0]*xfit + z[1], color='#C55A11', lw=2,
        label=f'y = {z[0]:.1f}x + {z[1]:.1f}')
ax.set_xlabel('PDOP')
ax.set_ylabel('sigma-0 (m)')
ax.set_title('Solution Quality Correlation:  sigma-0 vs PDOP',
             fontsize=12, fontweight='bold', color='#1B3A6B')
ax.grid(True, alpha=0.3)
ax.legend(loc='upper left', fontsize=9)

# ── PLOT 4: 3D error boxplot by n_sats ──
ax = axes[1, 1]
# Group by number of satellites
unique_counts = sorted(df['n_sats'].unique())
box_data = [df.loc[df['n_sats']==k, '3d_err'].values
            for k in unique_counts if (df['n_sats']==k).sum() > 5]
box_labels = [k for k in unique_counts if (df['n_sats']==k).sum() > 5]
bp = ax.boxplot(box_data, tick_labels=box_labels, patch_artist=True)
for patch in bp['boxes']:
    patch.set_facecolor('#E2EFDA')
    patch.set_edgecolor('#1E7145')
ax.set_xlabel('Number of satellites used')
ax.set_ylabel('3D error (m)')
ax.set_title('3D Error by Satellite Count',
             fontsize=12, fontweight='bold', color='#1B3A6B')
ax.grid(True, alpha=0.3)

# saVE
fig.suptitle('BILL00USA — Engineering Diagnostics',
             fontsize=14, fontweight='bold', color='#1B3A6B', y=1.00)
plt.tight_layout()
plt.savefig('../data/diagnostics_crctd.png', dpi=200,
            bbox_inches='tight', facecolor='white')
plt.close()
print('Saved diagnostics_crctd.png')
