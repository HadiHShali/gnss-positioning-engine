# plot_trajectory.py
# portfolio plots for gnss-positioning-engine spp output


# ========================================================== #
#                       Header + Data Loading                #
# ========================================================== #
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 

# Load trajectory csv
df = pd.read_csv('../data/trajectory_crctd.csv')
print(f'Loaded {len(df)} epochs')
print(df.head())

# Keep only converged fixes
df = df[df['converged']==1].reset_index(drop=True)
print(f'After Filtering Non-Converged Fixes: {len(df)} epochs')

#True Position for BILL00USA
# BILL00USA truth, converted from IGS XYZ (-2420420.501, -4737131.551, 3507827.665), WGS84
TRUE_LAT = 33.578243000    # deg
TRUE_LON = -117.064599000  # deg
TRUE_H   = 470.061         # m, ellipsoidal

# Convert Lat/Lon to Local East North Errors
R_Earth = 6378137.0 #m
Lat_rad = np.radians(TRUE_LAT)
df['east_err'] = np.radians(df['lon_deg'] - TRUE_LON) * R_Earth * np.cos(Lat_rad)
df['north_err'] = np.radians(df['lat_deg'] - TRUE_LAT) * R_Earth
df['height_err'] = df['height_m'] - TRUE_H
df['3d_err'] = np.sqrt(df['east_err']**2 + df['north_err']**2 + df['height_err']**2)    

#Time axis in hours
df['t_hours'] = df['t_gps'] / 3600.0

#Summary Statistics
print(f"Mean 3D Error: {df['3d_err'].mean():.2f} m")
print(f"RMS 3D Error: {np.sqrt((df['3d_err']**2).mean()):.2f} m")
print(f"Median: 3D Error: {df['3d_err'].median():.2f} m")


# ========================================================== #
#                       Plots                                #
# ========================================================== #
# Create a 2*2 subplot grid
fig, axes = plt.subplots(2, 2, figsize=(14, 10), facecolor='w', edgecolor='k')

#Plot 1: Horizontal Position Scatter Plot
ax = axes[0, 0]
ax.scatter(df['east_err'], df['north_err'],
           s=8, c='#2E75B6', alpha=0.4, edgecolors='none')
ax.scatter([0], [0], s=300, c='#1E7145', marker='*',
           edgecolors='black', linewidths=2, zorder=5, label='Truth')
ax.scatter([df['east_err'].mean()], [df['north_err'].mean()],
           s=200, c='#C55A11', marker='X', edgecolors='black',
           linewidths=2, zorder=5, label='Mean')
ax.axhline(0, color='#404040', lw=0.5)
ax.axvline(0, color='#404040', lw=0.5)
ax.set_xlabel('East error (m)')
ax.set_ylabel('North error (m)')
ax.set_title('Horizontal Scatter', fontsize=12, fontweight='bold', color='#1B3A6B')
ax.grid(True, alpha=0.3)
ax.set_aspect('equal')
ax.legend(loc='upper left', fontsize=9)

#Plot2 : 3D Error over Time
ax = axes[0, 1]
ax.plot(df['t_hours'], df['3d_err'], color='#2E75B6', lw=1.5, alpha=0.7)
mean_err = df['3d_err'].mean()
ax.axhline(mean_err, color='#C55A11', lw=2, linestyle='--',
           label=f'Mean = {mean_err:.1f} m')
ax.set_xlabel('Time (hours)')
ax.set_ylabel('3D error (m)')
ax.set_title('3D Position Error vs Time', fontsize=12, fontweight='bold', color='#1B3A6B')
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 24)
ax.legend(loc='upper right', fontsize=9)

# plot 3: Error Histogram
ax = axes[1, 0]
ax.hist(df['3d_err'], bins=40, color='#5B2C8D', alpha=0.7, edgecolor='black')
ax.axvline(df['3d_err'].mean(), color='#C55A11', lw=2, linestyle='--',
           label=f"Mean = {df['3d_err'].mean():.1f} m")
ax.axvline(df['3d_err'].median(), color='#1E7145', lw=2, linestyle='--',
           label=f"Median = {df['3d_err'].median():.1f} m")
ax.set_xlabel('3D error (m)')
ax.set_ylabel('Count')
ax.set_title(f"Error Distribution ({len(df)} epochs)",
             fontsize=12, fontweight='bold', color='#1B3A6B')
ax.grid(True, alpha=0.3)
ax.legend(loc='upper right', fontsize=9)

# plot 4: Satellite Count over Time
ax = axes[1, 1]
ax.plot(df['t_hours'], df['n_sats'], color='#1E7145', lw=1.5, alpha=0.7)
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Satellites used')
ax.set_title('Satellites in Fix vs Time',
             fontsize=12, fontweight='bold', color='#1B3A6B')
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 24)
ax.set_ylim(3, 13)

#save and show
# ── FINAL LAYOUT + SAVE ──
fig.suptitle('BILL00USA — SPP Processor Output — Jan 7, 2024',
             fontsize=14, fontweight='bold', color='#1B3A6B', y=1.00)
plt.tight_layout()
plt.savefig('../data/trajectory_analysis_crctd.png', dpi=200,
            bbox_inches='tight', facecolor='white')
plt.close()
print('Saved trajectory_analysis_crctd.png')
