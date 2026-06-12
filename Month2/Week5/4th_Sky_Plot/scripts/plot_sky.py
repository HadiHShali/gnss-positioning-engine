# ================================================================
# plot_sky.py - Generate a sky plot from visible_sats.csv
# ================================================================

from matplotlib.lines import Line2D
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np

#load data
print('loading visible satellites...')
df = pd.read_csv('./../data/visible_sats.csv')
print(f'{len(df)} satellites loaded.')
print(df)

# Generate the Figure
fig = plt.figure(figsize=(10, 10), facecolor='white')
ax = fig.add_subplot(111, projection='polar')

# Convert El/Az to polar coordinates
# Azimuth is degrees  -> theta in radians
theta = np.deg2rad(df['azimuth_deg'])
# Elevation: high = center, low = edge -> r = 90 - elevation
r = 90 - df['elevation_deg']

# color by satellite elevation (high = green, low=orange)
colors = []
for el in df['elevation_deg']:
    if el >= 45:
        colors.append('#1E7145') # Green: high, clear
    elif el >= 30:
        colors.append('#FFD400') # Yellow: medium, good
    else:
        colors.append('#FF6B35') #Orange: low, noisy

# Plot the satellites
scatter = ax.scatter(theta, r, c=colors, s=400, edgecolors='black', linewidth=2, zorder=5)

# Add PRN label at top of each point
for i, row in df.iterrows():
    ax.text(theta[i], r[i], str(row['prn']), fontsize=9, fontweight='bold', color ='white', ha='center', va='center', zorder=6)

# Configure Axes - compass style
ax.set_theta_zero_location('N')  # 0 degrees at the top (North)
ax.set_theta_direction(-1)       # Clockwise
ax.set_rlim(0, 90)              # Elevation from 0 to 90 degrees

# Elevation rings labels
ax.set_yticks([0, 30, 60, 90])
ax.set_yticklabels(['90° (Zenith)', '60°', '30°', '0° (Horizon)'], fontsize=9)

# compass direction  lables
ax.set_xticks(np.deg2rad([0, 45, 90, 135, 180, 225, 270, 315]))
ax.set_xticklabels(['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'], fontsize=14, fontweight='bold')

# style the grid
ax.grid(color='gray', linestyle='--', linewidth=0.5, zorder=0)
ax.set_facecolor('#F8F8F8') # light gray background

# Receiver marker at the center 
ax.plot(0, 0, '^', ms=15, color='#2E75B6', zorder=10,
        markeredgecolor='white', markeredgewidth=2
)

# Title and Metadata
n_sats = len(df)
ax.set_title(f'GPS Sky View — Memphis, U of M Campus\n' +
             f'Date: 2024-01-07  Time: 12:00 UTC  |  '+
             f'{n_sats} visible satellites',
             fontsize=13, fontweight='bold', color='#1B3A6B', pad=25)

# Legend
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#1E7145',
           markeredgecolor='black', markersize=12, label='High (>= 45 deg)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#FFB400',
           markeredgecolor='black', markersize=12, label='Medium (25-45 deg)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#FF6B35',
           markeredgecolor='black', markersize=12, label='Low (< 25 deg)'),
]
ax.legend(handles=legend_elements, loc='lower right',
          bbox_to_anchor=(1.25, 0), fontsize=10, framealpha=0.95,
          title='Signal quality', title_fontsize=10)
 
# ── SAVE ─────────────────────────────────────────────────────────────────
plt.tight_layout()
plt.savefig('./../data/sky_plot.png', dpi=150, bbox_inches='tight',
            facecolor='white')
plt.savefig('./../data/sky_plot.pdf', bbox_inches='tight')
 
print('\nSaved:')
print('  ./../data/sky_plot.png  (raster, for README)')
print('  ./../data/sky_plot.pdf  (vector, for portfolio)')
 
plt.show()
