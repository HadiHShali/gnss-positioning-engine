# ================================================================
# plot_sky_crctd.py - Generate a sky plot from visible_sats_crctd.csv
# ================================================================
#
# _crctd = corrected copy of plot_sky.py (Week 5 Day 4), reviewed during the
# Week 16 correction pass. The original is kept unchanged.
# Corrections (search for "CORRECTED"):
#   1. Elevation classes: the plot used 45/30 deg while the legend said 45/25 deg,
#      and the legend's yellow (#FFB400) differed from the plotted one (#FFD400).
#      Both now come from ONE table (ELEV_CLASSES), so they cannot disagree.
#      Example of the old mismatch: G21 at 27.4 deg was drawn orange ("Low")
#      while the legend put 25-45 deg in "Medium".
#   2. PRN labels: df.iterrows() converts each row to ONE dtype (float, because
#      the other columns are floats), so labels printed as "24.0". Now "G24".
#   3. Title: 12:00:00 GPS time = 11:59:42 UTC (the old title said 12:00 UTC).
#   4. Health: unhealthy satellites (sv_health != 0) are drawn gray, not
#      colored as usable.
#   5. File names: reads visible_sats_crctd.csv, writes sky_plot_crctd.png/.pdf,
#      so the original outputs are not overwritten.
#   Run from the scripts/ folder, like the original.

from matplotlib.lines import Line2D
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ── CORRECTED: one table defines the classes for BOTH the dots and the legend ──
# (lower elevation limit in deg, color, label), highest class first
ELEV_CLASSES = [
    (45.0, '#1E7145', 'High (>= 45 deg)'),     # green: short path, little multipath
    (30.0, '#FFD400', 'Medium (30-45 deg)'),   # yellow
    (0.0,  '#FF6B35', 'Low (< 30 deg)'),       # orange: long atmospheric path, multipath
]
UNHEALTHY_COLOR = '#9E9E9E'                    # gray: broadcast health flag says "do not use"

def class_color(elev_deg):
    """Return the color of the first class whose lower limit the elevation reaches."""
    for lower, color, _ in ELEV_CLASSES:
        if elev_deg >= lower:
            return color
    return ELEV_CLASSES[-1][1]

# load data
print('loading visible satellites...')
df = pd.read_csv('./../data/visible_sats_crctd.csv')   # CORRECTED: _crctd input
print(f'{len(df)} satellites loaded.')
print(df)

# Generate the Figure
fig = plt.figure(figsize=(10, 10), facecolor='white')
ax = fig.add_subplot(111, projection='polar')

# Convert El/Az to polar coordinates
theta = np.deg2rad(df['azimuth_deg'])   # azimuth [deg] -> angle [rad]
r = 90 - df['elevation_deg']            # zenith at the center, horizon at the edge

# CORRECTED: colors from the shared class table; unhealthy satellites in gray
colors = [UNHEALTHY_COLOR if h != 0 else class_color(el)
          for el, h in zip(df['elevation_deg'], df['sv_health'])]

# Plot the satellites
ax.scatter(theta, r, c=colors, s=400, edgecolors='black', linewidth=2, zorder=5)

# CORRECTED: integer PRN labels. Reading the PRN column directly keeps its
# integer dtype; iterrows() would have turned it into a float (24 -> "24.0").
for th, rr, prn in zip(theta, r, df['prn'].astype(int)):
    ax.text(th, rr, f'G{prn:02d}', fontsize=8, fontweight='bold', color='white',
            ha='center', va='center', zorder=6)

# Configure Axes - compass style
ax.set_theta_zero_location('N')  # 0 degrees at the top (North)
ax.set_theta_direction(-1)       # Clockwise
ax.set_rlim(0, 90)               # Elevation from 0 to 90 degrees

# Elevation rings labels
ax.set_yticks([0, 30, 60, 90])
ax.set_yticklabels(['90° (Zenith)', '60°', '30°', '0° (Horizon)'], fontsize=9)

# compass direction labels
ax.set_xticks(np.deg2rad([0, 45, 90, 135, 180, 225, 270, 315]))
ax.set_xticklabels(['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'], fontsize=14, fontweight='bold')

# style the grid
ax.grid(color='gray', linestyle='--', linewidth=0.5, zorder=0)
ax.set_facecolor('#F8F8F8')

# Receiver marker at the center
ax.plot(0, 0, '^', ms=15, color='#2E75B6', zorder=10,
        markeredgecolor='white', markeredgewidth=2)

# Title and Metadata
n_sats = len(df)
ax.set_title('GPS Sky View — Memphis, U of M Campus\n' +
             'Date: 2024-01-07  Time: 12:00:00 GPS (11:59:42 UTC)  |  ' +   # CORRECTED
             f'{n_sats} visible satellites',
             fontsize=13, fontweight='bold', color='#1B3A6B', pad=25)

# CORRECTED: legend built from the same class table as the dots
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=color,
           markeredgecolor='black', markersize=12, label=label)
    for _, color, label in ELEV_CLASSES
]
if (df['sv_health'] != 0).any():               # only show the gray entry if it is used
    legend_elements.append(
        Line2D([0], [0], marker='o', color='w', markerfacecolor=UNHEALTHY_COLOR,
               markeredgecolor='black', markersize=12, label='Unhealthy (do not use)'))
ax.legend(handles=legend_elements, loc='lower right',
          bbox_to_anchor=(1.25, 0), fontsize=10, framealpha=0.95,
          title='Elevation class', title_fontsize=10)

# ── SAVE ─────────────────────────────────────────────────────────────────
plt.tight_layout()
plt.savefig('./../data/sky_plot_crctd.png', dpi=150, bbox_inches='tight',
            facecolor='white')
plt.savefig('./../data/sky_plot_crctd.pdf', bbox_inches='tight')

print('\nSaved:')
print('  ./../data/sky_plot_crctd.png  (raster, for README)')
print('  ./../data/sky_plot_crctd.pdf  (vector, for portfolio)')

plt.show()
