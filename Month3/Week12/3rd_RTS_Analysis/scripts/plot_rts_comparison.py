import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
 
raw = pd.read_csv('../data/trajectory_crctd.csv')
raw = raw[raw['converged'] == 1].reset_index(drop=True)
rts = pd.read_csv('../data/trajectory_rts_crctd.csv')
 

TRUE_LAT = 33.578243000    # deg
TRUE_LON = -117.064599000  # deg
TRUE_H   = 470.061         # m, ellipsoidal

R_EARTH  = 6378137.0
lat_rad  = np.radians(TRUE_LAT)
 
def add_errors(df, lat_col, lon_col, h_col, prefix):
    df[prefix+'_east']  = np.radians(df[lon_col]-TRUE_LON)*R_EARTH*np.cos(lat_rad)
    df[prefix+'_north'] = np.radians(df[lat_col]-TRUE_LAT)*R_EARTH
    df[prefix+'_h']     = df[h_col]-TRUE_H
    df[prefix+'_3d']    = np.sqrt(df[prefix+'_east']**2 +
                                   df[prefix+'_north']**2 + df[prefix+'_h']**2)
    return df
 
raw = add_errors(raw, 'lat_deg', 'lon_deg', 'height_m', 'raw')
rts = add_errors(rts, 'kf_lat',  'kf_lon',  'kf_h',     'kf')
rts = add_errors(rts, 'rts_lat', 'rts_lon', 'rts_h',    'rts')
rts['t_hours'] = rts['t_gps'] / 3600.0
raw['t_hours'] = raw['t_gps'] / 3600.0
 
print(f"Raw   3D error: mean={raw['raw_3d'].mean():.2f}m  RMS={np.sqrt((raw['raw_3d']**2).mean()):.2f}m")
print(f"KF    3D error: mean={rts['kf_3d'].mean():.2f}m  RMS={np.sqrt((rts['kf_3d']**2).mean()):.2f}m")
print(f"RTS   3D error: mean={rts['rts_3d'].mean():.2f}m  RMS={np.sqrt((rts['rts_3d']**2).mean()):.2f}m")
 
kf_rms  = np.sqrt((rts['kf_3d']**2).mean())
rts_rms = np.sqrt((rts['rts_3d']**2).mean())
improvement = (kf_rms - rts_rms) / kf_rms * 100
print(f"RTS improvement over KF: {improvement:.1f}%  (literature analog: ~10%)")



fig, axes = plt.subplots(2, 2, figsize=(14, 10), facecolor='white')
 
# PLOT 1: Horizontal scatter — all 3
ax = axes[0, 0]
ax.scatter(raw['raw_east'], raw['raw_north'], s=5, c='#C55A11', alpha=0.2,
           edgecolors='none', label='Raw SPP')
ax.plot(rts['kf_east'], rts['kf_north'], color='#2E75B6', lw=1.2, alpha=0.8,
        label='KF (Week 11)')
ax.plot(rts['rts_east'], rts['rts_north'], color='#1E7145', lw=1.8,
        label='RTS smoothed (new)')
ax.scatter([0],[0], s=300, c='gold', marker='*', edgecolors='black',
           linewidths=1.5, zorder=10, label='Truth')
ax.axhline(0, color='#404040', lw=0.5); ax.axvline(0, color='#404040', lw=0.5)
ax.set_xlabel('East error (m)'); ax.set_ylabel('North error (m)')
ax.set_title('Horizontal Scatter — All 3 Stages', fontweight='bold', color='#1B3A6B')
ax.legend(fontsize=8); ax.grid(alpha=0.3); ax.set_aspect('equal')
 
# PLOT 2: 3D error vs time — KF vs RTS only (raw too noisy to see)
ax = axes[0, 1]
ax.plot(rts['t_hours'], rts['kf_3d'], color='#2E75B6', lw=1, alpha=0.6, label='KF')
ax.plot(rts['t_hours'], rts['rts_3d'], color='#1E7145', lw=1.5, label='RTS')
ax.set_xlabel('Time (hours)'); ax.set_ylabel('3D error (m)')
ax.set_title('KF vs RTS Error Over Time', fontweight='bold', color='#1B3A6B')
ax.legend(fontsize=9); ax.grid(alpha=0.3); ax.set_xlim(0,24)
 
# PLOT 3: Histograms — all 3
ax = axes[1, 0]
ax.hist(raw['raw_3d'], bins=40, color='#C55A11', alpha=0.4, edgecolor='black',
        label=f"Raw (RMS={np.sqrt((raw['raw_3d']**2).mean()):.1f}m)")
ax.hist(rts['kf_3d'], bins=40, color='#2E75B6', alpha=0.5, edgecolor='black',
        label=f"KF (RMS={kf_rms:.1f}m)")
ax.hist(rts['rts_3d'], bins=40, color='#1E7145', alpha=0.5, edgecolor='black',
        label=f"RTS (RMS={rts_rms:.1f}m)")
ax.set_xlabel('3D error (m)'); ax.set_ylabel('Count')
ax.set_title('Error Distribution — All 3 Stages', fontweight='bold', color='#1B3A6B')
ax.legend(fontsize=8); ax.grid(alpha=0.3)
 
# PLOT 4: Uncertainty comparison
ax = axes[1, 1]
ax.plot(rts['t_hours'], rts['kf_sigma'], color='#2E75B6', lw=1.5, label='KF 1-sigma')
ax.plot(rts['t_hours'], rts['rts_sigma'], color='#1E7145', lw=1.5, label='RTS 1-sigma')
ax.set_xlabel('Time (hours)'); ax.set_ylabel('Position uncertainty (m)')
ax.set_title('Confidence: KF vs RTS (RTS always <= KF)', fontweight='bold', color='#1B3A6B')
ax.set_yscale('log'); ax.legend(fontsize=9); ax.grid(alpha=0.3, which='both')
ax.set_xlim(0,24)
 
fig.suptitle('BILL00USA — Raw vs KF vs RTS-Smoothed (24 hours)',
             fontsize=14, fontweight='bold', color='#1B3A6B', y=1.00)
plt.tight_layout()
plt.savefig('../data/rts_comparison_crctd.png', dpi=200, bbox_inches='tight', facecolor='white')
print('Saved rts_comparison_crctd.png')
