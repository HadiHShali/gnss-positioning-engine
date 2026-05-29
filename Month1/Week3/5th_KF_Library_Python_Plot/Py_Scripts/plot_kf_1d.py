#plot_kf_1d.py - Reads KF CSV output in data directory and plots the results
# Run from the project root with: python scripts/plot_kf_1d.py

#Section1: Imports
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#Section2: Read the CSV data
#pd.read_csv() loads the file into a DataFrame (like a spreadsheet)
df = pd.read_csv('../data/kf_1d_output.csv')

#print to verify the data loaded correctly
print('Data Loaded:')
print(df)

#SECTION3: Prepare the data for plotting
epochs = df['epoch']
gps_alt = df['gps_alt']
kf_alt = df['kf_alt']
uncertainty = df['uncertainty']

#standard deviation = sqrt(variance) - for the uncertainty band
std = np.sqrt(uncertainty)

#true altitude
true_alt = 1191.0  # meters

#SECTION4: Create the plot
#create a figure with specifc size (width, height) in inches
plt.figure(figsize=(10, 6))

# plot 1: GPS true altitude as horizontal dashed line
plt.axhline(true_alt, color = 'navy', linestyle = '--', linewidth = 2,
             label = f'True altitude: {true_alt} m')

# plot2: GPS Measurements as orange circles
plt.plot(epochs, gps_alt, 'o', color = 'orange', markersize =10,
          label = 'GPS noisy measurements')


# plot3: KF estimates as blue line
plt.plot(epochs, kf_alt, 's-', color = 'green', linewidth = 2.5,markersize= 8,
         label = 'Kalman Filter estimates')

# plot4: Uncertainty band around KF estimates
plt.fill_between(epochs, kf_alt - std, kf_alt +std, alpha = 0.2, color = 'green',
                 label = 'KF uncertainty (±1 std)')

#SECTION5: Style the plot
plt.xlabel('Epoch', fontsize = 12)
plt.ylabel('Altitude (m)', fontsize = 12)
plt.title('1D Kalman Filter Altitude Estimation', fontsize = 14, fontweight = 'bold')
plt.grid(True, linestyle = '--', alpha = 0.5)
plt.legend()
plt.xticks(epochs)  # Show epoch numbers on x-axis

#SECTION6: SAve and display
plt.savefig('../data/kf_1d_plot.png', dpi = 300, bbox_inches = 'tight')  # Save the plot as a PNG file
print('Plot saved to ../data/kf_1d_plot.png')
plt.show()

