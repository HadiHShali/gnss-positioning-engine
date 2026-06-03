# ====================================================================#
# plot_gps_map.py - Read GPS CSV and create an interactive folium map #
# ====================================================================#

#run from project root: python scripts/plot_gps_map.py

#SECTION0: Headers or libraries
import pandas as pd
import math
import folium

# SECTION1: Load the data ---------------------------------------------------------
print('Loading GPS track...')
df =pd.read_csv('data/gps_track.csv')
print(f' {len(df)} GPS positions loaded.')
print(df.head())

#SECTION2: Prepare Coordinates ----------------------------------------------------
#Folium needs a list of [lat,lon] pairs
coordinates = list(zip(df['lat'], df['lon']))

#find the center of the track (for using in map centering)
center_lat = df['lat'].mean()
center_lon = df['lon'].mean()

#SECTION3: Compute statistics -----------------------------------------------------
def haversine_m(lat1, lon1, lat2, lon2):
	"""Distance between two GPS points in meters"""
	R = 6371000 #Earth Radius in meters
	lat1, lon1, lat2, lon2 = map(math.radians, [lat1, lon1, lat2, lon2])
	dlat = lat2 - lat1
	dlon = lon2 - lon1
	a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
	c = 2 * math.asin(math.sqrt(a))
	return R * c

total_distance = 0
for i in range(1, len(df)):
	total_distance += haversine_m(df['lat'].iloc[i-1], df['lon'].iloc[i-1],
									df['lat'].iloc[i], df['lon'].iloc[i])

#time duration
duration_sec  = len(df)
#Average speed in m/s
avg_speed_ms = total_distance / duration_sec if duration_sec>0 else 0

#print statistics
print(f'\nStatistics:')
print(f'Total Distance: {total_distance:.1f} in meter')
print(f'Duration: {duration_sec} in seconds or {duration_sec/60:.1f} in minutes')
print(f'Average Speed: {avg_speed_ms:.2f} in meter/sec or {avg_speed_ms*3.6:.1f} in km/h ')

#Section4: Create the map -------------------------------------------------------------
print('\nCreating map...')
m = folium.Map(
	location=[center_lat, center_lon],
	zoom_start = 15,
	tiles = 'OpenStreetMap' #default free ,ap tiles
	)

#SECTION5: Add the GPS Track (PolyLine) ------------------------------------------------
folium.PolyLine(
	coordinates,
	color = 'blue',
	weight = 4,
	opacity = 0.8,
	tooltip = f'Walk: {total_distance:.0f} m').add_to(m)


#SECTION6: Add Start and End Markers ----------------------------------------------------------------
#start: green Markers
#<b>text</b>: Makes text bold
#<br>: Inserts a line break (new line)
folium.Marker([df['lat'].iloc[0], df['lon'].iloc[0]],													
				popup=f"<b>START</b><br>Time: {df['time'].iloc[0]}<br>Alt: {df['alt'].iloc[0]:.1f} m",  
				icon = folium.Icon(color='green', icon = 'play')).add_to(m)
# END: red marker — offset by a tiny amount to avoid overlap
end_lat = df['lat'].iloc[-1] + 0.00005   # ~5 meters north
end_lon = df['lon'].iloc[-1] + 0.00005   # ~5 meters east
#End: red Markers
folium.Marker([end_lat, end_lon],
				popup=f"<b>End</b><br>Time: {df['time'].iloc[-1]}<br>Alt: {df['alt'].iloc[-1]:.1f} m",
				icon = folium.Icon(color='red', icon = 'stop')).add_to(m)


#SECTION7: Add a Statistics Box. ----------------------------------------------------------------------------
# We create the box using <div>..</div> command
# To style the box we use style="..." which is a CSS code to make the style of the box. Each CSS rule is a property: value; pair, like position: fixed;
# The <div> — a container:  <div> ... content goes here ... </div>


stats_html = f'''
<div style="position: fixed; top: 10px; right: 10px;
            background: white; padding: 10px; border: 2px solid #1B3A6B;
            border-radius: 5px; font-family: Arial; z-index: 1000;">
 <b>Walk Statistics</b><br>
Distance: {total_distance:.0f} m<br>
Duration: {duration_sec/60:.1f} min<br>
Avg speed: {avg_speed_ms:.2f} m/s<br>
Points: {len(df)}
</div>
'''

m.get_root().html.add_child(folium.Element(stats_html))


#SECTION8: SAve----------------------------------------------------------------------------------
output_path = 'data/gps_track_map.html'
m.save(output_path)
print(f'\nMap saved to: {output_path}')
print(f'Open it in your browser to see the interactive map!!')
