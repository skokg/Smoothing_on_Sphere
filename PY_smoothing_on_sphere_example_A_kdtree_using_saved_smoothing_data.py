# ------------------------------
# Import the needed libraries
# ------------------------------

# Import the PAD python library
# The PAD python library file (PY_smoothing_on_sphere_library.py) as well as the PAD C++ shared library file (smoothing_on_sphere_Cxx_shared_library.so) 
# need to be in the same folder as the script using them
from PY_smoothing_on_sphere_library import *

# Import netcdf library 
from netCDF4 import Dataset

# ------------------------------------------------------------------------------------------------------------------------
# Read the sample precipiation field that is defined in two-dimenstions on a regular 0.25deg lat/long grid 
# -----------------------------------------------------------------------------------------------------------------------

nc_file_id = Dataset("PY_smoothing_on_sphere_example_field1.nc", 'r') 
lon_netcdf = nc_file_id.variables["lon"][:].data
lat_netcdf = nc_file_id.variables["lat"][:].data
f_netcdf = nc_file_id.variables["precipitation"][:].data
nc_file_id.close()
 
# ------------------------------------------------------------------------------------------------------------------------
# convert data to lists of lats, lons, and field values for all grid points 
# ------------------------------------------------------------------------------------------------------------------------

f =np.reshape(f_netcdf,(-1))
lon = np.tile(lon_netcdf,f_netcdf.shape[0])
lat = np.repeat(lat_netcdf,f_netcdf.shape[1])

# ------------------------------------------------------------------------------------------------------------------------
# Calculate area size data
# ------------------------------------------------------------------------------------------------------------------------
Earth_radius= 6371.0*1000.0
dlat = 0.25 # resolution of lat/lon grid of the input field
area_size =  np.deg2rad(dlat)*Earth_radius*np.deg2rad(dlat)*Earth_radius*np.cos(np.deg2rad(lat))
area_size [area_size < 0] = 0  # fix small negative values that occur at the poles due to the float rounding error

# ------------------------------------------------------------------------------------------------------------------------
# Generate the kdtree-based smoothing data for smoothing kernel radiuses of 200 km and write it to the disk
# ------------------------------------------------------------------------------------------------------------------------

smoothing_kernel_radius_in_metres = 200*1000

# set the output folder for the smoothing data files
kdtree_smoothing_data_folder = bytes("smoothing_data_kdtree/", encoding='utf8')

# create the folder if it does not exist
os.makedirs(kdtree_smoothing_data_folder, exist_ok = True)

# Generate the smoothing data and write it to the disk
generate_smoothing_data_for_the_kdtree_based_approach_and_write_it_to_the_disk(lat, lon, smoothing_kernel_radius_in_metres, kdtree_smoothing_data_folder)

# ------------------------------------------------------------------------------------------------------------------------
# Use the kdtree-based smoothing data to calculate the smoothed field for the 200 km smoothing kernel radius
# ------------------------------------------------------------------------------------------------------------------------

# read the kdtree-based smoothing data from the disk into the memory
kdtree_smoothing_data_pointer = read_smoothing_data_for_the_kdtree_based_approach_from_binary_file_ctypes(kdtree_smoothing_data_folder, smoothing_kernel_radius_in_metres)

# calculate the smoothed field from the kdtree-based smoothing data
f_smoothed = smooth_field_using_smoothing_data_for_the_kdtree_approach(area_size, f, kdtree_smoothing_data_pointer)

# free the kdtree-based smoothing data memory
free_smoothing_data_memory_kdtree(kdtree_smoothing_data_pointer)
smoothing_data_pointer = None


# ------------------------------------------------------------------------------------------------------------------------
# Visualization 
# ------------------------------------------------------------------------------------------------------------------------

# Import matplotlib library 
import matplotlib
import matplotlib.pyplot

# Import cartopy library 
import cartopy

# reshape the smoothed field into a 2D shape of the original netcdf field
f_smoothed_2D = np.reshape(f_smoothed,f_netcdf.shape)

fig = matplotlib.pyplot.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=cartopy.crs.Orthographic(central_longitude=0, central_latitude=0, globe=None))
ax.set_global()
matplotlib.pyplot.title("Original precipiation field")
cmap_b = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',["white",(0.7, 0.7, 1.0), "blue"],15)
img = matplotlib.pyplot.imshow(f_netcdf, transform=cartopy.crs.PlateCarree(), interpolation='nearest', origin='lower', extent=(0, 360, -90, 90), cmap=cmap_b, norm = matplotlib.colors.LogNorm(vmin=0.01, vmax=100))
cb = fig.colorbar(img, extend='both', shrink=0.5)
cb.set_label('Precipitation (mm/6h)')
ax.coastlines(resolution='110m', color='grey', linestyle='-', alpha=1)
matplotlib.pyplot.show()
matplotlib.pyplot.close()

fig = matplotlib.pyplot.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=cartopy.crs.Orthographic(central_longitude=0, central_latitude=0, globe=None))
ax.set_global()
matplotlib.pyplot.title("Smoothed precipiation field")
cmap_b = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',["white",(0.7, 0.7, 1.0), "blue"],15)
img = matplotlib.pyplot.imshow(f_smoothed_2D, transform=cartopy.crs.PlateCarree(), interpolation='nearest', origin='lower', extent=(0, 360, -90, 90), cmap=cmap_b, norm = matplotlib.colors.LogNorm(vmin=0.01, vmax=100))
cb = fig.colorbar(img, extend='both', shrink=0.5)
cb.set_label('Precipitation (mm/6h)')
ax.coastlines(resolution='110m', color='grey', linestyle='-', alpha=1)
matplotlib.pyplot.show()
matplotlib.pyplot.close()

