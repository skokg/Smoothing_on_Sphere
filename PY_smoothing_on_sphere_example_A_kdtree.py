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
# Smooth the field using a k-d tree 
# ------------------------------------------------------------------------------------------------------------------------

# construct a k-d tree 
kdtree_pointer = construct_KdTree(lat, lon)

# smooth the field using a k-d tree and smoothing kernel radius of R = 200 km
smoothing_kernel_radius_in_metres = 200*1000
f_smoothed = smooth_field_using_KdTree(smoothing_kernel_radius_in_metres, kdtree_pointer, lat, lon, area_size, f)

# free k-d tree memory
free_KdTree_memory(kdtree_pointer)
kdtree_pointer = None

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

