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
# Read the sample precipiation field that is defined in two-dimenstions on a regular 0.25deg lat/long grid  -----------------------------------------------------------------------------------------------------------------------

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
# Ganerate some missing data in the rectangular region 10S-10N, 0E-20E
# ------------------------------------------------------------------------------------------------------------------------

lat_min = -10
lat_max = 10
lon_min = 0
lon_max = 20
f[np.logical_and.reduce((lat > lat_min, lat < lat_max, lon > lon_min, lon < lon_max))] = np.nan

# ------------------------------------------------------------------------------------------------------------------------
# Preparing the the f and area_size arrays 
# ------------------------------------------------------------------------------------------------------------------------

# in order for the missing data to be handled properly the area size of the coresponding points needs to be set to 0
area_size[np.isnan(f)] = 0

# the library does not allow the f array to be a masked array or to contain non-numeric values (e.g., not a number - nan), 
# thus we need to set the missing data points values to some arbitrary numeric value (the value itself is not important and will 
# not affect the calculation of the smoothed values, as long as the value is numeric). 
# Here we set the value of the missing data point values to 0.
f_numeric_only = f.copy()
f_numeric_only[np.isnan(f)] = 0

# ------------------------------------------------------------------------------------------------------------------------
# Smooth the field using a k-d tree 
# ------------------------------------------------------------------------------------------------------------------------

# construct a k-d tree 
kdtree_pointer = construct_KdTree(lat, lon)

# smooth the field using a k-d tree and smoothing kernel radius of R = 200 km
smoothing_kernel_radius_in_metres = 200*1000
f_smoothed = smooth_field_using_KdTree(smoothing_kernel_radius_in_metres, kdtree_pointer, lat, lon, area_size, f_numeric_only)

# free k-d tree memory
free_KdTree_memory(kdtree_pointer)

# the library will allways return 0 for the missing data points so we set the values of these points back to nan
f_smoothed[np.isnan(f)] = np.nan


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
ax.add_patch(matplotlib.patches.Rectangle( (lon_min, lat_min), (lon_max-lon_min), (lat_max - lat_min),color = "silver", transform=cartopy.crs.PlateCarree()))
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
ax.add_patch(matplotlib.patches.Rectangle( (lon_min, lat_min), (lon_max-lon_min), (lat_max - lat_min),color = "silver", transform=cartopy.crs.PlateCarree()))
ax.coastlines(resolution='110m', color='grey', linestyle='-', alpha=1)
matplotlib.pyplot.show()
matplotlib.pyplot.close()

