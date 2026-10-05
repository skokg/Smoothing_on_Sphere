"""Compare individual and multiple-field direct KD-tree smoothing."""

from pathlib import Path

import numpy as np
from smoothing_on_sphere_library import (
    construct_KdTree,
    free_KdTree_memory,
    smooth_field_using_KdTree,
    smooth_multiple_fields_simultaneously_using_KdTree,
)

from netCDF4 import Dataset

# Locate the sample files relative to this script.
data_folder = Path(__file__).resolve().parent

# Load the first field and the regular grid shared by both fields.
with Dataset(str(data_folder / "example_data" / "example_field1.nc"), "r") as nc_file_id:
    nc_file_id.set_auto_mask(False)
    lon_netcdf = np.asarray(nc_file_id.variables["lon"][:], dtype=np.float64)
    lat_netcdf = np.asarray(nc_file_id.variables["lat"][:], dtype=np.float64)
    f1_netcdf = np.asarray(nc_file_id.variables["precipitation"][:], dtype=np.float64)

# Load the second field on the same latitude/longitude grid.
with Dataset(str(data_folder / "example_data" / "example_field2.nc"), "r") as nc_file_id:
    nc_file_id.set_auto_mask(False)
    f2_netcdf = np.asarray(nc_file_id.variables["precipitation"][:], dtype=np.float64)

# Flatten each (latitude, longitude) field into the point order expected by the library.
f1 = np.asarray(f1_netcdf, dtype=np.float64).reshape(-1)
f2 = np.asarray(f2_netcdf, dtype=np.float64).reshape(-1)
# Match that order: longitude varies fastest, with one latitude per grid row.
lon = np.tile(lon_netcdf,f1_netcdf.shape[0])
lat = np.repeat(lat_netcdf,f1_netcdf.shape[1])

# Approximate each 0.25-degree grid cell's area; cells shrink toward the poles.
Earth_radius= 6371.0*1000.0
dlat = 0.25 # resolution of lat/lon grid of the input field
area_size = np.deg2rad(dlat)*Earth_radius*np.deg2rad(dlat)*Earth_radius*np.cos(np.deg2rad(lat))
area_size[area_size < 0] = 0  # fix small negative values that occur at the poles due to the float rounding error

# Put the fields in separate rows of one (number_of_fields, number_of_points) array.
# Both variants use the same grid and the same one-dimensional area vector.
f_multiple_fields = np.asarray([f1, f2], dtype=np.float64)
f1 = f_multiple_fields[0]
f2 = f_multiple_fields[1]

# Configure the native OpenMP team, console output and 200 km smoothing radius.
number_of_threads = 10
show_output_log = True  # Set False to hide native smoothing timings and progress.
smoothing_kernel_radius_in_metres = 200*1000
# Build the spatial tree once and reuse it for all three smoothing calls below.
kdtree = construct_KdTree(lat, lon, number_of_threads)
try:
    # First smooth each field individually using the same KD-tree.
    f1_separately_smoothed = smooth_field_using_KdTree(
        smoothing_kernel_radius_in_metres, kdtree, lat, lon,
        area_size, f1, number_of_threads, show_output_log)
    # The second call smooths only field 2; each call prints its own native timing.
    f2_separately_smoothed = smooth_field_using_KdTree(
        smoothing_kernel_radius_in_metres, kdtree, lat, lon,
        area_size, f2, number_of_threads, show_output_log)
    # Keep the individual results in the same row layout as the multiple-field result.
    f_separately_smoothed = np.stack([f1_separately_smoothed, f2_separately_smoothed])

    # Then smooth both fields in one call, sharing neighborhood searches and area sums, so it is faster than smoothing each field separately
    f_smoothed = smooth_multiple_fields_simultaneously_using_KdTree(
        smoothing_kernel_radius_in_metres, kdtree, lat, lon,
        area_size, f_multiple_fields, number_of_threads,
        show_output_log)
finally:
    # Release the native tree even if a smoothing call raises an exception.
    free_KdTree_memory(kdtree)
    kdtree = None

# Visualization: original, individually smoothed and simultaneously smoothed.
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# Use the same precipitation color scale for every panel.
cmap = colors.LinearSegmentedColormap.from_list(
    "precipitation", ["white", (0.7, 0.7, 1.0), "blue"], 15)
fig = plt.figure(figsize=(15, 9))
# One row per field; columns show the input and the two smoothing variants.
for index, original in enumerate((f1_netcdf, f2_netcdf)):
    # Restore the original two-dimensional grid shape for plotting.
    panels = (
        ("Original", original),
        ("Single-field smoothing", f_separately_smoothed[index].reshape(original.shape)),
        ("Multi-field smoothing", f_smoothed[index].reshape(original.shape)),
    )
    for column, (label, values) in enumerate(panels):
        ax = fig.add_subplot(2, 3, index * 3 + column + 1,
                            projection=ccrs.Orthographic(central_longitude=0, central_latitude=0))
        ax.set_global()
        ax.set_title(f"Field {index + 1}: {label}")
        # Data coordinates are geographic even though the map uses an orthographic view.
        image = ax.pcolormesh(
            lon_netcdf, lat_netcdf,
            values,
            transform=ccrs.PlateCarree(), shading="auto", cmap=cmap,
            norm=colors.LogNorm(vmin=0.01, vmax=100))
        ax.coastlines(resolution="110m", color="grey")
        fig.colorbar(image, ax=ax, extend="both", shrink=0.6).set_label("Precipitation (mm/6h)")
# Display the plots only; no smoothed fields are saved to disk.
fig.tight_layout()
plt.show()
plt.close(fig)

