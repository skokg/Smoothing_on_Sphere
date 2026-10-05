"""Generate overlap caches and compare individual and multiple-field smoothing."""

from pathlib import Path

import numpy as np
from smoothing_on_sphere_library import (
    generate_overlap_cache_data_and_write_to_disk,
    read_overlap_cache_data_from_binary_file,
    free_overlap_cache_data,
    smooth_field_using_overlap_detection,
    smooth_multiple_fields_simultaneously_using_overlap_detection,
)

from netCDF4 import Dataset

# Locate the example inputs and cache files relative to this script.
data_folder = Path(__file__).resolve().parent

# Load the first precipitation field and the grid shared by both fields.
with Dataset(str(data_folder / "example_data" / "example_field1.nc"), "r") as nc_file_id:
    nc_file_id.set_auto_mask(False)
    lon_netcdf = np.asarray(nc_file_id.variables["lon"][:], dtype=np.float64)
    lat_netcdf = np.asarray(nc_file_id.variables["lat"][:], dtype=np.float64)
    f1_netcdf = np.asarray(nc_file_id.variables["precipitation"][:], dtype=np.float64)

# Load the second field on the same latitude/longitude grid.
with Dataset(str(data_folder / "example_data" / "example_field2.nc"), "r") as nc_file_id:
    nc_file_id.set_auto_mask(False)
    f2_netcdf = np.asarray(nc_file_id.variables["precipitation"][:], dtype=np.float64)

# Flatten the two-dimensional fields into the point order used by the library.
f1 = np.asarray(f1_netcdf, dtype=np.float64).reshape(-1)
f2 = np.asarray(f2_netcdf, dtype=np.float64).reshape(-1)
# Longitude varies fastest; repeat each latitude for its full grid row.
lon = np.tile(lon_netcdf,f1_netcdf.shape[0])
lat = np.repeat(lat_netcdf,f1_netcdf.shape[1])

# Approximate areas for the regular 0.25-degree grid, accounting for latitude.
Earth_radius= 6371.0*1000.0
dlat = 0.25 # resolution of lat/lon grid of the input field
area_size = np.deg2rad(dlat)*Earth_radius*np.deg2rad(dlat)*Earth_radius*np.cos(np.deg2rad(lat))
area_size[area_size < 0] = 0  # fix small negative values that occur at the poles due to the float rounding error

# One row per field; both fields share the same one-dimensional area vector.
f_multiple_fields = np.asarray([f1, f2], dtype=np.float64)

# Configure native threading 
number_of_threads = 10
show_output_log = True  # Set False to hide cache paths, timings and progress; radius warnings always print.
# Create the output directory for reusable geometry caches.
overlap_cache_data_folder = data_folder / "overlap_cache_data"
overlap_cache_data_folder.mkdir(parents=True, exist_ok=True)

# Generate and save the 100/200 km cache data, replacing existing files.
radii_in_meters = [100_000, 200_000]
# Passing both radii together shares geometry preparation; field values are not cached.
generate_overlap_cache_data_and_write_to_disk(
    lat, lon, radii_in_meters, overlap_cache_data_folder, number_of_threads, show_output_log)

# Select only the 200 km cache for smoothing
smoothing_kernel_radius_in_metres = 200_000

# Load the 200 km cache once and use it for both variants.
overlap_cache_data = read_overlap_cache_data_from_binary_file(
    overlap_cache_data_folder, smoothing_kernel_radius_in_metres, show_output_log)
try:
    # Variant 1: two independent single-field calls.
    f1_separately_smoothed = smooth_field_using_overlap_detection(
        area_size, f_multiple_fields[0], overlap_cache_data, number_of_threads, show_output_log)
    # Reuse the loaded cache for field 2 without another neighborhood search.
    f2_separately_smoothed = smooth_field_using_overlap_detection(
        area_size, f_multiple_fields[1], overlap_cache_data, number_of_threads, show_output_log)
    # Arrange individual results in the same row layout as the multiple-field result.
    f_separately_smoothed = np.stack([f1_separately_smoothed, f2_separately_smoothed])

    # Variant 2: one call for both fields,
    # sharing area denominators and buffers; each field still uses OpenMP.
    f_smoothed = smooth_multiple_fields_simultaneously_using_overlap_detection(
        area_size, f_multiple_fields, overlap_cache_data, number_of_threads, show_output_log)
finally:
    # Release native cache memory even if a smoothing call raises an exception.
    free_overlap_cache_data(overlap_cache_data)
    overlap_cache_data = None

# Visualization after both variants finish: original, individual and batch.
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# Use the same precipitation scale for every panel.
cmap = colors.LinearSegmentedColormap.from_list(
    "precipitation", ["white", (0.7, 0.7, 1.0), "blue"], 15)
fig = plt.figure(figsize=(15, 9))
# One row per field; columns show the input and the two smoothing variants.
for index, original in enumerate((f1_netcdf, f2_netcdf)):
    # Restore the original two-dimensional shape for plotting.
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
        # Supply geographic data coordinates to the orthographic map.
        image = ax.pcolormesh(
            lon_netcdf, lat_netcdf,
            values,
            transform=ccrs.PlateCarree(), shading="auto", cmap=cmap,
            norm=colors.LogNorm(vmin=0.01, vmax=100))
        ax.coastlines(resolution="110m", color="grey")
        fig.colorbar(image, ax=ax, extend="both", shrink=0.6).set_label("Precipitation (mm/6h)")
# Display the plots only; the caches are saved, but smoothed fields are not.
fig.tight_layout()
plt.show()
plt.close(fig)



