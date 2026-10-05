"""Smooth a precipitation field with missing data using the direct KD-tree path."""

from pathlib import Path

import numpy as np
from netCDF4 import Dataset
from smoothing_on_sphere_library import (
    construct_KdTree,
    free_KdTree_memory,
    smooth_field_using_KdTree,
)

# Locate the example inputs and cache files relative to this script.
data_folder = Path(__file__).resolve().parent
field_filename = data_folder / "example_data" / "example_field1.nc"
# Configure the direct smoothing call and its native console output.
number_of_threads = 10
show_output_log = True  # Set False to hide native timings and progress.
smoothing_kernel_radius_in_metres = 200_000

# Read the sample field and coordinates as ordinary NumPy arrays.
with Dataset(str(field_filename), "r") as dataset:
    dataset.set_auto_mask(False)
    lon_netcdf = np.asarray(dataset.variables["lon"][:], dtype=np.float64)
    lat_netcdf = np.asarray(dataset.variables["lat"][:], dtype=np.float64)
    original_field = np.asarray(dataset.variables["precipitation"][:], dtype=np.float64)

# Keep the grid shape so flattened results can be restored for plotting.
expected_shape = (lat_netcdf.size, lon_netcdf.size)

# Match flattened field order: longitude varies fastest within each latitude row.
lon = np.tile(lon_netcdf, lat_netcdf.size)
lat = np.repeat(lat_netcdf, lon_netcdf.size)
# Copy the field before replacing values in the artificial missing-data region.
f_numeric_only = np.asarray(original_field, dtype=np.float64).reshape(-1).copy()

# Add an artificial missing-data region: 10S–10N, 0E–20E.
lat_min, lat_max = -10, 10
lon_min, lon_max = 0, 20
# Express longitude in 0–360 degrees so the region works with either input convention.
longitude_degrees_east = np.mod(lon, 360.0)
missing = (
    (lat > lat_min) & (lat < lat_max)
    & (longitude_degrees_east > lon_min) & (longitude_degrees_east < lon_max))

# Approximate areas for this regular 0.25-degree latitude/longitude example grid.
earth_radius = 6_371_000.0
grid_spacing_in_degrees = 0.25
area_size = (
    (np.deg2rad(grid_spacing_in_degrees) * earth_radius) ** 2
    * np.cos(np.deg2rad(lat)))
# Clamp tiny negative areas caused by floating-point rounding near the poles.
np.maximum(area_size, 0.0, out=area_size)

# The native API requires finite numeric arrays. Zero area excludes a point
# from both the weighted numerator and the neighborhood denominator.
area_size[missing] = 0.0
f_numeric_only[missing] = 0.0

# Construct the spatial tree once from all grid points.
kdtree = construct_KdTree(lat, lon, number_of_threads)
try:
    # Zero-area points are excluded from neighborhood contributions.
    f_smoothed = smooth_field_using_KdTree(
        smoothing_kernel_radius_in_metres, kdtree, lat, lon,
        area_size, f_numeric_only, number_of_threads,
        show_output_log)
finally:
    # Always release the native tree, including when smoothing fails.
    free_KdTree_memory(kdtree)
    kdtree = None

# Show excluded points as NaN without changing the finite smoothing inputs.
display_original = np.where(missing, np.nan, f_numeric_only).reshape(expected_shape)
display_smoothed = np.where(missing, np.nan, f_smoothed).reshape(expected_shape)

import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import cartopy.crs as ccrs

# Use one color scale for the input and smoothed field; NaNs appear white.
cmap = colors.LinearSegmentedColormap.from_list(
    "precipitation", ["white", (0.7, 0.7, 1.0), "blue"], 15)
cmap.set_bad("white")
fig = plt.figure(figsize=(14, 7))
# Show the artificial missing-data region before and after smoothing.
for column, (title, values) in enumerate((
        ("Input precipitation with missing data", display_original),
        ("Smoothed precipitation (200 km)", display_smoothed))):
    ax = fig.add_subplot(
        1, 2, column + 1,
        projection=ccrs.Orthographic(central_longitude=0, central_latitude=0))
    ax.set_global()
    ax.set_title(title)
    image = ax.pcolormesh(
        lon_netcdf, lat_netcdf, values,
        transform=ccrs.PlateCarree(), shading="auto", cmap=cmap,
        norm=colors.LogNorm(vmin=0.01, vmax=100))
    # Outline the same excluded region on both maps.
    ax.add_patch(Rectangle(
        (lon_min, lat_min), lon_max - lon_min, lat_max - lat_min,
        fill=False, edgecolor="grey", transform=ccrs.PlateCarree()))
    ax.coastlines(resolution="110m", color="grey")
    fig.colorbar(image, ax=ax, extend="both", shrink=0.6).set_label("Precipitation (mm/6h)")
# Display the figures without writing smoothed fields to disk.
fig.tight_layout()
plt.show()
plt.close(fig)

