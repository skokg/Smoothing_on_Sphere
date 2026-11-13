# Smoothing on Sphere

Smooth fields defined on a sphere using area-weighted spherical neighbourhoods.
The package provides direct KD-tree smoothing, reusable KD-tree caches and
overlap-detection caches through a NumPy interface to C++ with OpenMP.

The distribution and Python import name are `smoothing_on_sphere`.
GitHub repository: [skokg/Smoothing_on_Sphere](https://github.com/skokg/Smoothing_on_Sphere).

## Installation

Python 3.10 or newer is required. Install the package with:

```sh
pip install smoothing_on_sphere
```

To install from source:

```sh
pip install --no-binary=smoothing_on_sphere smoothing_on_sphere
```

Compiling from source requires a C++ compiler and OpenMP development/runtime libraries: GCC or Clang on Linux, MSVC Build Tools with `clang-cl` and LLVM OpenMP on Windows, or Xcode Command Line Tools with LLVM OpenMP on macOS.
On Windows, set `LIBOMP_PREFIX` to the LLVM OpenMP prefix for building and importing source installations.

To install with conda:

```sh
conda install -c conda-forge smoothing_on_sphere
```

Alternatively, download the source from the
[GitHub repository](https://github.com/skokg/Smoothing_on_Sphere) and run `make` in the
project's root directory on Linux. This requires GNU Make, a C++ compiler and OpenMP.
Running `make` creates `libsmoothing_on_sphere.so` in the `src` folder.
Copy `smoothing_on_sphere.py` and `libsmoothing_on_sphere.so` from `src` into the folder
containing your Python script. NumPy must also be installed.

## Usage

Supply one-dimensional latitude/longitude arrays in degrees, cell areas (or
consistent relative area weights), and finite field values in matching point
order. Smoothing radii are in metres. For missing data, use zero area and a finite
field placeholder; missing target cells receive zero output. Do not pass NaNs.

This example builds a reproducible global field from 200 random rectangles,
each adding 1 wherever it covers the grid, so overlapping rectangles add up.
Direct KD-tree and overlap-cache smoothing use radii of 200 and 1000 km.

```python
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import smoothing_on_sphere as smoothing

# Global 1-degree cell centres: 180 latitude rows and 360 longitude columns.
lon_grid, lat_grid = np.meshgrid(np.arange(-179.5, 180., 1.), np.arange(-89.5, 90., 1.))
field = np.zeros_like(lat_grid)
rng = np.random.default_rng(42)  # Fixed seed: generate the same field each run.

# Add 200 rectangles, each 5-30 degrees high and 5-60 degrees wide.
for height, width in rng.integers([5, 5], [31, 61], size=(200, 2)):
    field[tuple(slice(start, start + size) for start, size in zip(rng.integers(0, np.array(field.shape) - [height, width] + 1), (height, width)))] += 1  # Overlaps accumulate.

# Flatten all inputs in the same order; equal-angle cell areas vary as cos(latitude).
lat, lon = lat_grid.ravel(), lon_grid.ravel()
area = np.cos(np.deg2rad(lat))  # Relative area weights are sufficient.
values = field.ravel()
threads = 4  # Adjust to the number of CPU threads you want to use.
direct_smoothed_fields, overlap_smoothed_fields = [field], [field]  # First column shows the same original field.

# Reuse one tree for both radii; each calculation starts from the original field.
tree = smoothing.construct_KdTree(lat, lon, threads)
try:
    for radius_km in (200, 1000):
        radius = radius_km * 1000.  # The smoothing functions take metres.
        # Direct KD-tree smoothing needs no cache, but repeats neighbourhood searches and is typically slower than cached smoothing.
        direct_smoothed_fields.append(smoothing.smooth_field_using_KdTree(radius, tree, lat, lon, area, values, threads).reshape(field.shape))

        # Overlap smoothing needs a cache: preparation takes time and memory, but reuse makes smoothing faster.
        # Construct only this radius's cache in memory; reuse it for other fields on the same grid.
        cache = smoothing.generate_overlap_cache_data_in_memory(lat, lon, radius, threads)
        try:
            overlap_smoothed_fields.append(smoothing.smooth_field_using_overlap_detection(area, values, cache, threads).reshape(field.shape))
        finally:
            smoothing.free_overlap_cache_data(cache)  # Release before the next radius.
finally:
    smoothing.free_KdTree_memory(tree)

# Two rows (methods), three columns (original, 200 km, 1000 km), in one figure.
fig = plt.figure(figsize=(15, 9))
for row, (method, fields) in enumerate((("Direct KD-tree", direct_smoothed_fields), ("Overlap", overlap_smoothed_fields))):
    for column, (label, data) in enumerate(zip(("Original", "Smoothed (200 km)", "Smoothed (1000 km)"), fields)):
        ax = fig.add_subplot(2, 3, row * 3 + column + 1, projection=ccrs.Orthographic(central_longitude=0, central_latitude=0))
        ax.set_global()
        ax.set_title(f"{method}: {label}")
        # Data coordinates are geographic; all panels share the original field's scale.
        image = ax.pcolormesh(lon_grid, lat_grid, data, transform=ccrs.PlateCarree(), shading="auto", cmap="Blues", vmin=0, vmax=field.max())
        ax.coastlines(resolution="110m", color="grey")
        fig.colorbar(image, ax=ax, shrink=0.6).set_label("Value (overlap count)")
fig.tight_layout()
plt.show()  # Display only; no caches or smoothed fields are saved to disk.
plt.close(fig)
```

## Further examples

Download the examples from the [GitHub repository](https://github.com/skokg/Smoothing_on_Sphere) or obtain them from the source archive; the scripts are in `examples/` and their NetCDF sample data are in `examples/sample_data/`.

The examples additionally require netCDF4, Matplotlib and Cartopy.
They also demonstrate loading NetCDF fields, smoothing multiple fields simultaneously, using KD-tree caches, handling missing data, and optional cache saving/loading through commented commands.

* `Example_A_kdtree_direct.py`: direct single-field and multiple-field smoothing.
* `Example_B_kdtree_using_saved_cache.py`: one in-memory KD-tree cache and reuse.
* `Example_C_overlap_using_saved_cache.py`: one in-memory overlap cache and reuse.
* `Example_D_with_missing_data.py`: direct smoothing with missing data.

## Function descriptions

The following functions are available through `import smoothing_on_sphere as smoothing`.

### `smoothing.construct_KdTree(lat, lon, number_of_threads)`

**Inputs:** `lat` and `lon` are matching, nonempty one-dimensional NumPy arrays
of finite coordinates in degrees; masked arrays and NaNs are not accepted.
`number_of_threads` is a positive integer controlling native OpenMP threads.

**Output:** A KD-tree handle for these points; do not copy it.

**Functionality:** Builds the spatial tree used by direct smoothing. Reuse it for
different fields and radii on the same grid and in the same point order, then call
`free_KdTree_memory()`, preferably in a `finally` block.

### `smoothing.smooth_field_using_KdTree(smoothing_kernel_radius_in_metres, kdtree_pointer, lat, lon, area_size, f, number_of_threads, show_output_log=True)`

**Inputs:** A live tree handle and matching nonempty one-dimensional NumPy arrays
`lat`, `lon`, `area_size` and `f`, shaped `(number_of_points,)`. Coordinates are in
degrees and must match the tree's grid and point order; `area_size` contains cell
areas or nonnegative relative weights. Arrays must contain finite real numbers;
masked arrays and NaNs are not accepted. The radius is a finite, nonnegative scalar in metres,
capped at half the Earth's circumference. `number_of_threads` is a positive
OpenMP thread count. `show_output_log=True` enables timing/progress messages;
radius warnings always appear.

**Output:** A new `float64` array with the same shape as `f`, in the original point order.

**Functionality:** Searches the KD-tree for each target's neighbourhood and computes
its area-weighted mean. No smoothing cache is required; searches are repeated on
each call. Use zero area and a finite field placeholder for missing cells; their
output is zero. Other cells receive the area-weighted mean within the radius,
or zero if the neighbourhood has no positive total area.

### `smoothing.smooth_multiple_fields_simultaneously_using_KdTree(smoothing_kernel_radius_in_metres, kdtree_pointer, lat, lon, area_size, f, number_of_threads, show_output_log=True)`

**Inputs:** A live tree handle, matching one-dimensional `lat`, `lon` and `area_size`
arrays, and `f` shaped `(number_of_fields, number_of_points)`. Coordinates are in
degrees and must match the tree's grid and point order. All fields share the same
cell areas or nonnegative relative weights and zero-area missing-data mask.
Arrays must be nonempty NumPy arrays of finite real numbers; masked arrays and
NaNs are not accepted. The radius is a finite, nonnegative scalar in metres, capped at half the Earth's
circumference. `number_of_threads` is a positive OpenMP thread count.
`show_output_log=True` enables timing/progress messages; radius warnings always appear.

**Output:** A new `float64` array with the same two-dimensional shape as `f`,
preserving field and point order.

**Functionality:** Reuses each neighbourhood search and area denominator across
the fields, avoiding separate searches for every field. Use zero area and finite
field placeholders for missing cells; their output is zero in every field.
Other outputs are area-weighted neighbourhood means, or zero where the total area
is not positive.

### `smoothing.free_KdTree_memory(kdtree_pointer)`

**Input:** A live handle returned by `construct_KdTree()`, with no other call
using it concurrently.

**Output:** `None`.

**Functionality:** Releases the tree's native memory and invalidates the handle.
Call it in a `finally` block when finished; do not copy handles, use a freed handle,
or free it again.

### `smoothing.generate_kdtree_cache_data_in_memory(lat, lon, smoothing_kernel_radius_in_metres, number_of_threads, show_output_log=True)`

**Inputs:** Matching, nonempty one-dimensional NumPy arrays `lat` and `lon` of finite
coordinates in degrees; masked arrays and NaNs are not accepted. The radius is a finite,
nonnegative scalar in metres, capped at half the Earth's circumference.
`number_of_threads` is a positive OpenMP thread count. `show_output_log=True`
enables construction timing and memory messages; radius warnings always appear.

**Output:** A KD-tree cache handle for that grid, point order
and single radius; do not copy it.

**Functionality:** Prepares neighbourhood information entirely in
memory, without creating files. With logging enabled, successful construction
reports estimated cache memory in bytes and MiB. Reuse the cache for different fields and area weights on
the same grid and in the same point order, then call `free_kdtree_cache_data()`,
preferably in a `finally` block.

### `smoothing.smooth_field_using_kdtree_cache_data(area_size, f, kdtree_cache_data, number_of_threads, show_output_log=True)`

**Inputs:** Matching nonempty one-dimensional NumPy arrays `area_size` and `f`,
shaped `(number_of_points,)`, and a live KD-tree cache handle for the same grid
and point order. `area_size` contains cell areas or nonnegative relative weights.
Arrays must contain finite real numbers; masked arrays and NaNs are not accepted.
`number_of_threads`
is a positive OpenMP thread count; `show_output_log=True` enables timing messages.

**Output:** A new `float64` array with the same shape as `f`, in the original point order.

**Functionality:** Uses the cached neighbourhood information instead of repeating
radius searches. The grid, point order and radius are fixed by the cache. Use zero
area and a finite field placeholder for missing cells; their output is zero.
Other cells receive the area-weighted mean within the cached radius, or zero if
the neighbourhood has no positive total area.

### `smoothing.smooth_multiple_fields_simultaneously_using_kdtree_cache_data(area_size, fields, kdtree_cache_data, number_of_threads, show_output_log=True)`

**Inputs:** A one-dimensional `area_size` vector, `fields` shaped
`(number_of_fields, number_of_points)`, and a live KD-tree cache handle for the same
grid and point order. All fields share the same cell areas or nonnegative relative
weights and zero-area missing-data mask. Arrays must be nonempty NumPy arrays of
finite real numbers; masked arrays and NaNs are not accepted. `number_of_threads` is a positive OpenMP
thread count; `show_output_log=True` enables timing messages.

**Output:** A new `float64` array with the same shape as `fields`, preserving field
and point order.

**Functionality:** Reuses the cache and area calculations across the fields.
The cache fixes the radius. Use finite field placeholders at zero-area missing
cells; their output is zero in every field. Other outputs are area-weighted
neighbourhood means, or zero where the total area is not positive.

### `smoothing.save_kdtree_cache(kdtree_cache_data, filename, show_output_log=True)`

**Inputs:** A live KD-tree cache handle, a destination filename as a string or
`pathlib.Path`, and optional logging flag. The parent directory must already exist;
`show_output_log=True` prints the absolute saved filename.

**Output:** `None`; writes a binary cache file.

**Functionality:** Saves the existing cache, including its radius, without rebuilding
it. The filename is used exactly as supplied. An existing file is replaced only
after writing and closing succeed; saving leaves the cache available for reuse.

### `smoothing.load_kdtree_cache(filename, show_output_log=True)`

**Inputs:** A saved cache filename as a string or `pathlib.Path`, and optional
logging flag. `show_output_log=True` prints the cache path and loading time.

**Output:** A new KD-tree cache handle; do not copy it.
Release it with `free_kdtree_cache_data()`, preferably in a `finally` block.

**Functionality:** Loads and validates the binary cache, reading its radius and
point count from the contents. The filename need not encode the radius. The loaded
cache can be used for smoothing on the same grid and in the same point order, at
the stored radius, or saved again under another filename.

### `smoothing.free_kdtree_cache_data(kdtree_cache_data)`

**Input:** A live KD-tree cache handle returned by generation or loading, with no
other call using it concurrently.

**Output:** `None`.

**Functionality:** Releases the cache's native memory and invalidates the handle.
Call it in a `finally` block when finished; do not copy handles, use a freed handle,
or free it again.

### `smoothing.generate_overlap_cache_data_in_memory(lat, lon, smoothing_kernel_radius_in_metres, number_of_threads, show_output_log=True)`

**Inputs:** Matching, nonempty one-dimensional NumPy arrays `lat` and `lon` of finite
coordinates in degrees; masked arrays and NaNs are not accepted. The radius is a finite,
nonnegative scalar in metres, capped at half the Earth's circumference.
`number_of_threads` is a positive OpenMP thread count. `show_output_log=True`
enables construction progress, timing and memory messages; radius warnings always appear.

**Output:** An overlap cache handle for that grid, point order
and single radius; do not copy it.

**Functionality:** Prepares an in-memory cache for overlap smoothing, without creating files.
With logging enabled, successful construction reports estimated cache memory in bytes and MiB. Reuse it for different
fields and area weights on the same grid and in the same point order, then call
`free_overlap_cache_data()`, preferably in a `finally` block.

### `smoothing.smooth_field_using_overlap_detection(area_size, f, overlap_cache_data, number_of_threads, show_output_log=True)`

**Inputs:** Matching nonempty one-dimensional NumPy arrays `area_size` and `f`,
shaped `(number_of_points,)`, and a live overlap cache handle for the same grid
and point order. `area_size` contains cell areas or nonnegative relative weights.
Arrays must contain finite real numbers; masked arrays and NaNs are not accepted.
`number_of_threads`
is a positive OpenMP thread count; `show_output_log=True` enables timing messages.

**Output:** A new `float64` array with the same shape as `f`, in the original point order.

**Functionality:** Computes area-weighted means using cached additions and removals
between neighbourhoods, reusing sums along the smoothing sequence instead of
performing new radius searches. The cache fixes the radius. Use zero area and a
finite field placeholder for missing cells; their output is zero. Other cells
receive the area-weighted mean within the cached radius, or zero if the
neighbourhood has no positive total area.

### `smoothing.smooth_multiple_fields_simultaneously_using_overlap_detection(area_size, fields, overlap_cache_data, number_of_threads, show_output_log=True)`

**Inputs:** A one-dimensional `area_size` vector, `fields` shaped
`(number_of_fields, number_of_points)`, and a live overlap cache handle for the same
grid and point order. All fields share the same cell areas or nonnegative relative
weights and zero-area missing-data mask. Arrays must be nonempty NumPy arrays of
finite real numbers; masked arrays and NaNs are not accepted. `number_of_threads` is a positive OpenMP
thread count; `show_output_log=True` enables timing messages.

**Output:** A new `float64` array with the same shape as `fields`, preserving field
and point order.

**Functionality:** Reuses the overlap cache and area calculations across the fields.
The cache fixes the radius. Use finite field placeholders at zero-area missing
cells; their output is zero in every field. Other outputs are area-weighted
neighbourhood means, or zero where the total area is not positive.

### `smoothing.save_overlap_cache(overlap_cache_data, filename, show_output_log=True)`

**Inputs:** A live overlap cache handle, a destination filename as a string or
`pathlib.Path`, and optional logging flag. The parent directory must already exist;
`show_output_log=True` prints the absolute saved filename.

**Output:** `None`; writes a binary cache file.

**Functionality:** Saves the existing cache, including its radius, without rebuilding
it. The filename is used exactly as supplied. An existing file is replaced only
after writing and closing succeed; saving leaves the cache available for reuse.

### `smoothing.load_overlap_cache(filename, show_output_log=True)`

**Inputs:** A saved cache filename as a string or `pathlib.Path`, and optional
logging flag. `show_output_log=True` prints the cache path and loading time.

**Output:** A new overlap cache handle; do not copy it.
Release it with `free_overlap_cache_data()`, preferably in a `finally` block.

**Functionality:** Loads and validates the binary cache, reading its radius and
point count from the contents. The filename need not encode the radius. The loaded
cache can be used for smoothing on the same grid and in the same point order, at
the stored radius, or saved again under another filename.

### `smoothing.free_overlap_cache_data(overlap_cache_data)`

**Input:** A live overlap cache handle returned by generation or loading, with no
other call using it concurrently.

**Output:** `None`.

**Functionality:** Releases the cache's native memory and invalidates the handle.
Call it in a `finally` block when finished; do not copy handles, use a freed handle,
or free it again.

## Author and license

Gregor Skok, Faculty of Mathematics and Physics, University of Ljubljana, Slovenia.
Email: Gregor.Skok@fmf.uni-lj.si.

MIT License; see `LICENSE`, including the modified third-party KD-tree notice.
The release wheel configuration bundles LLVM OpenMP on Linux, Windows and macOS.
Local Linux source builds can also use GCC with the separately installed GNU runtime.
The LLVM OpenMP notices and license text are included under `licenses/` in source archives
and in the installed distribution's `.dist-info/licenses/` directory.
See `licenses/THIRD_PARTY_NOTICES.txt` for the applicable files and upstream sources.

## Reference

Skok, G. and Kosovelj, K. (2025): Smoothing and spatial verification of global fields,
Geoscientific Model Development, 18, 7417–7433.
[https://doi.org/10.5194/gmd-18-7417-2025](https://doi.org/10.5194/gmd-18-7417-2025).
