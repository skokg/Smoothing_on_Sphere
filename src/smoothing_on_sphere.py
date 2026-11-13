"""Sphere smoothing with radii capped in C++ at half the Earth's circumference.

Radius warnings always appear in native output. show_output_log controls
absolute cache filenames, timing and progress messages.
"""

import ctypes
import operator
import os
import sys
import tempfile

import numpy as np

# Windows wheel repair inserts the bundled OpenMP DLL search path here.
# delvewheel: patch

_library_folder = os.path.dirname(os.path.abspath(__file__))
_library_names = {
    "linux": "libsmoothing_on_sphere.so",
    "win32": "smoothing_on_sphere.dll",
    "darwin": "libsmoothing_on_sphere.dylib",
}
if sys.platform not in _library_names:
    raise ImportError("smoothing_on_sphere supports Linux, Windows and macOS.")
_library_path = os.path.join(_library_folder, _library_names[sys.platform])
_compile_advice = (
    "Install a compatible smoothing_on_sphere wheel, or build from source with "
    "a C++ compiler and OpenMP runtime. For manual Linux compilation, run `make` "
    "in the repository root, then use smoothing_on_sphere.py and "
    "libsmoothing_on_sphere.so from src together."
)
if sys.platform == "win32":
    _compile_advice += (
        " For Windows source installations, set LIBOMP_PREFIX at import time "
        "to the LLVM OpenMP prefix containing bin/libomp.dll."
    )
if not os.path.isfile(_library_path):
    raise ImportError(
        f"Smoothing shared library not found: {_library_path}. {_compile_advice}"
    )

# Repaired wheels arrange their bundled dependencies at the marker above.
# Try that first so an unrelated LIBOMP_PREFIX cannot interfere with a wheel.
# Source installations can fall back to an explicitly configured runtime.
# Keep the registration alive for any later native dependency loads.
_openmp_dll_directory = None
try:
    libc = ctypes.CDLL(_library_path)
except OSError as exc:
    _openmp_prefix = os.environ.get("LIBOMP_PREFIX") if sys.platform == "win32" else None
    if not _openmp_prefix:
        raise ImportError(
            f"Cannot load smoothing shared library: {_library_path}. {_compile_advice}"
        ) from exc
    _openmp_bin = os.path.abspath(os.path.join(_openmp_prefix, "bin"))
    try:
        _openmp_dll_directory = os.add_dll_directory(_openmp_bin)
        libc = ctypes.CDLL(_library_path)
    except OSError as retry_exc:
        if _openmp_dll_directory is not None:
            _openmp_dll_directory.close()
            _openmp_dll_directory = None
        raise ImportError(
            f"Cannot load smoothing shared library: {_library_path} using "
            f"the configured OpenMP DLL directory {_openmp_bin}. {_compile_advice}"
        ) from retry_exc

ND_POINTER_1D = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags=("C_CONTIGUOUS", "ALIGNED"))
ND_POINTER_2D = np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags=("C_CONTIGUOUS", "ALIGNED"))

try:
    libc.smoothing_wrapper_abi_version.argtypes = []
    libc.smoothing_wrapper_abi_version.restype = ctypes.c_int
except AttributeError as exc:
    raise ImportError("Reinstall smoothing_on_sphere: native wrapper ABI 3 required.") from exc
if libc.smoothing_wrapper_abi_version() != 3:
    raise ImportError("Incompatible smoothing native library: wrapper ABI 3 required.")

libc.smoothing_get_last_error.argtypes = []
libc.smoothing_get_last_error.restype = ctypes.c_char_p

def _raise_native_error(function_name):
    # Read immediately: another C-interface call, including cleanup, clears it.
    message = libc.smoothing_get_last_error()
    detail = message.decode("utf-8", errors="replace") if message else "Unknown native error"
    raise RuntimeError(function_name + ": " + detail)

def _check_status(result, function, arguments):
    if result != 0:
        _raise_native_error(function.__name__)
    return result

def _check_pointer(result, function, arguments):
    if not result:
        _raise_native_error(function.__name__)
    return result

def _bind(name, arguments, returns_handle=False):
    function = getattr(libc, name)
    function.argtypes = arguments
    function.restype = ctypes.c_void_p if returns_handle else ctypes.c_int
    function.errcheck = _check_pointer if returns_handle else _check_status

# Explicit free calls are still required. Do not free a handle concurrently
# with another use, copy its ownership metadata, or access its private pointer.
class _NativeHandle:
    __slots__ = ("_pointer", "_kind", "_number_of_points")

    def __init__(self, kind, number_of_points):
        self._pointer = None
        self._kind = kind
        self._number_of_points = number_of_points

    def __copy__(self):
        raise TypeError("Native handles must not be copied")

    def __deepcopy__(self, memo):
        raise TypeError("Native handles must not be copied")

def _new_handle(kind, number_of_points, function, *arguments):
    # Allocate ownership metadata before allocating the native resource.
    handle = _NativeHandle(kind, number_of_points)
    handle._pointer = function(*arguments)
    return handle

def _require_handle(handle, kind, number_of_points=None):
    if not isinstance(handle, _NativeHandle) or handle._kind != kind:
        raise TypeError("Expected a " + kind + " handle returned by this library")
    if handle._pointer is None:
        raise ValueError("The " + kind + " handle has already been freed")
    if (number_of_points is not None and handle._number_of_points is not None
            and number_of_points != handle._number_of_points):
        raise ValueError("Point count differs from the " + kind + " handle")
    return handle._pointer

def _positive_int(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(name + " must be an integer, not a boolean")
    try:
        value = operator.index(value)
    except TypeError:
        raise TypeError(name + " must be an integer") from None
    if not 1 <= value <= np.iinfo(np.int32).max:
        raise ValueError(name + " must be between 1 and INT32_MAX")
    return value

def _thread_count(number_of_threads):
    return _positive_int(number_of_threads, "number_of_threads")

def _check_array(array, array_name, dimensions):
    if type(array) is not np.ndarray:
        raise TypeError(array_name + " must be a numpy.ndarray, not a masked array")
    if array.ndim != dimensions:
        raise ValueError(array_name + " must have " + str(dimensions) + " dimensions")
    if array.size == 0:
        raise ValueError(array_name + " must not be empty")
    if array.dtype.kind not in "iuf":
        raise TypeError(array_name + " must contain real numeric values")
    if not np.all(np.isfinite(array)):
        raise ValueError(array_name + " must contain only finite values")
    return True

def check_array(array, array_name):
    return _check_array(array, array_name, 1)

def check_array_2D(array, array_name):
    return _check_array(array, array_name, 2)

def _float_array(array, array_name, dimensions=1):
    _check_array(array, array_name, dimensions)
    array = np.require(array, dtype=np.float64, requirements=["C", "A"])
    if not np.all(np.isfinite(array)):
        raise ValueError(array_name + " cannot be represented as finite float64 values")
    return array

def _coordinates(lat, lon):
    lat = _float_array(lat, "lat")
    lon = _float_array(lon, "lon")
    if lat.shape != lon.shape:
        raise ValueError("lat and lon must have the same shape")
    _positive_int(lat.size, "number_of_points")
    return lat, lon

def _fields(area_size, f, dimensions=1):
    area_size = _float_array(area_size, "area_size", dimensions)
    f = _float_array(f, "f", dimensions)
    if area_size.shape != f.shape:
        raise ValueError("area_size and f must have the same shape")
    _positive_int(f.shape[-1], "number_of_points")
    return area_size, f

def _radius(value):
    if np.ma.isMaskedArray(value) or np.ndim(value) != 0:
        raise ValueError("smoothing_kernel_radius_in_metres must be an unmasked scalar")
    value = float(value)
    if not np.isfinite(value) or value < 0:
        raise ValueError("smoothing_kernel_radius_in_metres must be finite and nonnegative")
    return value

def _filename(path):
    path = os.fsencode(os.fspath(path))
    if not path or b"\0" in path or os.path.basename(path) in (b"", b".", b".."):
        raise ValueError("Cache filename must be a nonempty file path without null bytes")
    return path

# Tree construction, direct smoothing and cleanup
_bind("construct_KdTree_ctypes", [ND_POINTER_1D, ND_POINTER_1D, ctypes.c_size_t, ctypes.c_int], True)

def construct_KdTree(lat, lon, number_of_threads):
    lat, lon = _coordinates(lat, lon)
    number_of_threads = _thread_count(number_of_threads)
    return _new_handle("tree", lat.size, libc.construct_KdTree_ctypes,
                       lat, lon, lat.size, number_of_threads)

_bind("free_KdTree_memory_ctypes", [ctypes.c_void_p])

def free_KdTree_memory(kdtree_pointer):
    pointer = _require_handle(kdtree_pointer, "tree")
    libc.free_KdTree_memory_ctypes(pointer)
    kdtree_pointer._pointer = None

_bind("smooth_field_using_KdTree_ctypes", [ctypes.c_double, ctypes.c_void_p,
      ND_POINTER_1D, ND_POINTER_1D, ND_POINTER_1D, ND_POINTER_1D,
      ctypes.c_size_t, ND_POINTER_1D, ctypes.c_int, ctypes.c_bool])

def smooth_field_using_KdTree(smoothing_kernel_radius_in_metres, kdtree_pointer,
                             lat, lon, area_size, f, number_of_threads, show_output_log=True):
    """Smooth one field; show_output_log=False suppresses native timing/progress."""
    lat, lon = _coordinates(lat, lon)
    area_size, f = _fields(area_size, f)
    if lat.shape != f.shape:
        raise ValueError("lat, lon, area_size and f must have the same shape")
    pointer = _require_handle(kdtree_pointer, "tree", lat.size)
    radius = _radius(smoothing_kernel_radius_in_metres)
    number_of_threads = _thread_count(number_of_threads)
    f_smoothed = np.zeros(f.shape, dtype=np.float64)
    libc.smooth_field_using_KdTree_ctypes(radius, pointer, lat, lon, area_size,
                                        f, f.size, f_smoothed, number_of_threads, bool(show_output_log))
    return f_smoothed

_bind("smooth_multiple_fields_simultaneously_using_KdTree_ctypes", [ctypes.c_double,
      ctypes.c_void_p, ND_POINTER_1D, ND_POINTER_1D, ND_POINTER_1D, ND_POINTER_2D,
      ctypes.c_size_t, ctypes.c_size_t, ND_POINTER_2D, ctypes.c_int, ctypes.c_bool])

def smooth_multiple_fields_simultaneously_using_KdTree(smoothing_kernel_radius_in_metres,
        kdtree_pointer, lat, lon, area_size, f, number_of_threads, show_output_log=True):
    """Smooth fields together using one shared 1D area_size vector.

    f has shape (number_of_fields, number_of_points). All fields use identical
    points and areas, including a common zero-area missing-value mask.
    show_output_log=False suppresses native timing/progress.
    """
    lat, lon = _coordinates(lat, lon)
    area_size = _float_array(area_size, "area_size")
    f = _float_array(f, "f", 2)
    if f.shape[1] != area_size.size:
        raise ValueError("Each field must have one value per area element")
    if f.shape[1] != lat.size:
        raise ValueError("The second field dimension must match the number of coordinates")
    pointer = _require_handle(kdtree_pointer, "tree", lat.size)
    radius = _radius(smoothing_kernel_radius_in_metres)
    number_of_threads = _thread_count(number_of_threads)
    f_smoothed = np.zeros(f.shape, dtype=np.float64)
    libc.smooth_multiple_fields_simultaneously_using_KdTree_ctypes(radius, pointer,
        lat, lon, area_size, f, lat.size, f.shape[0], f_smoothed, number_of_threads, bool(show_output_log))
    return f_smoothed

# KD-tree cache data: S2 spatial ordering, collapsed subtrees and 8-bit deltas
# Generation, loading and smoothing accept show_output_log=True for native output.


_bind("generate_kdtree_cache_data_in_memory_ctypes",
      [ND_POINTER_1D, ND_POINTER_1D, ctypes.c_size_t, ctypes.c_double, ctypes.c_int, ctypes.c_bool], True)

def generate_kdtree_cache_data_in_memory(lat, lon, smoothing_kernel_radius_in_metres, number_of_threads, show_output_log=True):
    """Build one radius's KD-tree cache in memory and log its retained storage size."""
    lat, lon = _coordinates(lat, lon)
    return _new_handle("kdtree_cache", lat.size, libc.generate_kdtree_cache_data_in_memory_ctypes,
                       lat, lon, lat.size, _radius(smoothing_kernel_radius_in_metres),
                       _thread_count(number_of_threads), bool(show_output_log))

def _write_cache_to_disk(cache, kind, filename, function, show_output_log):
    """Publish a complete cache file; preserve any existing file on write failure."""
    pointer = _require_handle(cache, kind)
    destination = os.path.abspath(os.fsdecode(_filename(filename)))
    folder = os.path.dirname(destination)
    # The caller creates the destination directory. Staging on the same filesystem
    # lets os.replace publish the complete file atomically, including on Windows.
    with tempfile.TemporaryDirectory(prefix=".smoothing-cache-", dir=folder) as staging:
        temporary_file = os.path.join(staging, "cache.tmp")
        function(pointer, _filename(temporary_file))
        os.replace(temporary_file, destination)
    if show_output_log:
        print("----- Cache saved: " + os.path.abspath(destination), flush=True)


_bind("write_kdtree_cache_data_to_disk_ctypes", [ctypes.c_void_p, ctypes.c_char_p])


def save_kdtree_cache(kdtree_cache_data, filename, show_output_log=True):
    """Save a cache to the exact filename; replace an existing file only on success."""
    _write_cache_to_disk(kdtree_cache_data, "kdtree_cache", filename,
                        libc.write_kdtree_cache_data_to_disk_ctypes, show_output_log)


_bind("read_kdtree_cache_data_from_binary_file_ctypes",
      [ctypes.c_char_p, ctypes.c_bool], True)

def load_kdtree_cache(filename, show_output_log=True):
    """Load the specified file; its stored radius and point count are validated in C++."""
    # The native cache owns its point count and validates it on each smoothing call.
    return _new_handle("kdtree_cache", None,
        libc.read_kdtree_cache_data_from_binary_file_ctypes,
        _filename(filename), bool(show_output_log))

_bind("smooth_field_using_kdtree_cache_data_ctypes",
      [ND_POINTER_1D, ND_POINTER_1D, ctypes.c_size_t, ctypes.c_void_p, ND_POINTER_1D, ctypes.c_int, ctypes.c_bool])

def smooth_field_using_kdtree_cache_data(
        area_size, f, kdtree_cache_data, number_of_threads, show_output_log=True):
    area_size, f = _fields(area_size, f)
    kdtree_cache_data_pointer = _require_handle(kdtree_cache_data, "kdtree_cache")
    number_of_threads = _thread_count(number_of_threads)
    # Native smoothing writes every element before returning successfully.
    f_smoothed = np.empty(f.shape, dtype=np.float64)
    libc.smooth_field_using_kdtree_cache_data_ctypes(
        area_size, f, f.size, kdtree_cache_data_pointer, f_smoothed, number_of_threads, bool(show_output_log))
    return f_smoothed

_bind("smooth_multiple_fields_simultaneously_using_kdtree_cache_data_ctypes",
      [ND_POINTER_1D, ND_POINTER_2D, ctypes.c_size_t, ctypes.c_size_t,
       ctypes.c_void_p, ND_POINTER_2D, ctypes.c_int, ctypes.c_bool])

def smooth_multiple_fields_simultaneously_using_kdtree_cache_data(
        area_size, fields, kdtree_cache_data, number_of_threads, show_output_log=True):
    """Smooth a batch of fields, sharing area calculations and native buffers.

    fields has shape (number_of_fields, number_of_points) in original point order.
    Every field must use the same points and the supplied 1D area_size, including
    its zero-area missing-value mask. Fields are processed one at a time; each uses OpenMP.
    The cache can come from disk or generate_kdtree_cache_data_in_memory().
    """
    area_size = _float_array(area_size, "area_size")
    fields = _float_array(fields, "fields", 2)
    if fields.shape[1] != area_size.size:
        raise ValueError("Each field must have one value per area element")
    kdtree_cache_data_pointer = _require_handle(kdtree_cache_data, "kdtree_cache", area_size.size)
    number_of_threads = _thread_count(number_of_threads)
    # Native processing writes all output elements before a successful return.
    smoothed_fields = np.empty(fields.shape, dtype=np.float64)
    libc.smooth_multiple_fields_simultaneously_using_kdtree_cache_data_ctypes(
        area_size, fields, area_size.size, fields.shape[0], kdtree_cache_data_pointer,
        smoothed_fields, number_of_threads, bool(show_output_log))
    return smoothed_fields

_bind("free_kdtree_cache_data_ctypes", [ctypes.c_void_p])

def free_kdtree_cache_data(kdtree_cache_data):
    kdtree_cache_data_pointer = _require_handle(kdtree_cache_data, "kdtree_cache")
    libc.free_kdtree_cache_data_ctypes(kdtree_cache_data_pointer)
    kdtree_cache_data._pointer = None

# Overlap cache data: S2 spatial ordering, point differences and 8-bit deltas
# Generation, loading and smoothing accept show_output_log=True for native output.
_bind("generate_overlap_cache_data_in_memory_ctypes",
      [ND_POINTER_1D, ND_POINTER_1D, ctypes.c_size_t, ctypes.c_double, ctypes.c_int, ctypes.c_bool], True)


def generate_overlap_cache_data_in_memory(lat, lon, smoothing_kernel_radius_in_metres,
                                         number_of_threads, show_output_log=True):
    """Build one radius's overlap cache in memory and log its retained storage size."""
    lat, lon = _coordinates(lat, lon)
    return _new_handle("overlap_cache", lat.size, libc.generate_overlap_cache_data_in_memory_ctypes,
                       lat, lon, lat.size, _radius(smoothing_kernel_radius_in_metres),
                       _thread_count(number_of_threads), bool(show_output_log))


_bind("write_overlap_cache_data_to_disk_ctypes", [ctypes.c_void_p, ctypes.c_char_p])


def save_overlap_cache(overlap_cache_data, filename, show_output_log=True):
    """Save a cache to the exact filename; replace an existing file only on success."""
    _write_cache_to_disk(overlap_cache_data, "overlap_cache", filename,
                        libc.write_overlap_cache_data_to_disk_ctypes, show_output_log)


_bind("read_overlap_cache_data_from_binary_file_ctypes",
      [ctypes.c_char_p, ctypes.c_bool], True)

def load_overlap_cache(filename, show_output_log=True):
    """Load the specified file; its stored radius and point count are validated in C++."""
    return _new_handle("overlap_cache", None,
        libc.read_overlap_cache_data_from_binary_file_ctypes, _filename(filename), bool(show_output_log))

_bind("free_overlap_cache_data_ctypes", [ctypes.c_void_p])

def free_overlap_cache_data(overlap_cache_data):
    overlap_cache_data_pointer = _require_handle(overlap_cache_data, "overlap_cache")
    libc.free_overlap_cache_data_ctypes(overlap_cache_data_pointer)
    overlap_cache_data._pointer = None

_bind("smooth_field_using_overlap_detection_ctypes",
      [ND_POINTER_1D, ND_POINTER_1D, ctypes.c_size_t, ctypes.c_void_p, ND_POINTER_1D, ctypes.c_int, ctypes.c_bool])

def smooth_field_using_overlap_detection(area_size, f, overlap_cache_data, number_of_threads, show_output_log=True):
    area_size, f = _fields(area_size, f)
    overlap_cache_data_pointer = _require_handle(overlap_cache_data, "overlap_cache", f.size)
    number_of_threads = _thread_count(number_of_threads)
    # Native smoothing writes every element before returning successfully.
    f_smoothed = np.empty(f.shape, dtype=np.float64)
    libc.smooth_field_using_overlap_detection_ctypes(
        area_size, f, f.size, overlap_cache_data_pointer, f_smoothed, number_of_threads, bool(show_output_log))
    return f_smoothed

_bind("smooth_multiple_fields_simultaneously_using_overlap_detection_ctypes",
      [ND_POINTER_1D, ND_POINTER_2D, ctypes.c_size_t, ctypes.c_size_t,
       ctypes.c_void_p, ND_POINTER_2D, ctypes.c_int, ctypes.c_bool])

def smooth_multiple_fields_simultaneously_using_overlap_detection(
        area_size, fields, overlap_cache_data, number_of_threads, show_output_log=True):
    """Smooth multiple fields in one call, sharing area denominators and native buffers.

    fields has shape (number_of_fields, number_of_points). All rows use the
    supplied area_size, including its common zero-area missing-value mask.
    Fields are processed one at a time; each uses OpenMP.
    """
    area_size = _float_array(area_size, "area_size")
    fields = _float_array(fields, "fields", 2)
    if fields.shape[1] != area_size.size:
        raise ValueError("Each field must have one value per area element")
    overlap_cache_data_pointer = _require_handle(overlap_cache_data, "overlap_cache", area_size.size)
    number_of_threads = _thread_count(number_of_threads)
    # The native routine writes every output element before returning successfully.
    smoothed_fields = np.empty(fields.shape, dtype=np.float64)
    libc.smooth_multiple_fields_simultaneously_using_overlap_detection_ctypes(
        area_size, fields, area_size.size, fields.shape[0], overlap_cache_data_pointer,
        smoothed_fields, number_of_threads, bool(show_output_log))
    return smoothed_fields




