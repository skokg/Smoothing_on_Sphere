#include <iostream>
#include <algorithm>
#include <vector>
#include <sstream>
#include <chrono>
#include <cstring>
#include <cmath>
#include <cstdint>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <limits>
#include <memory>

using namespace std;

#include "CU_smoothing_on_sphere_code.cc"

// C++ errors are caught at the C boundary. Error reporting itself must not
// allocate or throw, including when handling std::bad_alloc.
namespace
{
thread_local char smoothing_last_error[2048] = {};

void smoothing_clear_error() noexcept
{
    smoothing_last_error[0] = '\0';
}

void smoothing_record_error(const char *function, const char *message) noexcept
{
    const char *parts[] = {function, ": ", message};
    size_t position = 0;
    for (const char *part : parts)
    {
        if (part == nullptr) continue;
        while (*part != '\0' && position < sizeof(smoothing_last_error) - 1)
            smoothing_last_error[position++] = *part++;
    }
    smoothing_last_error[position] = '\0';
}
}

// Borrowed string, valid until the next operation on the same native thread.
// Reading the message does not clear it.
extern "C" const char *smoothing_get_last_error() noexcept
{
    return smoothing_last_error;
}

#define SMOOTHING_C_API_CATCH(failure_result) \
    catch (const std::exception &exception) { smoothing_record_error(__func__, exception.what()); return failure_result; } \
    catch (...) { smoothing_record_error(__func__, "Unknown C++ exception in smoothing library."); return failure_result; }

extern "C" void *construct_KdTree_ctypes(const double * const lat, const double * const lon, const size_t size, const int number_of_threads) noexcept
try
{
    smoothing_clear_error();
    smoothing_kdtree::validate_number_of_threads(number_of_threads);
    const auto begin = std::chrono::high_resolution_clock::now();
    auto kdtree_points = generate_vector_of_kdtree_points_from_lat_lon_points_provided_as_arrays(lat, lon, size, number_of_threads);
    unique_ptr<kdtree::KdTree<3>> tree(new kdtree::KdTree<3>());
    if (!kdtree_points.empty() && !tree->buildKdTree(kdtree_points))
        throw runtime_error("Cannot construct smoothing KD-tree.");
    cout << "----- KD-tree construction " << std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - begin).count() << " s" << endl;
    return tree.release();
}
SMOOTHING_C_API_CATCH(nullptr)

extern "C" int free_KdTree_memory_ctypes(void * const kdtree_void_pointer) noexcept
try
{
    smoothing_clear_error();
    delete static_cast<kdtree::KdTree<3> *>(kdtree_void_pointer);
    return 0;
}
SMOOTHING_C_API_CATCH(-1)

extern "C" int smooth_field_using_KdTree_ctypes(const double r_kernel_in_metres, const void * const kdtree_void_pointer, const double * const lat, const double * const lon, const double * const area_size, const double * const f,  const size_t number_of_points, double * const f_smoothed, const int number_of_threads, const bool show_output_log = true) noexcept
try
{
    smoothing_clear_error();
    smoothing_kdtree::validate_number_of_threads(number_of_threads);
    smoothing_kdtree::require_buffer(kdtree_void_pointer, 1);
    const auto &tree = *static_cast<const kdtree::KdTree<3> *>(kdtree_void_pointer);
    smooth_field_using_kd_tree(r_kernel_in_metres, tree, lat, lon, area_size, f, number_of_points, f_smoothed, number_of_threads, show_output_log);
    return 0;
}
SMOOTHING_C_API_CATCH(-1)

extern "C" int smooth_multiple_fields_simultaneously_using_KdTree_ctypes(const double r_kernel_in_metres, const void * const kdtree_void_pointer, const double * const lat, const double * const lon, const double * const area_size, const double * const f_multiple_fields_2D_numpy_array,  const size_t number_of_points, const size_t number_of_fields, double * const f_smoothed_multiple_fields_2D_numpy_array, const int number_of_threads, const bool show_output_log = true) noexcept
try
{
    smoothing_clear_error();
    smoothing_kdtree::validate_number_of_threads(number_of_threads);
    smoothing_kdtree::require_buffer(kdtree_void_pointer, 1);
    const auto &tree = *static_cast<const kdtree::KdTree<3> *>(kdtree_void_pointer);
    smoothing_kdtree::validate_point_count(number_of_points);
    validate_smoothing_radius(r_kernel_in_metres);
    if (tree.count_number_of_nonnull_nodes() != number_of_points) throw invalid_argument("Point count differs from smoothing KD-tree.");
    if (number_of_points == 0 || number_of_fields == 0) return 0;
    if (number_of_fields > static_cast<size_t>(numeric_limits<ptrdiff_t>::max()) / sizeof(double) / number_of_points)
        throw length_error("Multiple-field array dimensions are too large.");
    smoothing_kdtree::require_buffer(area_size, number_of_points);
    smoothing_kdtree::require_buffer(f_multiple_fields_2D_numpy_array, number_of_points);
    smoothing_kdtree::require_buffer(f_smoothed_multiple_fields_2D_numpy_array, number_of_points);
    vector<const double *> fields(number_of_fields);
    vector<double *> outputs(number_of_fields);
    for (size_t i=0; i<number_of_fields; ++i)
    {
        fields[i] = f_multiple_fields_2D_numpy_array + i * number_of_points;
        outputs[i] = f_smoothed_multiple_fields_2D_numpy_array + i * number_of_points;
    }
    smooth_field_using_kd_tree_multiple_fields_simultaneously(r_kernel_in_metres, tree, lat, lon,
        area_size, fields.data(), number_of_points, number_of_fields, outputs.data(), number_of_threads, show_output_log);
    return 0;
}
SMOOTHING_C_API_CATCH(-1)


extern "C" int generate_overlap_cache_data_and_write_to_disk_ctypes(const double * const lat, const double * const lon, const size_t number_of_points, const double * const smoothing_kernel_radius_in_metres, const size_t smoothing_kernel_radius_in_metres_size, const char * const overlap_cache_data_folder, const size_t starting_point, const int number_of_threads, const bool show_output_log = true) noexcept
try
{
    smoothing_clear_error();
    smoothing_kdtree::validate_number_of_threads(number_of_threads);
    smoothing_kdtree::require_buffer(overlap_cache_data_folder, 1);
    smoothing_kdtree::require_buffer(smoothing_kernel_radius_in_metres, smoothing_kernel_radius_in_metres_size);
    if (smoothing_kernel_radius_in_metres_size > static_cast<size_t>(numeric_limits<ptrdiff_t>::max()) / sizeof(double))
        throw length_error("Smoothing radius array is too large.");
    vector<double> radii;
    if (smoothing_kernel_radius_in_metres_size != 0)
        radii.assign(smoothing_kernel_radius_in_metres,
                     smoothing_kernel_radius_in_metres + smoothing_kernel_radius_in_metres_size);
    generate_overlap_cache_data_and_write_to_disk(
        lat, lon, number_of_points, radii, string(overlap_cache_data_folder), starting_point, number_of_threads, show_output_log);
    return 0;
}
SMOOTHING_C_API_CATCH(-1)

extern "C" void *read_overlap_cache_data_from_binary_file_ctypes(const char * const overlap_cache_data_folder, const double smoothing_kernel_radius_in_metres, const bool show_output_log = true) noexcept
try
{
    smoothing_clear_error();
    smoothing_kdtree::require_buffer(overlap_cache_data_folder, 1);
    const double radius = cap_smoothing_radius(smoothing_kernel_radius_in_metres);
    const string cache_filename = make_overlap_cache_filename(string(overlap_cache_data_folder), radius);
    OverlapCacheData *overlap_cache_data_pointer = nullptr;
    read_overlap_cache_data_from_binary_file(cache_filename, overlap_cache_data_pointer, show_output_log);
    unique_ptr<OverlapCacheData> overlap_cache_data(overlap_cache_data_pointer);
    if (overlap_cache_data->radius != radius)
        throw runtime_error("Overlap cache radius mismatch.");
    return overlap_cache_data.release();
}
SMOOTHING_C_API_CATCH(nullptr)

extern "C" int free_overlap_cache_data_ctypes(void * overlap_cache_data_void_pointer) noexcept
try
{
    smoothing_clear_error();
    auto overlap_cache_data_pointer = static_cast<OverlapCacheData *>(overlap_cache_data_void_pointer);
    free_overlap_cache_data(overlap_cache_data_pointer);
    return 0;
}
SMOOTHING_C_API_CATCH(-1)

extern "C" int smooth_field_using_overlap_detection_ctypes(const double * const area_size, const double * const f, const size_t number_of_points, const void * const overlap_cache_data_void_pointer, double * const f_smoothed, const int number_of_threads, const bool show_output_log = true) noexcept
try
{
    smoothing_clear_error();
    smoothing_kdtree::validate_number_of_threads(number_of_threads);
    smooth_field_using_overlap_detection(area_size, f, number_of_points,
        static_cast<const OverlapCacheData *>(overlap_cache_data_void_pointer), f_smoothed, number_of_threads, show_output_log);
    return 0;
}
SMOOTHING_C_API_CATCH(-1)

extern "C" int smooth_multiple_fields_simultaneously_using_overlap_detection_ctypes(
    const double * const area_size, const double * const fields,
    const size_t number_of_points, const size_t number_of_fields,
    const void * const overlap_cache_data_void_pointer, double * const smoothed_fields,
    const int number_of_threads, const bool show_output_log = true) noexcept
try
{
    smoothing_clear_error();
    smooth_multiple_fields_simultaneously_using_overlap_detection(
        area_size, fields, number_of_points, number_of_fields,
        static_cast<const OverlapCacheData *>(overlap_cache_data_void_pointer),
        smoothed_fields, number_of_threads, show_output_log);
    return 0;
}
SMOOTHING_C_API_CATCH(-1)

extern "C" int generate_kdtree_cache_data_and_write_to_disk_ctypes(const double * const lat, const double * const lon, const size_t number_of_points, const double smoothing_kernel_radius_in_metres, const char * const kdtree_cache_data_folder, const int number_of_threads, const bool show_output_log = true) noexcept
try
{
    smoothing_clear_error();
    smoothing_kdtree::validate_number_of_threads(number_of_threads);
    smoothing_kdtree::require_buffer(kdtree_cache_data_folder, 1);
    generate_kdtree_cache_data_and_write_to_disk(
        lat, lon, number_of_points, smoothing_kernel_radius_in_metres, string(kdtree_cache_data_folder), number_of_threads, show_output_log);
    return 0;
}
SMOOTHING_C_API_CATCH(-1)

extern "C" int generate_kdtree_cache_data_for_multiple_radii_ctypes(
    const double *lat,const double *lon,const size_t n,const double *radii,
    const size_t radius_count,const char *kdtree_cache_data_folder,const int threads, const bool show_output_log = true) noexcept
try
{
    smoothing_clear_error();
    smoothing_kdtree::require_buffer(kdtree_cache_data_folder,1);
    if (radius_count==0) throw invalid_argument("At least one smoothing radius is required.");
    smoothing_kdtree::require_buffer(radii,radius_count);
    saved_s2::generate_cache_data_and_write_batched(
        lat,lon,n,vector<double>(radii,radii+radius_count),string(kdtree_cache_data_folder),threads, show_output_log);
    return 0;
}
SMOOTHING_C_API_CATCH(-1)

extern "C" void *read_kdtree_cache_data_from_binary_file_ctypes(const char * const kdtree_cache_data_folder, const double smoothing_kernel_radius_in_metres, const bool show_output_log = true) noexcept
try
{
    smoothing_clear_error();
    smoothing_kdtree::require_buffer(kdtree_cache_data_folder, 1);
    const double radius = cap_smoothing_radius(smoothing_kernel_radius_in_metres);
    auto kdtree_cache_data = saved_s2::read_cache_data_from_binary_file(saved_s2::make_cache_filename(string(kdtree_cache_data_folder), radius), show_output_log);
    if (kdtree_cache_data->radius != radius)
        throw runtime_error("KD-tree cache radius mismatch.");
    return kdtree_cache_data.release();
}
SMOOTHING_C_API_CATCH(nullptr)

extern "C" int smooth_field_using_kdtree_cache_data_ctypes(const double * const area_size, const double * const f, const size_t number_of_points, const char * const kdtree_cache_data_pointer, double * const f_smoothed, const int number_of_threads, const bool show_output_log = true) noexcept
try
{
    smoothing_clear_error();
    smoothing_kdtree::validate_number_of_threads(number_of_threads);
    smooth_field_using_kdtree_cache_data(
        area_size, f, number_of_points, kdtree_cache_data_pointer, f_smoothed, number_of_threads, show_output_log);
    return 0;
}
SMOOTHING_C_API_CATCH(-1)

extern "C" int smooth_multiple_fields_simultaneously_using_kdtree_cache_data_ctypes(
    const double *area_size,const double *fields,const size_t number_of_points,
    const size_t number_of_fields,const void *kdtree_cache_data_void_pointer,
    double *smoothed_fields,const int number_of_threads, const bool show_output_log = true) noexcept
try
{
    smoothing_clear_error();
    smooth_multiple_fields_simultaneously_using_kdtree_cache_data(
        area_size,fields,number_of_points,number_of_fields,
        static_cast<const KdTreeCacheData *>(kdtree_cache_data_void_pointer),
        smoothed_fields,number_of_threads, show_output_log);
    return 0;
}
SMOOTHING_C_API_CATCH(-1)

extern "C" int free_kdtree_cache_data_ctypes(void * kdtree_cache_data_void_pointer) noexcept
try
{
    smoothing_clear_error();
    auto kdtree_cache_data_pointer = static_cast<char *>(kdtree_cache_data_void_pointer);
    free_kdtree_cache_data(kdtree_cache_data_pointer);
    return 0;
}
SMOOTHING_C_API_CATCH(-1)

extern "C" void *generate_kdtree_cache_data_in_memory_ctypes(const double *lat, const double *lon,
    size_t count, double radius, int number_of_threads, const bool show_output_log = true) noexcept
try {
    smoothing_clear_error();
    return saved_s2::generate_cache_data_in_memory(lat, lon, count, radius, number_of_threads, show_output_log).release();
}
SMOOTHING_C_API_CATCH(nullptr)

#undef SMOOTHING_C_API_CATCH




