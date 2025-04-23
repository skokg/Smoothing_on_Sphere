// Compile command for a linux system
// g++ -fopenmp -O2 -Wall -Wno-unused-result -Wno-unknown-pragmas -shared -o smoothing_on_sphere_Cxx_shared_library.so -fPIC CC_smoothing_on_sphere_python_lib.cc

// the NUMBER_OF_THREADS specifies how many cores should be utilized via OpenMP parallel computation
#define NUMBER_OF_THREADS 10

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <sstream>
#include <sys/stat.h>
#include <sys/resource.h>
#include <string.h>
#include <random>
#include <chrono>
#include <cstring>

using namespace std;

#define BAD_DATA_FLOAT -9999

// da prav dela error(..) - da prav displaya line number
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__)
#define FU __PRETTY_FUNCTION__
#define TT "\t"
#define tabt "\t"
#define ERRORIF(x) if (x) error(AT,FU, #x)

const rlim_t kStackSize = 1000 * 1024 * 1024;   // min stack size = 16 MB

long random_number_seed=-1;

#include "CU_utils_subset.cc"
#include "CU_smoothing_on_sphere_code.cc"


extern "C" void free_mem_double_array(double * const a)
	{
	delete[] a;
	}

extern "C"  void * construct_KdTree_ctypes(const double * const lat, const double * const lon, const size_t size)
	{

	vector <kdtree::Point_str> kdtree_points = generate_vector_of_kdtree_points_from_lat_lon_points_provided_as_arrays(lat, lon, size);

	// create kdtree
	auto begin2 = std::chrono::high_resolution_clock::now();
	kdtree::KdTree *kdtree = new kdtree::KdTree;
	kdtree->buildKdTree(kdtree_points);
	cout << "----- kdtree construction " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin2).count() * 1e-9 << " s" <<  endl;

	return((void*)kdtree);
	}

extern "C"  void free_KdTree_memory_ctypes(void * const kdtree_void_pointer)
	{
	kdtree::KdTree *kdtree = (kdtree::KdTree*) kdtree_void_pointer;
	kdtree->free_memory();
	}

extern "C"  void smooth_field_using_KdTree_ctypes(const double r_kernel_in_metres, const void * const kdtree_void_pointer, const double * const lat, const double * const lon, const double * const area_size, const double * const f,  const size_t number_of_points, double * const f_smoothed)
	{
	const kdtree::KdTree * const kdtree = (const kdtree::KdTree * const) kdtree_void_pointer;

	smooth_field_using_kd_tree(r_kernel_in_metres, *kdtree, lat, lon, area_size, f, number_of_points, f_smoothed);
	}

extern "C"  void smooth_multiple_fields_simultaneously_using_KdTree_ctypes(const double r_kernel_in_metres, const void * const kdtree_void_pointer, const double * const lat, const double * const lon, const double * const area_size_multiple_fields_2D_numpy_array, const double * const f_multiple_fields_2D_numpy_array,  const size_t number_of_points, const size_t number_of_fields, double * const f_smoothed_multiple_fields_2D_numpy_array)
	{
	const kdtree::KdTree * const kdtree = (const kdtree::KdTree * const) kdtree_void_pointer;

	vector <double*> area_size_pointers;
	for (size_t iff=0; iff < number_of_fields; iff++)
		area_size_pointers.push_back((double *)area_size_multiple_fields_2D_numpy_array + iff*number_of_points);

	vector <double*> f_pointers;
	for (size_t iff=0; iff < number_of_fields; iff++)
		f_pointers.push_back((double *)f_multiple_fields_2D_numpy_array + iff*number_of_points);

	vector <double*> f_smoothed_pointers;
	for (size_t iff=0; iff < number_of_fields; iff++)
		f_smoothed_pointers.push_back((double *)f_smoothed_multiple_fields_2D_numpy_array + iff*number_of_points);

	smooth_field_using_kd_tree_multiple_fields_simultaneously(r_kernel_in_metres, *kdtree, lat, lon, &area_size_pointers[0], &f_pointers[0], number_of_points, number_of_fields, &f_smoothed_pointers[0]);
	}

extern "C"  void generate_smoothing_data_for_the_overlap_detection_and_write_it_to_disk_ctypes(const double * const lat, const double * const lon, const size_t number_of_points, const double * const smoothing_kernel_radius_in_metres, const size_t smoothing_kernel_radius_in_metres_size, const char * const output_folder_char, const size_t starting_point)
	{

	vector <double> smoothing_kernel_radius_in_metres_vector;
	for (size_t il=0; il < smoothing_kernel_radius_in_metres_size; il++)
		smoothing_kernel_radius_in_metres_vector.push_back(smoothing_kernel_radius_in_metres[il]);

	//for (size_t il=0; il < smoothing_kernel_radius_in_metres_vector.size(); il++)
	//	cout << smoothing_kernel_radius_in_metres_vector[il] << endl;

	string output_folder = output_folder_char;

	generate_smoothing_data_for_the_overlap_detection_and_write_it_to_disk(lat, lon, number_of_points, smoothing_kernel_radius_in_metres_vector, output_folder, starting_point);
	}

extern "C"  void * read_smoothing_data_from_binary_file_ctypes(const char * const output_folder_char, const double smoothing_kernel_radius_in_metres, const size_t number_of_points)
	{
	uint32_t number_of_points_uint32_t = number_of_points;

	string output_folder = output_folder_char;

	string fname = output_folder + "out_smoothing_information_r_"+output_leading_zero_string(round(smoothing_kernel_radius_in_metres),8)+"_m.bin";

	uint32_t **data_pointer = nullptr;
	Read_smoothing_data_from_binary_file(fname, data_pointer, number_of_points_uint32_t);

	return((void*) data_pointer);
	}

extern "C"  void free_smoothing_data_memory_ctypes(void * data_pointer_void, const size_t number_of_points)
	{
	uint32_t** data_pointer = (uint32_t **)data_pointer_void;
	free_smoothing_data_memory(data_pointer, (uint32_t)number_of_points);
	data_pointer_void = nullptr;
	}

extern "C"  void smooth_field_using_overlap_detection_ctypes(const double * const area_size, const double * const f, const size_t number_of_points, const void * const data_pointer_void, double * const f_smoothed)
	{
	smooth_field_using_overlap_detection(area_size, f, number_of_points, (uint32_t **)data_pointer_void,  f_smoothed);
	}


