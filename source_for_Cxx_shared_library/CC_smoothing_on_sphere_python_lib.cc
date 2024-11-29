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


extern "C" void free_mem_double_array(double* a)
	{
	delete[] a;
	}

extern "C"  void * construct_KdTree_ctypes(const double *lat, const double *lon, size_t size)
	{

	vector <kdtree::Point_str> kdtree_points = generate_vector_of_kdtree_points_from_lat_lon_points_provided_as_arrays(lat, lon, size);

	// create kdtree
	auto begin2 = std::chrono::high_resolution_clock::now();
	kdtree::KdTree *kdtree = new kdtree::KdTree;
	kdtree->buildKdTree_and_do_not_change_the_kdtree_points_vector(kdtree_points);
	cout << "----- kdtree construction " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin2).count() * 1e-9 << " s" <<  endl;

	return((void*)kdtree);
	}

extern "C"  void free_KdTree_memory_ctypes(void *kdtree_void_pointer)
	{
	kdtree::KdTree *kdtree = (kdtree::KdTree*) kdtree_void_pointer;
	kdtree->free_memory();
	}

extern "C"  void smooth_field_using_KdTree_ctypes(const double r_kernel_in_metres, const void *kdtree_void_pointer, const double *lat, const double *lon, const double *area_size, const double *f,  size_t number_of_points, double *f_smoothed)
	{
	kdtree::KdTree *kdtree = (kdtree::KdTree*) kdtree_void_pointer;

	smooth_field_using_kd_tree(r_kernel_in_metres, *kdtree, lat, lon, area_size, f, number_of_points, f_smoothed);
	}

extern "C"  void generate_smoothing_data_for_the_overlap_detection_and_write_it_to_disk_ctypes(const double *lat, const double *lon, size_t number_of_points, const double *smoothing_kernel_radius_in_metres, size_t smoothing_kernel_radius_in_metres_size, const char* output_folder_char, size_t starting_point)
	{

	vector <double> smoothing_kernel_radius_in_metres_vector;
	for (size_t il=0; il < smoothing_kernel_radius_in_metres_size; il++)
		smoothing_kernel_radius_in_metres_vector.push_back(smoothing_kernel_radius_in_metres[il]);

	//for (size_t il=0; il < smoothing_kernel_radius_in_metres_vector.size(); il++)
	//	cout << smoothing_kernel_radius_in_metres_vector[il] << endl;

	string output_folder = output_folder_char;

	generate_smoothing_data_for_the_overlap_detection_and_write_it_to_disk(lat, lon, number_of_points, smoothing_kernel_radius_in_metres_vector, output_folder, starting_point);
	}

extern "C"  void * read_smoothing_data_from_binary_file_ctypes(const char* output_folder_char, const double smoothing_kernel_radius_in_metres, const size_t number_of_points)
	{
	uint32_t number_of_points_uint32_t = number_of_points;

	string output_folder = output_folder_char;

	string fname = output_folder + "out_smoothing_information_r_"+output_leading_zero_string(round(smoothing_kernel_radius_in_metres),8)+"_m.bin";

	uint32_t **data_pointer = nullptr;
	Read_smoothing_data_from_binary_file(fname, data_pointer, number_of_points_uint32_t);

	return((void*) data_pointer);
	}

extern "C"  void free_smoothing_data_memory_ctypes(void *data_pointer_void, size_t number_of_points)
	{
	uint32_t** data_pointer = (uint32_t **)data_pointer_void;
	free_smoothing_data_memory(data_pointer, (uint32_t)number_of_points);
	data_pointer_void = nullptr;
	}

extern "C"  void smooth_field_using_overlap_detection_ctypes(const double *area_size, const double *f, size_t number_of_points, const void *data_pointer_void, double *f_smoothed)
	{
	smooth_field_using_overlap_detection(area_size, f, number_of_points, (uint32_t **)data_pointer_void,  f_smoothed);
	}


