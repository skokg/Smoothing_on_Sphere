#include "CU_kdtree_with_index.cc"
#include <exception>
#include <utility>
#include <locale>
#include <iomanip>
#include <cerrno>
#include <cstdlib>
#ifdef _WIN32
#include <direct.h>
#else
#include <unistd.h>
#endif

// Shared KD-tree types and input validation.
namespace smoothing_kdtree
{
using Point = kdtree::Point_str<3>;
using Node = kdtree::KdTreeNode<3>;
using Index = kdtree::IndexType;
using Coordinates = std::array<kdtree::PointType, 3>;

// Per-call OpenMP team size; no global runtime setting is changed.
inline void validate_number_of_threads(const int number_of_threads)
{
    if (number_of_threads <= 0)
        throw std::invalid_argument("Number of threads must be positive.");
}

inline void validate_point_count(const size_t count)
{
    if (count > static_cast<size_t>(std::numeric_limits<int32_t>::max()))
        throw std::length_error("Smoothing supports at most INT32_MAX points.");
}

inline void require_buffer(const void *pointer, const size_t count)
{
    if (count != 0 && pointer == nullptr)
        throw std::invalid_argument("A nonempty smoothing array has a null pointer.");
}

inline void validate_lat_lon(const double *lat, const double *lon, const size_t count)
{
    validate_point_count(count);
    require_buffer(lat, count);
    require_buffer(lon, count);
    // Check before entering an OpenMP region, where exceptions cannot escape.
    for (size_t i = 0; i < count; ++i)
        if (!std::isfinite(lat[i]) || !std::isfinite(lon[i]))
            throw std::invalid_argument("Latitude and longitude must be finite.");
}
}

// Utilities used by smoothing and its C interface.

// Release the vector's storage as well as its elements.
template<typename T>
void free_vector_memory(std::vector<T>& vec)
	{
	std::vector<T>().swap(vec);
	}

double deg2rad(double kot)
	{
	return (kot*M_PI/180.0);
	}

void spherical_to_cartesian_coordinates(double lat, double lon, double r, double &x, double &y, double &z)
	{
	x=r*cos(lon) * cos(lat);
	y=r*sin(lon) * cos(lat);
	z=r*sin(lat);
	}

// Pad an integer with leading zeroes to the requested width.
string output_leading_zero_string(long il,int mest)
	{
	int ix;

	ostringstream s1;
	s1.str("");
	int stevke;
	string s;

	if (il==0) stevke=1;
	else stevke=(int)floor(log10((double)il))+1;

	s1 << il;
	s=s1.str();
	s1.str("");
	for (ix=0;ix< mest-stevke;ix++)
		{
		s1 << "0" << s;
		s=s1.str();
		s1.str("");
		}
	return(s);
	}

double round_to_digits(double x, int digits)
	{
	return(round(pow(10,(double)digits)*x)/pow(10,(double)digits));
	}

double Earth_radius= 6371*1000;

void validate_smoothing_radius(const double radius)
	{
	if (!isfinite(radius) || radius < 0)
		throw invalid_argument("Smoothing radius must be finite and nonnegative.");
	}

double maximum_smoothing_radius()
	{
	return acos(-1.0) * Earth_radius;
	}

double cap_smoothing_radius(const double radius)
	{
	validate_smoothing_radius(radius);
	const double maximum_radius = maximum_smoothing_radius();
	if (radius <= maximum_radius) return radius;
	ostringstream warning;
	warning.imbue(std::locale::classic());
	warning.precision(numeric_limits<double>::max_digits10);
	warning << "WARNING: Smoothing radius " << radius << " m exceeds half the Earth's circumference; capped to "
		<< maximum_radius << " m.";
	cerr << warning.str() << endl;
	return maximum_radius;
	}

// Stored records cannot be repaired by changing their header radius.
void validate_cache_radius(const double radius)
	{
	validate_smoothing_radius(radius);
	if (radius > maximum_smoothing_radius())
		throw runtime_error("Cache radius exceeds half the Earth's circumference; regenerate the cache.");
	}

string make_overlap_cache_filename(const string &folder, const double requested_radius)
	{
	const double radius = cap_smoothing_radius(requested_radius);
	return folder + "overlap_cache_data_r_" +
		output_leading_zero_string(static_cast<long>(round(radius)), 8) + "_m.bin";
	}

namespace smoothing_io
	{
	using File = unique_ptr<FILE, decltype(&fclose)>;

	string absolute_path(const string &name)
		{
#ifdef _WIN32
		unique_ptr<char, decltype(&std::free)> resolved(_fullpath(nullptr, name.c_str(), 0), &std::free);
		if (!resolved) throw runtime_error("Cannot resolve cache file path: " + name);
		return string(resolved.get());
#else
		if (!name.empty() && name[0] == '/') return name;
		vector<char> directory(256);
		while (!getcwd(directory.data(), directory.size()))
			{
			if (errno != ERANGE) throw runtime_error("Cannot resolve cache file path: " + name);
			if (directory.size() > directory.max_size() / 2)
				throw length_error("Working directory path is too long.");
			directory.resize(directory.size() * 2);
			}
		return string(directory.data()) + "/" + name;
#endif
		}

	File open(const string &name, const char *mode)
		{
		File file(fopen(name.c_str(), mode), &fclose);
		if (!file) throw runtime_error("Cannot open cache data file: " + name);
		return file;
		}

	File open_cache(const string &name, const char *mode, const bool show_output_log)
		{
		const string full_path = absolute_path(name);
		if (show_output_log)
			cout << (strcmp(mode, "rb") == 0 ? "Reading file: " : "Opening file for write: ")
				<< full_path << endl;
		return open(full_path, mode);
		}

	void read(FILE *file, void *data, const size_t element_size, const size_t count)
		{
		if (count != 0 && fread(data, element_size, count, file) != count)
			throw runtime_error("Cannot read complete cache data.");
		}

	void write(FILE *file, const void *data, const size_t element_size, const size_t count)
		{
		if (count != 0 && fwrite(data, element_size, count, file) != count)
			throw runtime_error("Cannot write complete cache data.");
		}

	void require_end(FILE *file)
		{
		if (fgetc(file) != EOF || ferror(file))
			throw runtime_error("Unexpected trailing data or read failure in smoothing cache.");
		}

	void close(File &file)
		{
		if (fclose(file.release()) != 0)
			throw runtime_error("Cannot close cache data file.");
		}
	}

double great_circle_distance_to_euclidian_distance(const double great_circle_distance)
	{
	const double radius = cap_smoothing_radius(great_circle_distance);
	// At the maximum radius every point is inside. An infinite search bound
	// also includes antipodes whose float coordinates round beyond the diameter.
	if (radius == maximum_smoothing_radius()) return numeric_limits<float>::infinity();
	return(2*Earth_radius*sin(radius/(2*Earth_radius)));
	}

// Geometry helpers for KD-tree radius traversal.
namespace smoothing_kdtree
{
inline float distance_squared(const Point &a, const Point &b) noexcept
{
    float sum = 0;
    for (size_t axis = 0; axis < 3; ++axis)
    {
        const float delta = b.coords[axis] - a.coords[axis];
        sum += delta * delta;
    }
    return sum;
}

inline bool bounds_intersect(const Coordinates &center, const float radius_squared,
                             const Coordinates &minimum, const Coordinates &maximum) noexcept
{
    float sum = 0;
    for (size_t axis = 0; axis < 3; ++axis)
    {
        const float nearest = std::max(minimum[axis], std::min(center[axis], maximum[axis]));
        const float delta = center[axis] - nearest;
        sum += delta * delta;
    }
    return sum <= radius_squared;
}

inline bool bounds_inside(const Coordinates &center, const float radius_squared,
                          const Coordinates &minimum, const Coordinates &maximum) noexcept
{
    float sum = 0;
    for (size_t axis = 0; axis < 3; ++axis)
    {
        const float delta = std::max(std::fabs(center[axis] - minimum[axis]),
                                     std::fabs(center[axis] - maximum[axis]));
        sum += delta * delta;
    }
    return sum <= radius_squared;
}
}

void generate_fxareasize_and_areasize_Bounding_Box_data(const smoothing_kdtree::Node * const node, const double * const f, const double * const area_size,  vector <double> &fxareasize_BB_data, vector <double> &area_size_BB_data)
	{
	if (node == nullptr) return;
	size_t index = node->val.index;
	double fxareasize_BB = f[index]*area_size[index];
	double areasize_BB = area_size[index];
	if(node->left_child() != nullptr)
		{
		generate_fxareasize_and_areasize_Bounding_Box_data(node->left_child(), f, area_size, fxareasize_BB_data, area_size_BB_data);
		fxareasize_BB+= fxareasize_BB_data[node->left_child()->val.index];
		areasize_BB+= area_size_BB_data[node->left_child()->val.index];
		}
	if(node->right_child() != nullptr)
		{
		generate_fxareasize_and_areasize_Bounding_Box_data(node->right_child(), f, area_size, fxareasize_BB_data, area_size_BB_data);
		fxareasize_BB+= fxareasize_BB_data[node->right_child()->val.index];
		areasize_BB+= area_size_BB_data[node->right_child()->val.index];
		}

	fxareasize_BB_data[index] =fxareasize_BB;
	area_size_BB_data[index] =areasize_BB;
	}

void generate_fxareasize_and_areasize_Bounding_Box_data_multiple_fields_simultaneously(const smoothing_kdtree::Node * const node, const double * const * const f, const double * const area_size,  const size_t number_of_fields, vector <vector <double>> &fxareasize_BB_data, vector <double> &area_size_BB_data  )
	{
	if (node == nullptr) return;
	size_t index = node->val.index;
	vector <double> fxareasize_BB (number_of_fields,0);
	double areasize_BB = area_size[index];
	for (unsigned long iff=0; iff < number_of_fields; iff++)
		{
		fxareasize_BB[iff] = f[iff][index]*area_size[index];
		}
	if(node->left_child() != nullptr)
		{
		generate_fxareasize_and_areasize_Bounding_Box_data_multiple_fields_simultaneously(node->left_child(), f, area_size, number_of_fields, fxareasize_BB_data, area_size_BB_data);
		areasize_BB += area_size_BB_data[node->left_child()->val.index];
		for (unsigned long iff=0; iff < number_of_fields; iff++)
			{
			fxareasize_BB[iff] += fxareasize_BB_data[iff][node->left_child()->val.index];
			}
		}
	if(node->right_child() != nullptr)
		{
		generate_fxareasize_and_areasize_Bounding_Box_data_multiple_fields_simultaneously(node->right_child(), f, area_size, number_of_fields, fxareasize_BB_data, area_size_BB_data);
		areasize_BB += area_size_BB_data[node->right_child()->val.index];
		for (unsigned long iff=0; iff < number_of_fields; iff++)
			{
			fxareasize_BB[iff] += fxareasize_BB_data[iff][node->right_child()->val.index];
			}
		}

	for (unsigned long iff=0; iff < number_of_fields; iff++)
		{
		fxareasize_BB_data[iff][index] =fxareasize_BB[iff];
		}
	area_size_BB_data[index] = areasize_BB;
	}

void get_fxareasize_and_areasize_sums_of_points_in_radius( const smoothing_kdtree::Node * const node, const kdtree::KdTree<3> &kdtree, const smoothing_kdtree::Point& point, const float distance, const float distance_sqaured, const vector <double> &fxareasize, const double * const areasize, const vector <double> &fxareasize_BB_data, const vector <double> &area_size_BB_data, double &fxareasize_sum, double &areasize_sum)
	{

	if(node == nullptr) return;

	if (smoothing_kdtree::bounds_intersect( point.coords, distance_sqaured, node->coords_min, node->coords_max ))
		{
		float dist_sqr = smoothing_kdtree::distance_squared(point, node->val);
		bool inside_sphere = false;
		if(dist_sqr <= distance_sqaured)
			inside_sphere = true;

		if (inside_sphere && smoothing_kdtree::bounds_inside( point.coords, distance_sqaured, node->coords_min, node->coords_max))
			{
			fxareasize_sum+=fxareasize_BB_data[node->val.index];
			areasize_sum+=area_size_BB_data[node->val.index];
			}

		else
			{
			if (inside_sphere)
				{
				fxareasize_sum+=fxareasize[node->val.index];
				areasize_sum+=areasize[node->val.index];
				}

			get_fxareasize_and_areasize_sums_of_points_in_radius(node->right_child(), kdtree, point, distance, distance_sqaured, fxareasize, areasize, fxareasize_BB_data, area_size_BB_data,fxareasize_sum, areasize_sum );
			get_fxareasize_and_areasize_sums_of_points_in_radius(node->left_child(), kdtree, point, distance, distance_sqaured, fxareasize, areasize, fxareasize_BB_data, area_size_BB_data,fxareasize_sum, areasize_sum );
			}
		}
	}

void get_fxareasize_and_areasize_sums_of_points_in_radius_multiple_fields_simultaneously( const smoothing_kdtree::Node * const node, const kdtree::KdTree<3> &kdtree, const smoothing_kdtree::Point& point, const float distance, const float distance_sqaured, const vector <vector <double>> &fxareasize, const double * const areasize, const vector <vector <double>> &fxareasize_BB_data, const vector <double> &area_size_BB_data, vector <double> &fxareasize_sum, double &areasize_sum, const size_t number_of_fields)
	{

	if(node == nullptr) return;

	if (smoothing_kdtree::bounds_intersect( point.coords, distance_sqaured, node->coords_min, node->coords_max ))
		{
		float dist_sqr = smoothing_kdtree::distance_squared(point, node->val);
		bool inside_sphere = false;
		if(dist_sqr <= distance_sqaured)
			inside_sphere = true;

		if (inside_sphere && smoothing_kdtree::bounds_inside( point.coords, distance_sqaured, node->coords_min, node->coords_max))
			{
			for (unsigned long iff=0; iff < number_of_fields; iff++)
				{
				fxareasize_sum[iff]+=fxareasize_BB_data[iff][node->val.index];
				}
			areasize_sum+=area_size_BB_data[node->val.index];
			}

		else
			{
			if (inside_sphere)
				{
				for (unsigned long iff=0; iff < number_of_fields; iff++)
					{
					fxareasize_sum[iff]+=fxareasize[iff][node->val.index];
					}
				areasize_sum+=areasize[node->val.index];
				}

			get_fxareasize_and_areasize_sums_of_points_in_radius_multiple_fields_simultaneously(node->right_child(), kdtree, point, distance, distance_sqaured, fxareasize, areasize, fxareasize_BB_data, area_size_BB_data,fxareasize_sum, areasize_sum, number_of_fields );
			get_fxareasize_and_areasize_sums_of_points_in_radius_multiple_fields_simultaneously(node->left_child(), kdtree, point, distance, distance_sqaured, fxareasize, areasize, fxareasize_BB_data, area_size_BB_data,fxareasize_sum, areasize_sum, number_of_fields );
			}
		}
	}

vector <smoothing_kdtree::Point> generate_vector_of_kdtree_points_from_lat_lon_points_provided_as_arrays(const double * const lat, const double * const lon, const size_t size, const int number_of_threads)
	{
    smoothing_kdtree::validate_number_of_threads(number_of_threads);
	smoothing_kdtree::validate_lat_lon(lat, lon, size);
	vector <smoothing_kdtree::Point> kdtree_points(size);
	#pragma omp parallel for num_threads(number_of_threads)
	for (size_t il=0; il < size; il++)
		{
		double x,y,z;
		spherical_to_cartesian_coordinates(deg2rad(lat[il]), deg2rad(lon[il]), Earth_radius, x, y, z);
		kdtree_points[il].coords = {{(float)x, (float)y, (float)z}};
		// The count was checked before entering the parallel region.
		kdtree_points[il].index = static_cast<smoothing_kdtree::Index>(il);
		}
	return kdtree_points;
	}

namespace smoothing_kdtree
{
// Use for loops that allocate inside workers: rethrow on the calling thread,
// where the exported C entry point can apply the library's error policy.
template<typename Function>
void parallel_for(const size_t count, const int chunk_size, const int number_of_threads, const Function &function)
{
    validate_number_of_threads(number_of_threads);
    std::exception_ptr failure;
    #pragma omp parallel for num_threads(number_of_threads) schedule(dynamic, chunk_size)
    for (size_t i = 0; i < count; ++i)
    {
        try { function(i); }
        catch (...)
        {
            #pragma omp critical (smoothing_worker_failure)
            {
                if (!failure) failure = std::current_exception();
            }
        }
    }
    if (failure) std::rethrow_exception(failure);
}

template<typename Function>
void overlap_accumulation_for(const size_t count, const int number_of_threads,
    const Function &function)
{
    validate_number_of_threads(number_of_threads);
    std::exception_ptr failure;
    #pragma omp parallel for num_threads(number_of_threads) schedule(static)
    for (size_t i = 0; i < count; ++i)
    {
        try { function(i); }
        catch (...)
        {
            #pragma omp critical (smoothing_worker_failure)
            {
                if (!failure) failure = std::current_exception();
            }
        }
    }
    if (failure) std::rethrow_exception(failure);
}
}

// Match the cached KD-tree frontier policy. Independent subtrees finish before
// upper nodes are combined, keeping every sum in point/left/right order.
template<typename SubtreeFunction, typename CombineFunction>
void direct_parallel_subtree_sums(const smoothing_kdtree::Node *root,
    const size_t n, const int threads, SubtreeFunction subtree, CombineFunction combine)
{
    unsigned depth=0;
    while ((size_t(1)<<depth)<size_t(threads)*4 && (n>>(depth+1))>=4096) ++depth;
    if (threads==1 || !depth) { subtree(root); return; }
    vector<const smoothing_kdtree::Node*> frontier, upper;
    vector<pair<const smoothing_kdtree::Node*,unsigned>> pending;
    if (root) pending.emplace_back(root,0);
    while (!pending.empty()) {
        const auto entry=pending.back(); pending.pop_back();
        const auto *node=entry.first;
        if (entry.second==depth) { frontier.push_back(node); continue; }
        upper.push_back(node);
        if (node->left_child()) pending.emplace_back(node->left_child(),entry.second+1);
        if (node->right_child()) pending.emplace_back(node->right_child(),entry.second+1);
    }
    smoothing_kdtree::parallel_for(frontier.size(),1,threads,[&](size_t i) {
        subtree(frontier[i]);
    });
    for (auto it=upper.rbegin(); it!=upper.rend(); ++it) combine(*it);
}

void smooth_field_using_kd_tree(const double r_kernel_in_metres, const kdtree::KdTree<3> &kdtree, const double * const lat, const double * const lon, const double * const area_size, const double * const f, const size_t number_of_points, double * const f_smoothed, const int number_of_threads, const bool show_output_log = true)
	{
    smoothing_kdtree::validate_number_of_threads(number_of_threads);

	smoothing_kdtree::validate_lat_lon(lat, lon, number_of_points);
	const double radius = cap_smoothing_radius(r_kernel_in_metres);
	if (kdtree.count_number_of_nonnull_nodes() != number_of_points) throw invalid_argument("Point count differs from smoothing KD-tree.");
	smoothing_kdtree::require_buffer(area_size, number_of_points);
	smoothing_kdtree::require_buffer(f, number_of_points);
	smoothing_kdtree::require_buffer(f_smoothed, number_of_points);
	if (number_of_points == 0) return;
	auto begin3 = std::chrono::high_resolution_clock::now();
	vector <double> fxareasize (number_of_points,0);
	#pragma omp parallel for num_threads(number_of_threads) schedule(static)
	for (unsigned long il=0; il < number_of_points; il++)
		fxareasize[il] = f[il]*area_size[il];

	vector <double> fxareasize_BB_data (number_of_points,0);
	vector <double> area_size_BB_data (number_of_points,0);

	direct_parallel_subtree_sums(kdtree.root, number_of_points, number_of_threads,
		[&](const smoothing_kdtree::Node *node) {
			generate_fxareasize_and_areasize_Bounding_Box_data(node, f, area_size,
				fxareasize_BB_data, area_size_BB_data);
		},
		[&](const smoothing_kdtree::Node *node) {
			const size_t id=node->val.index;
			double value=f[id]*area_size[id], weight=area_size[id];
			if (node->left_child()) {
				value+=fxareasize_BB_data[node->left_child()->val.index];
				weight+=area_size_BB_data[node->left_child()->val.index];
			}
			if (node->right_child()) {
				value+=fxareasize_BB_data[node->right_child()->val.index];
				weight+=area_size_BB_data[node->right_child()->val.index];
			}
			fxareasize_BB_data[id]=value; area_size_BB_data[id]=weight;
		});

	float r_Tunel_Distance = great_circle_distance_to_euclidian_distance(radius);
	float r_Tunel_Distance_squared = r_Tunel_Distance * r_Tunel_Distance;

	if (show_output_log) {
		cout << "---0 %" << endl;
	}
	size_t completed_points = 0;
	unsigned next_progress_percent = 10;
	const auto report_progress = [&]() {
		if (!show_output_log) return;
		size_t completed;
		#pragma omp atomic capture
		{ ++completed_points; completed = completed_points; }
		const size_t milestone = completed * 10 / number_of_points;
		if (milestone == (completed - 1) * 10 / number_of_points) return;
		#pragma omp critical (direct_smoothing_progress)
		{
			while (next_progress_percent <= milestone * 10) {
				cout << "---" << next_progress_percent << " %" << endl;
				next_progress_percent += 10;
			}
		}
	};

	smoothing_kdtree::parallel_for(number_of_points, 200, number_of_threads, [&](const size_t il)
		{
		double fxareasize_sum = 0;
		double areasize_sum = 0;

		double x,y,z;
		spherical_to_cartesian_coordinates(deg2rad(lat[il]), deg2rad(lon[il]), Earth_radius, x, y, z);

		smoothing_kdtree::Point point;
		point.coords = {{(float)x, (float)y, (float)z}};
		point.index = static_cast<smoothing_kdtree::Index>(il);
		get_fxareasize_and_areasize_sums_of_points_in_radius(kdtree.root, kdtree, point, r_Tunel_Distance, r_Tunel_Distance_squared, fxareasize, area_size, fxareasize_BB_data, area_size_BB_data, fxareasize_sum, areasize_sum);

		if (area_size[il] > 0 && areasize_sum > 0)
			f_smoothed[il] = fxareasize_sum/areasize_sum;
		else
			f_smoothed[il] = 0;

		report_progress();
		});
	const auto smoothing_end = std::chrono::high_resolution_clock::now();
	if (show_output_log)
		cout << "----- Smoothing "
			<< std::chrono::duration<double>(smoothing_end-begin3).count() << " s" << endl;

	}

void smooth_field_using_kd_tree_multiple_fields_simultaneously(const double r_kernel_in_metres, const kdtree::KdTree<3> &kdtree, const double * const lat, const double * const lon, const double * const area_size, const double * const * const f, const size_t number_of_points, const size_t number_of_fields, double * const * const f_smoothed, const int number_of_threads, const bool show_output_log = true)
	{
    smoothing_kdtree::validate_number_of_threads(number_of_threads);

	smoothing_kdtree::validate_point_count(number_of_points);
	const double radius = cap_smoothing_radius(r_kernel_in_metres);
	if (kdtree.count_number_of_nonnull_nodes() != number_of_points) throw invalid_argument("Point count differs from smoothing KD-tree.");
	if (number_of_points == 0 || number_of_fields == 0) return;
	smoothing_kdtree::validate_lat_lon(lat, lon, number_of_points);
	smoothing_kdtree::require_buffer(area_size, number_of_points);
	smoothing_kdtree::require_buffer(f, number_of_fields);
	smoothing_kdtree::require_buffer(f_smoothed, number_of_fields);
	for (size_t iff=0; iff<number_of_fields; ++iff)
		{
		smoothing_kdtree::require_buffer(f[iff], number_of_points);
		smoothing_kdtree::require_buffer(f_smoothed[iff], number_of_points);
		}
	auto begin3 = std::chrono::high_resolution_clock::now();
	vector <vector <double>> fxareasize (number_of_fields, vector <double> (number_of_points,0));
	for (unsigned long iff=0; iff < number_of_fields; iff++)
		{
		#pragma omp parallel for num_threads(number_of_threads) schedule(static)
		for (unsigned long il=0; il < number_of_points; il++)
			fxareasize[iff][il] = f[iff][il]*area_size[il];
		}

	vector <vector <double>> fxareasize_BB_data (number_of_fields, vector <double> (number_of_points,0));
	// Shared subtree area sums and query denominators are computed once.
	vector <double> area_size_BB_data (number_of_points,0);

	direct_parallel_subtree_sums(kdtree.root, number_of_points, number_of_threads,
		[&](const smoothing_kdtree::Node *node) {
			generate_fxareasize_and_areasize_Bounding_Box_data_multiple_fields_simultaneously(
				node, f, area_size, number_of_fields, fxareasize_BB_data, area_size_BB_data);
		},
		[&](const smoothing_kdtree::Node *node) {
			const size_t id=node->val.index;
			double weight=area_size[id];
			if (node->left_child()) weight+=area_size_BB_data[node->left_child()->val.index];
			if (node->right_child()) weight+=area_size_BB_data[node->right_child()->val.index];
			area_size_BB_data[id]=weight;
			for (size_t iff=0; iff<number_of_fields; ++iff) {
				double value=f[iff][id]*area_size[id];
				if (node->left_child()) value+=fxareasize_BB_data[iff][node->left_child()->val.index];
				if (node->right_child()) value+=fxareasize_BB_data[iff][node->right_child()->val.index];
				fxareasize_BB_data[iff][id]=value;
			}
		});

	float r_Tunel_Distance = great_circle_distance_to_euclidian_distance(radius);
	float r_Tunel_Distance_squared = r_Tunel_Distance * r_Tunel_Distance;

	if (show_output_log) {
		cout << "---0 %" << endl;
	}
	size_t completed_points = 0;
	unsigned next_progress_percent = 10;
	const auto report_progress = [&]() {
		if (!show_output_log) return;
		size_t completed;
		#pragma omp atomic capture
		{ ++completed_points; completed = completed_points; }
		const size_t milestone = completed * 10 / number_of_points;
		if (milestone == (completed - 1) * 10 / number_of_points) return;
		#pragma omp critical (direct_smoothing_progress)
		{
			while (next_progress_percent <= milestone * 10) {
				cout << "---" << next_progress_percent << " %" << endl;
				next_progress_percent += 10;
			}
		}
	};

	smoothing_kdtree::parallel_for(number_of_points, 1000, number_of_threads, [&](const size_t il)
		{
		vector <double> fxareasize_sum (number_of_fields,0);
		double areasize_sum = 0;

		double x,y,z;
		spherical_to_cartesian_coordinates(deg2rad(lat[il]), deg2rad(lon[il]), Earth_radius, x, y, z);

		smoothing_kdtree::Point point;
		point.coords = {{(float)x, (float)y, (float)z}};
		point.index = static_cast<smoothing_kdtree::Index>(il);
		get_fxareasize_and_areasize_sums_of_points_in_radius_multiple_fields_simultaneously(kdtree.root, kdtree, point, r_Tunel_Distance, r_Tunel_Distance_squared, fxareasize, area_size, fxareasize_BB_data, area_size_BB_data, fxareasize_sum, areasize_sum, number_of_fields);

		for (unsigned long iff=0; iff < number_of_fields; iff++)
			{
			if (area_size[il] > 0 && areasize_sum > 0)
				f_smoothed[iff][il] = fxareasize_sum[iff]/areasize_sum;
			else
				f_smoothed[iff][il] = 0;
			}

		report_progress();
		});
	const auto smoothing_end = std::chrono::high_resolution_clock::now();
	if (show_output_log)
		cout << "----- Smoothing " << number_of_fields << " fields simultaneously "
			<< std::chrono::duration<double>(smoothing_end-begin3).count() << " s" << endl;

	}

struct overlap_sequence_link
	{
	uint32_t current_original_index;
	uint32_t previous_sequence_index;
	};
static_assert(sizeof(overlap_sequence_link) == 2 * sizeof(uint32_t),
	"Overlap sequence links must be stored compactly.");

struct OverlapCacheData
	{
	double radius = 0;
	vector<uint32_t> spatial_to_original;
	vector<overlap_sequence_link> sequence;
	vector<unique_ptr<uint8_t[]>> records;
	};

namespace overlap_s2
	{
	const char cache_magic[8] = {'O','V','S','2','D','L','T','1'};
	const uint32_t endian_marker = 0x01020304;

	#if defined(__GNUC__) || defined(__clang__)
	__attribute__((always_inline)) inline
	#else
	inline
	#endif
	uint32_t read_u32_le(const uint8_t *data)
		{
		return uint32_t(data[0]) | (uint32_t(data[1]) << 8) |
			(uint32_t(data[2]) << 16) | (uint32_t(data[3]) << 24);
		}

	void validate_record(const uint8_t *record, const size_t size, const size_t number_of_points)
		{
		if (size < 16) throw runtime_error("Truncated overlap S2-delta record header.");
		uint32_t header[4];
		std::memcpy(header, record, sizeof(header));
		if (header[0] >= number_of_points || header[1] >= number_of_points ||
			header[2] > number_of_points || header[3] > number_of_points - header[2])
			throw runtime_error("Invalid overlap S2-delta record header.");

		size_t position = 16;
		for (size_t channel=0; channel<2; ++channel)
			{
			uint32_t used = 0;
			uint32_t previous = 0;
			bool first = true;
			const uint32_t count = header[2 + channel];
			while (used < count)
				{
				if (size - position < 5)
					throw runtime_error("Truncated overlap S2-delta block.");
				const unsigned extra = record[position++];
				const uint32_t base = read_u32_le(record + position);
				position += 4;
				if (extra >= count - used || extra > size - position ||
					base >= number_of_points || (!first && base <= previous))
					throw runtime_error("Invalid overlap S2-delta block.");
				first = false;
				uint32_t index = base;
				for (unsigned i=0; i<extra; ++i)
					{
					const unsigned delta = record[position++];
					if (delta == 0 || uint64_t(index) + delta >= number_of_points)
						throw runtime_error("Invalid overlap S2 delta.");
					index += delta;
					}
				previous = index;
				used += extra + 1;
				}
			}
		if (position != size)
			throw runtime_error("Trailing overlap S2-delta record bytes.");
		}
	}

void read_overlap_cache_data_from_binary_file(const string &cache_filename,
		OverlapCacheData *&overlap_cache_data_pointer, const bool show_output_log = true)
	{
	const auto loading_begin = std::chrono::steady_clock::now();
	if (overlap_cache_data_pointer != nullptr) throw invalid_argument("Output overlap cache data pointer must be null.");
	auto file = smoothing_io::open_cache(cache_filename, "rb", show_output_log);
	char magic[8];
	uint32_t endian = 0, file_point_count = 0;
	smoothing_io::read(file.get(), magic, 1, sizeof(magic));
	if (std::memcmp(magic, overlap_s2::cache_magic, sizeof(magic)) != 0)
		throw runtime_error("Not a version-1 S2-delta overlap cache; regenerate it.");
	smoothing_io::read(file.get(), &endian, sizeof(endian), 1);
	if (endian != overlap_s2::endian_marker)
		throw runtime_error("Incompatible overlap-cache byte order.");
	smoothing_io::read(file.get(), &file_point_count, sizeof(file_point_count), 1);
	smoothing_kdtree::validate_point_count(file_point_count);
	const uint32_t number_of_points = file_point_count;

	unique_ptr<OverlapCacheData> overlap_cache_data(new OverlapCacheData());
	smoothing_io::read(file.get(), &overlap_cache_data->radius, sizeof(overlap_cache_data->radius), 1);
	validate_cache_radius(overlap_cache_data->radius);
	overlap_cache_data->spatial_to_original.resize(number_of_points);
	smoothing_io::read(file.get(), overlap_cache_data->spatial_to_original.data(), sizeof(uint32_t), number_of_points);
	vector<uint8_t> permutation_seen(number_of_points, 0);
	for (const uint32_t original : overlap_cache_data->spatial_to_original)
		{
		if (original >= number_of_points || permutation_seen[original])
			throw runtime_error("Invalid overlap-cache S2 permutation.");
		permutation_seen[original] = 1;
		}

	overlap_cache_data->records.resize(number_of_points);
	overlap_cache_data->sequence.resize(number_of_points);
	vector<uint32_t> sequence_position_by_spatial_index(
		number_of_points, numeric_limits<uint32_t>::max());
	for (size_t i=0; i<number_of_points; ++i)
		{
		uint64_t record_size = 0;
		smoothing_io::read(file.get(), &record_size, sizeof(record_size), 1);
		const uint64_t maximum_size = 16 + uint64_t(number_of_points) * 5;
		if (record_size < 16 || record_size > maximum_size ||
			record_size > numeric_limits<size_t>::max())
			throw runtime_error("Invalid overlap-cache record size.");
		const size_t loaded_record_size = static_cast<size_t>(record_size);
		overlap_cache_data->records[i].reset(new uint8_t[loaded_record_size]);
		smoothing_io::read(file.get(), overlap_cache_data->records[i].get(), 1, loaded_record_size);
		overlap_s2::validate_record(overlap_cache_data->records[i].get(), loaded_record_size, number_of_points);
		uint32_t header[4];
		std::memcpy(header, overlap_cache_data->records[i].get(), sizeof(header));
		const uint32_t current = header[0], previous = header[1];
		if (sequence_position_by_spatial_index[current] != numeric_limits<uint32_t>::max() ||
			(i == 0 && previous != current) ||
			(i != 0 && sequence_position_by_spatial_index[previous] ==
				numeric_limits<uint32_t>::max()))
			throw runtime_error("Invalid overlap-cache sequence dependency.");
		const uint32_t previous_sequence_index = i == 0 ? 0 :
			sequence_position_by_spatial_index[previous];
		overlap_cache_data->sequence[i] = overlap_sequence_link{
			overlap_cache_data->spatial_to_original[current], previous_sequence_index};
		sequence_position_by_spatial_index[current] = static_cast<uint32_t>(i);
		}
	smoothing_io::require_end(file.get());
	smoothing_io::close(file);
	if (show_output_log)
		cout << "----- Overlap cache loading time "
			<< std::chrono::duration<double>(std::chrono::steady_clock::now() - loading_begin).count()
			<< " s" << endl;
	overlap_cache_data_pointer = overlap_cache_data.release();
	}

void free_overlap_cache_data(OverlapCacheData *&overlap_cache_data_pointer)
	{
	if (overlap_cache_data_pointer == nullptr) return;
	delete overlap_cache_data_pointer;
	overlap_cache_data_pointer = nullptr;
	}

namespace smoothing_kdtree
{
// One immutable spatial tree tracks both sets used while constructing the
// overlap sequence. Assigning a point changes only membership counters; node
// positions, child links and bounding boxes remain unchanged.
class OverlapSequenceKdTree
{
    using NodeIndex = uint32_t;
    enum : NodeIndex { invalid_node = std::numeric_limits<NodeIndex>::max() };

    struct Group
    {
        Point val;
        Coordinates coords_max{};
        Coordinates coords_min{};
        NodeIndex left = invalid_node;
        NodeIndex right = invalid_node;
        NodeIndex parent = invalid_node;
        uint32_t subtree_size = 1;
        uint32_t assigned_count = 0;
        bool assigned = false;
    };

    std::vector<Group> nodes_;
    std::vector<NodeIndex> node_from_point_index_;
    NodeIndex root_ = invalid_node;

    static float distance_squared_to_bounds(const Point &point, const Group &node) noexcept
    {
        float sum = 0;
        for (size_t axis = 0; axis < 3; ++axis)
        {
            const float nearest = std::max(node.coords_min[axis],
                                           std::min(point.coords[axis], node.coords_max[axis]));
            const float delta = point.coords[axis] - nearest;
            sum += delta * delta;
        }
        return sum;
    }

    template<bool Assigned>
    static bool subtree_has_requested_state(const Group &node) noexcept
    {
        return Assigned ? node.assigned_count != 0
                        : node.assigned_count != node.subtree_size;
    }

    void update_node_data(const NodeIndex node_index)
    {
        Group &node = nodes_[node_index];
        node.coords_min = node.val.coords;
        node.coords_max = node.val.coords;
        node.subtree_size = 1;
        const NodeIndex children[2] = {node.left, node.right};
        for (const NodeIndex child_index : children)
            if (child_index != invalid_node)
            {
                const Group &child = nodes_[child_index];
                if (node.subtree_size > std::numeric_limits<uint32_t>::max() - child.subtree_size)
                    throw std::length_error("Overlap-sequence KD-tree is too large.");
                node.subtree_size += child.subtree_size;
                for (size_t axis = 0; axis < 3; ++axis)
                {
                    node.coords_min[axis] = std::min(node.coords_min[axis], child.coords_min[axis]);
                    node.coords_max[axis] = std::max(node.coords_max[axis], child.coords_max[axis]);
                }
            }
    }

    NodeIndex build_helper(std::vector<Point> &points, const int32_t left_border,
                           const int32_t right_border, const int32_t depth,
                           const NodeIndex parent)
    {
        if (left_border > right_border) return invalid_node;
        const int32_t axis = depth % 3;
        int32_t middle = left_border + (right_border - left_border) / 2;
        std::nth_element(points.begin() + left_border, points.begin() + middle,
                         points.begin() + right_border + 1,
                         [axis](const Point &a, const Point &b)
                         { return a.coords[axis] < b.coords[axis]; });

        // Preserve the original KD-tree rule: values equal to the split value
        // belong to the node or its right subtree, never its left subtree.
        if (middle - left_border > 1)
        {
            int32_t last_taken_index = middle;
            if (points[middle - 1].coords[axis] == points[middle].coords[axis])
                --last_taken_index;
            for (int32_t i = middle - 2; i >= left_border; --i)
                if (points[i].coords[axis] == points[middle].coords[axis])
                {
                    std::swap(points[i], points[last_taken_index - 1]);
                    --last_taken_index;
                }
        }
        while (middle > left_border &&
               points[middle - 1].coords[axis] == points[middle].coords[axis])
            --middle;

        if (nodes_.size() >= static_cast<size_t>(invalid_node))
            throw std::length_error("Overlap-sequence KD-tree is too large.");
        const NodeIndex node_index = static_cast<NodeIndex>(nodes_.size());
        nodes_.emplace_back();
        nodes_[node_index].val = points[middle];
        nodes_[node_index].parent = parent;
        const size_t point_index = nodes_[node_index].val.index;
        if (point_index >= node_from_point_index_.size() ||
            node_from_point_index_[point_index] != invalid_node)
            throw std::invalid_argument("Invalid point index in overlap-sequence KD-tree.");
        node_from_point_index_[point_index] = node_index;

        nodes_[node_index].left = build_helper(points, left_border, middle - 1,
                                               depth + 1, node_index);
        nodes_[node_index].right = build_helper(points, middle + 1, right_border,
                                                depth + 1, node_index);
        update_node_data(node_index);
        return node_index;
    }

    // The caller has already established that this subtree contains at least
    // one point in the requested state. Assigned is a template parameter so
    // state selection is resolved at compile time throughout the recursion.
    template<size_t Axis, bool Assigned>
    void find_nearest_helper(const NodeIndex node_index, const Point &point,
                             float &best_distance, NodeIndex &nearest) const
    {
        static_assert(Axis < 3, "Overlap-sequence KD-tree axis is invalid.");
        const Group &node = nodes_[node_index];

        if (node.assigned == Assigned)
        {
            const float candidate_distance = distance_squared(point, node.val);
            if (candidate_distance < best_distance)
            {
                best_distance = candidate_distance;
                nearest = node_index;
            }
        }

        constexpr size_t next_axis = (Axis + 1) % 3;
        const bool left_first = point.coords[Axis] < node.val.coords[Axis];
        NodeIndex first = left_first ? node.left : node.right;
        NodeIndex second = left_first ? node.right : node.left;
        float first_bound = first != invalid_node &&
                            subtree_has_requested_state<Assigned>(nodes_[first])
                          ? distance_squared_to_bounds(point, nodes_[first])
                          : std::numeric_limits<float>::infinity();
        float second_bound = second != invalid_node &&
                             subtree_has_requested_state<Assigned>(nodes_[second])
                           ? distance_squared_to_bounds(point, nodes_[second])
                           : std::numeric_limits<float>::infinity();
        if (second_bound < first_bound)
        {
            std::swap(first, second);
            std::swap(first_bound, second_bound);
        }
        if (first != invalid_node && first_bound < best_distance)
            find_nearest_helper<next_axis, Assigned>(first, point, best_distance, nearest);
        if (second != invalid_node && second_bound < best_distance)
            find_nearest_helper<next_axis, Assigned>(second, point, best_distance, nearest);
    }

    template<bool Assigned>
    const Point *find_nearest(const Point &point) const
    {
        for (const float coordinate : point.coords)
            if (!std::isfinite(coordinate))
                throw std::invalid_argument("Overlap-sequence query coordinates must be finite.");
        if (root_ == invalid_node ||
            !subtree_has_requested_state<Assigned>(nodes_[root_]))
            return nullptr;
        float best_distance = std::numeric_limits<float>::infinity();
        NodeIndex nearest = invalid_node;
        find_nearest_helper<0, Assigned>(root_, point, best_distance, nearest);
        return nearest == invalid_node ? nullptr : &nodes_[nearest].val;
    }

    template<size_t Axis>
    void find_nearest_predecessor_helper(const NodeIndex node_index,
                                         const Point &point,
                                         const uint32_t sequence_position,
                                         float &best_distance,
                                         NodeIndex &nearest) const
    {
        static_assert(Axis < 3, "Overlap-sequence KD-tree axis is invalid.");
        const Group &node = nodes_[node_index];

        // After prepare_predecessor_queries(), assigned_count stores the
        // minimum sequence position in this subtree and subtree_size stores
        // this node's sequence position. A predecessor has position < the
        // query position.
        if (node.subtree_size < sequence_position)
        {
            const float candidate_distance = distance_squared(point, node.val);
            if (candidate_distance < best_distance)
            {
                best_distance = candidate_distance;
                nearest = node_index;
            }
        }

        constexpr size_t next_axis = (Axis + 1) % 3;
        const bool left_first = point.coords[Axis] < node.val.coords[Axis];
        NodeIndex first = left_first ? node.left : node.right;
        NodeIndex second = left_first ? node.right : node.left;
        float first_bound = first != invalid_node &&
                            nodes_[first].assigned_count < sequence_position
                          ? distance_squared_to_bounds(point, nodes_[first])
                          : std::numeric_limits<float>::infinity();
        float second_bound = second != invalid_node &&
                             nodes_[second].assigned_count < sequence_position
                           ? distance_squared_to_bounds(point, nodes_[second])
                           : std::numeric_limits<float>::infinity();
        if (second_bound < first_bound)
        {
            std::swap(first, second);
            std::swap(first_bound, second_bound);
        }
        if (first != invalid_node && first_bound < best_distance)
            find_nearest_predecessor_helper<next_axis>(
                first, point, sequence_position, best_distance, nearest);
        if (second != invalid_node && second_bound < best_distance)
            find_nearest_predecessor_helper<next_axis>(
                second, point, sequence_position, best_distance, nearest);
    }

public:
    explicit OverlapSequenceKdTree(const std::vector<Point> &points)
    {
        validate_point_count(points.size());
        nodes_.reserve(points.size());
        node_from_point_index_.assign(points.size(), invalid_node);
        if (points.empty()) return;
        std::vector<Point> working_points(points);
        root_ = build_helper(working_points, 0,
                             static_cast<int32_t>(working_points.size() - 1),
                             0, invalid_node);
    }

    const Point *find_nearest_unassigned(const Point &point) const
    {
        return find_nearest<false>(point);
    }

    void prepare_predecessor_queries(
        const std::vector<uint32_t> &sequence_position_by_point_index)
    {
        if (sequence_position_by_point_index.size() != node_from_point_index_.size())
            throw std::invalid_argument("Overlap-sequence position count mismatch.");

        // All points have now been assigned, so the membership counters are
        // no longer needed. Reuse them to avoid adding memory to every node.
        // Nodes are created before their descendants, so reverse storage order
        // computes each subtree minimum after both children are ready.
        for (Group &node : nodes_)
        {
            const size_t point_index = node.val.index;
            if (point_index >= sequence_position_by_point_index.size())
                throw std::logic_error("Invalid overlap-sequence point index.");
            node.subtree_size = sequence_position_by_point_index[point_index];
            node.assigned_count = node.subtree_size;
        }
        for (size_t i=nodes_.size(); i-- > 0;)
        {
            Group &node = nodes_[i];
            if (node.left != invalid_node)
                node.assigned_count = std::min(
                    node.assigned_count, nodes_[node.left].assigned_count);
            if (node.right != invalid_node)
                node.assigned_count = std::min(
                    node.assigned_count, nodes_[node.right].assigned_count);
        }
    }

    const Point *find_nearest_predecessor(const Point &point,
                                          const uint32_t sequence_position) const
    {
        if (sequence_position == 0)
            throw std::invalid_argument("The first sequence point has no predecessor.");
        for (const float coordinate : point.coords)
            if (!std::isfinite(coordinate))
                throw std::invalid_argument("Overlap-sequence query coordinates must be finite.");
        if (root_ == invalid_node ||
            nodes_[root_].assigned_count >= sequence_position)
            return nullptr;
        float best_distance = std::numeric_limits<float>::infinity();
        NodeIndex nearest = invalid_node;
        find_nearest_predecessor_helper<0>(
            root_, point, sequence_position, best_distance, nearest);
        return nearest == invalid_node ? nullptr : &nodes_[nearest].val;
    }

    void mark_assigned(const size_t point_index)
    {
        if (point_index >= node_from_point_index_.size())
            throw std::out_of_range("Overlap-sequence point index is outside the tree.");
        NodeIndex node_index = node_from_point_index_[point_index];
        if (node_index == invalid_node || nodes_[node_index].assigned)
            throw std::logic_error("Overlap-sequence point is missing or already assigned.");
        nodes_[node_index].assigned = true;
        while (node_index != invalid_node)
        {
            Group &node = nodes_[node_index];
            if (node.assigned_count >= node.subtree_size)
                throw std::logic_error("Invalid overlap-sequence membership count.");
            ++node.assigned_count;
            node_index = node.parent;
        }
    }

    bool has_unassigned_points() const noexcept
    {
        return root_ != invalid_node &&
               nodes_[root_].assigned_count != nodes_[root_].subtree_size;
    }

    void free_memory() noexcept
    {
        std::vector<Group>().swap(nodes_);
        std::vector<NodeIndex>().swap(node_from_point_index_);
        root_ = invalid_node;
    }
};
}

namespace saved_s2
	{
	uint64_t spherical_s2_style_key(double x, double y, double z);
	void sort_keys(vector<pair<uint64_t, uint32_t>> &keys, int threads);
	}

namespace overlap_s2
	{
	void encode_record(const uint32_t current, const uint32_t previous,
		const vector<uint32_t> &added, const vector<uint32_t> &removed,
		const size_t number_of_points, vector<uint8_t> &record)
		{
		if (added.size() > numeric_limits<uint32_t>::max() ||
			removed.size() > numeric_limits<uint32_t>::max())
			throw length_error("Overlap record contains too many point references.");
		if (std::adjacent_find(added.begin(), added.end()) != added.end() ||
			std::adjacent_find(removed.begin(), removed.end()) != removed.end())
			throw logic_error("Overlap record contains duplicate point references.");
		for (size_t added_index=0, removed_index=0;
			added_index<added.size() && removed_index<removed.size();)
			{
			if (added[added_index] == removed[removed_index])
				throw logic_error("Overlap added and removed sets are not disjoint.");
			if (added[added_index] < removed[removed_index]) ++added_index;
			else ++removed_index;
			}
		const uint32_t header[4] = {current, previous,
			static_cast<uint32_t>(added.size()), static_cast<uint32_t>(removed.size())};
		size_t size = 0;
		for (int pass=0; pass<2; ++pass)
			{
			size_t position = 16;
			uint8_t *output = pass ? record.data() : nullptr;
			if (output) std::memcpy(output, header, sizeof(header));
			for (size_t channel=0; channel<2; ++channel)
				{
				const vector<uint32_t> &indices = channel == 0 ? added : removed;
				for (size_t begin=0; begin<indices.size();)
					{
					const uint32_t base = indices[begin];
					size_t end = begin + 1;
					while (end < indices.size() && end - begin < 256 &&
						indices[end] > indices[end - 1] &&
						indices[end] - indices[end - 1] <= 255)
						++end;
					if (output)
						{
						output[position] = static_cast<uint8_t>(end - begin - 1);
						for (unsigned byte=0; byte<4; ++byte)
							output[position + 1 + byte] = static_cast<uint8_t>(base >> (8 * byte));
						for (size_t i=begin + 1; i<end; ++i)
							output[position + 4 + i - begin] =
								static_cast<uint8_t>(indices[i] - indices[i - 1]);
						}
					position += 5 + end - begin - 1;
					begin = end;
					}
				}
			if (!pass)
				{
				size = position;
				record.resize(size);
				}
			else if (position != size)
				throw runtime_error("Overlap S2-delta record size mismatch.");
			}
		validate_record(record.data(), size, number_of_points);
		}
	}

void generate_overlap_cache_data_and_write_to_disk(const double * const lat, const double * const lon, const size_t number_of_points, const vector <double> &requested_radii, const string overlap_cache_data_folder, const size_t starting_point, const int number_of_threads, const bool show_output_log = true)
	{
    smoothing_kdtree::validate_number_of_threads(number_of_threads);

	smoothing_kdtree::validate_lat_lon(lat, lon, number_of_points);
	if ((number_of_points == 0 && starting_point != 0) ||
		(number_of_points != 0 && starting_point >= number_of_points))
		throw invalid_argument("Starting point is outside the smoothing grid.");
	vector<double> smoothing_kernel_radius_in_metres;
	vector<string> cache_filenames;
	for (double radius : requested_radii)
		{
		const double capped_radius = cap_smoothing_radius(radius);
		smoothing_kernel_radius_in_metres.push_back(capped_radius);
		cache_filenames.push_back(make_overlap_cache_filename(overlap_cache_data_folder, capped_radius));
		}
	vector<string> sorted_filenames = cache_filenames;
	std::sort(sorted_filenames.begin(), sorted_filenames.end());
	if (std::adjacent_find(sorted_filenames.begin(), sorted_filenames.end()) != sorted_filenames.end())
		throw invalid_argument("Smoothing radii must have distinct overlap cache filenames after capping and rounding.");
	if (smoothing_kernel_radius_in_metres.empty()) return;
	// Create the S2 ordering before either overlap KD-tree is built.
	vector <smoothing_kdtree::Point> kdtree_points = generate_vector_of_kdtree_points_from_lat_lon_points_provided_as_arrays(lat, lon, number_of_points, number_of_threads);
	vector<pair<uint64_t, uint32_t>> spatial_keys(number_of_points);
	smoothing_kdtree::parallel_for(number_of_points, 200, number_of_threads, [&](const size_t i)
		{
		spatial_keys[i] = make_pair(saved_s2::spherical_s2_style_key(
			kdtree_points[i].coords[0], kdtree_points[i].coords[1], kdtree_points[i].coords[2]),
			static_cast<uint32_t>(i));
		});
	saved_s2::sort_keys(spatial_keys, number_of_threads);
	vector<uint32_t> spatial_to_original(number_of_points);
	vector<uint32_t> original_to_spatial(number_of_points);
	vector<smoothing_kdtree::Point> spatial_points(number_of_points);
	smoothing_kdtree::parallel_for(number_of_points, 200, number_of_threads, [&](const size_t spatial_index)
		{
		const uint32_t original_index = spatial_keys[spatial_index].second;
		spatial_to_original[spatial_index] = original_index;
		original_to_spatial[original_index] = static_cast<uint32_t>(spatial_index);
		spatial_points[spatial_index] = kdtree_points[original_index];
		spatial_points[spatial_index].index = static_cast<smoothing_kdtree::Index>(spatial_index);
		});
	kdtree_points.swap(spatial_points);
	free_vector_memory(spatial_keys);

	if (show_output_log)
		cout << "----- Generating smoothing sequence for the points " << endl;

	// ------------------------------------------------
	// Generate a smoothing sequence for the S2-ordered points
	// ------------------------------------------------

	const uint32_t sequence_starting_point_index = number_of_points == 0 ? 0 :
		original_to_spatial[static_cast<uint32_t>(starting_point)];
	free_vector_memory(original_to_spatial);

	const auto sequence_total_begin = std::chrono::high_resolution_clock::now();
	vector <uint32_t> iterative_sequence_current_smoothing_point;
	vector <uint32_t> iterative_sequence_previous_smoothing_point;
	iterative_sequence_current_smoothing_point.reserve(number_of_points);
	vector<uint32_t> sequence_position_by_point_index(
		number_of_points, std::numeric_limits<uint32_t>::max());

	uint32_t current_point_index = sequence_starting_point_index;

	// One immutable KD-tree represents both assigned and unassigned points.
	smoothing_kdtree::OverlapSequenceKdTree sequence_tree(kdtree_points);

	const size_t output_stepX = max<size_t>(1, number_of_points / 10);

	uint32_t number_of_unassigned_points=static_cast<uint32_t>(number_of_points);
	while (sequence_tree.has_unassigned_points())
		{

		if (show_output_log && number_of_unassigned_points % output_stepX == 0)
			{
			cout << "---" << round_to_digits( (double)(number_of_points - number_of_unassigned_points)/(double)number_of_points * 100, 1) << " %" << endl;
			}

		uint32_t ind1;

		// special case for the first point
		if (number_of_unassigned_points == number_of_points)
			{
			ind1 = 	current_point_index;
			}

		else
			{
			// the closest unassigned point to the current point
			auto node1 = sequence_tree.find_nearest_unassigned(kdtree_points[current_point_index]);
			if (node1 == nullptr) throw logic_error("Sequence KD-tree has no unassigned point.");
			ind1=node1->index;

			}

		// Move the point logically by updating membership counts in the same tree.
		sequence_tree.mark_assigned(ind1);

		// Add the point order. Nearest assigned predecessors are independent
		// once this complete order is known and are calculated in parallel below.
		sequence_position_by_point_index[ind1] =
			static_cast<uint32_t>(iterative_sequence_current_smoothing_point.size());
		iterative_sequence_current_smoothing_point.push_back(ind1);
		number_of_unassigned_points--;

		current_point_index = ind1;
		}

	sequence_tree.prepare_predecessor_queries(sequence_position_by_point_index);
	free_vector_memory(sequence_position_by_point_index);
	iterative_sequence_previous_smoothing_point.resize(number_of_points);
	if (number_of_points != 0)
	{
		iterative_sequence_previous_smoothing_point[0] =
			iterative_sequence_current_smoothing_point[0];
		smoothing_kdtree::parallel_for(number_of_points - 1, 200, number_of_threads,
			[&](const size_t offset)
			{
			const size_t sequence_position = offset + 1;
			const uint32_t current =
				iterative_sequence_current_smoothing_point[sequence_position];
			const auto predecessor = sequence_tree.find_nearest_predecessor(
				kdtree_points[current], static_cast<uint32_t>(sequence_position));
			if (predecessor == nullptr)
				throw std::logic_error("Overlap-sequence predecessor is missing.");
			iterative_sequence_previous_smoothing_point[sequence_position] =
				predecessor->index;
			});
	}
	const double parallel_sequence_seconds =
		std::chrono::duration_cast<std::chrono::nanoseconds>(
			std::chrono::high_resolution_clock::now() - sequence_total_begin).count() * 1e-9;
	if (show_output_log)
		cout << "----- Total smoothing-sequence preparation "
			<< parallel_sequence_seconds << " s" << endl;
	sequence_tree.free_memory();

	// Reuse the Cartesian points in a read-only KD-tree for cache differences.
	kdtree::KdTree<3> overlap_cache_tree;
	if (!kdtree_points.empty() &&
		!overlap_cache_tree.buildKdTree_and_do_not_change_the_kdtree_points_vector(kdtree_points))
		throw runtime_error("Cannot construct overlap-cache KD-tree.");

	// Generate overlap records from the completed smoothing sequence.

	if (show_output_log)
		cout << "----- Overlap detection and writing overlap cache data to disk " << endl;

	vector <float> kernel_tunnel_distance_radius_in_metres;
	for (unsigned long il=0; il < smoothing_kernel_radius_in_metres.size(); il++)
		kernel_tunnel_distance_radius_in_metres.push_back(
			static_cast<float>(great_circle_distance_to_euclidian_distance(
				static_cast<double>(smoothing_kernel_radius_in_metres[il]))));

	vector<smoothing_io::File> pFile;
	pFile.reserve(smoothing_kernel_radius_in_metres.size());

	for (unsigned long ir=0; ir < smoothing_kernel_radius_in_metres.size(); ir++)
		{
		auto file = smoothing_io::open_cache(cache_filenames[ir], "wb", show_output_log);
		const uint32_t file_point_count = static_cast<uint32_t>(number_of_points);
		smoothing_io::write(file.get(), overlap_s2::cache_magic, 1, sizeof(overlap_s2::cache_magic));
		smoothing_io::write(file.get(), &overlap_s2::endian_marker, sizeof(uint32_t), 1);
		smoothing_io::write(file.get(), &file_point_count, sizeof(uint32_t), 1);
		smoothing_io::write(file.get(), &smoothing_kernel_radius_in_metres[ir], sizeof(double), 1);
		smoothing_io::write(file.get(), spatial_to_original.data(), sizeof(uint32_t), number_of_points);
		pFile.push_back(std::move(file));
		}

	const size_t radius_count = smoothing_kernel_radius_in_metres.size();
	// Large reusable batches substantially reduce OpenMP-region, container and
	// stdio-call overhead without changing record order or encoded bytes.
	const size_t cache_block_size = 8192;
	const size_t block_capacity = min(cache_block_size,
		iterative_sequence_current_smoothing_point.size());
	vector<vector<vector<uint8_t>>> block_records(
		block_capacity, vector<vector<uint8_t>>(radius_count));
	vector<vector<uint8_t>> compressed_write_buffers(radius_count);
	const auto overlap_records_begin = std::chrono::high_resolution_clock::now();
	unsigned next_cache_progress_percent = 10;
	if (show_output_log && !iterative_sequence_current_smoothing_point.empty())
		cout << "---0 %" << endl;
	for (size_t block_begin=0; block_begin < iterative_sequence_current_smoothing_point.size();
		 block_begin += cache_block_size)
		{
		const size_t block_end = min(iterative_sequence_current_smoothing_point.size(),
			block_begin + cache_block_size);
		const size_t block_count = block_end - block_begin;
		const ptrdiff_t parallel_block_count = static_cast<ptrdiff_t>(block_count);
		exception_ptr failure;

		#pragma omp parallel for num_threads(number_of_threads) schedule(static)
		for (ptrdiff_t offset=0; offset < parallel_block_count; ++offset)
			{
			try
				{
				const size_t local_index = static_cast<size_t>(offset);
				const size_t sequence_index = block_begin + local_index;
				const uint32_t current_smoothing_point =
					iterative_sequence_current_smoothing_point[sequence_index];
				const uint32_t previous_smoothing_point =
					iterative_sequence_previous_smoothing_point[sequence_index];
				for (size_t ir=0; ir < radius_count; ++ir)
					{
					vector<uint32_t> added;
					vector<uint32_t> removed;
					if (sequence_index == 0)
						{
						const auto cluster = overlap_cache_tree.findNearestNodeCluster(
							kdtree_points[current_smoothing_point],
							kernel_tunnel_distance_radius_in_metres[ir]);
						added.reserve(cluster.size());
						for (const auto *node : cluster)
							added.push_back(node->val.index);
						}
					else
						{
						auto differences = overlap_cache_tree.findNodeClusterDifferences(
							kdtree_points[current_smoothing_point],
							kernel_tunnel_distance_radius_in_metres[ir],
							kdtree_points[previous_smoothing_point],
							kernel_tunnel_distance_radius_in_metres[ir]);
						added = std::move(differences.first);
						removed = std::move(differences.second);
						}
					std::sort(added.begin(), added.end());
					std::sort(removed.begin(), removed.end());
					overlap_s2::encode_record(
						current_smoothing_point, previous_smoothing_point,
						added, removed, number_of_points, block_records[local_index][ir]);
					}
				}
			catch (...)
				{
				#pragma omp critical (overlap_cache_worker_failure)
					{
					if (!failure) failure = current_exception();
					}
				}
			}
		if (failure) rethrow_exception(failure);

		// Assemble each radius into one contiguous byte stream. The size prefix
		// and payload bytes remain in exactly the same dependency order.
		for (size_t ir=0; ir < radius_count; ++ir)
			{
			size_t bytes = 0;
			for (size_t local_index=0; local_index < block_count; ++local_index)
				{
				const size_t record_size = block_records[local_index][ir].size();
				if (bytes > numeric_limits<size_t>::max() - sizeof(uint64_t) ||
					record_size > numeric_limits<size_t>::max() - sizeof(uint64_t) - bytes)
					throw length_error("Overlap write batch is too large.");
				bytes += sizeof(uint64_t) + record_size;
				}
			vector<uint8_t> &write_buffer = compressed_write_buffers[ir];
			write_buffer.resize(bytes);
			size_t position = 0;
			for (size_t local_index=0; local_index < block_count; ++local_index)
				{
				const vector<uint8_t> &record = block_records[local_index][ir];
				const uint64_t record_size = record.size();
				std::memcpy(write_buffer.data() + position, &record_size, sizeof(record_size));
				position += sizeof(record_size);
				if (!record.empty())
					std::memcpy(write_buffer.data() + position, record.data(), record.size());
				position += record.size();
				}
			smoothing_io::write(pFile[ir].get(), write_buffer.data(), 1, write_buffer.size());
			}

		for (size_t local_index=0; local_index < block_count; ++local_index)
			{
			const size_t sequence_index = block_begin + local_index;
			const size_t completed_records = sequence_index + 1;
			while (next_cache_progress_percent <= 100 &&
				completed_records * 100 >=
				iterative_sequence_current_smoothing_point.size() * next_cache_progress_percent)
				{
				if (show_output_log)
					cout << "---" << next_cache_progress_percent << " %" << endl;
				next_cache_progress_percent += 10;
				}
			}
		}
	if (show_output_log)
		cout << "----- Overlap detection and writing cache to disk time "
			<< std::chrono::duration_cast<std::chrono::nanoseconds>(
				std::chrono::high_resolution_clock::now() - overlap_records_begin).count() * 1e-9
			<< " s" << endl;

	for (unsigned long ir=0; ir < smoothing_kernel_radius_in_metres.size(); ir++)
		smoothing_io::close(pFile[ir]);
	overlap_cache_tree.free_memory();
	free_vector_memory(kdtree_points);
	}

#if defined(__GNUC__) || defined(__clang__)
#define SMOOTHING_NOINLINE __attribute__((noinline))
#else
#define SMOOTHING_NOINLINE
#endif

// Keep the compressed-record accumulation out of the OpenMP worker.  With GCC,
// inlining this loop into the outlined worker can spill both running sums to the
// stack for every reference instead of retaining them in floating-point registers.
static SMOOTHING_NOINLINE void accumulate_compressed_overlap_record(
	const uint8_t * const record,
	const double * const f_x_area, const double * const spatial_area,
	const size_t sequence_index,
	double * const f_x_area_sum, double * const area_sum)
	{
	uint32_t header[4];
	std::memcpy(header, record, sizeof(header));
	size_t position = 16;
	double partial_f_x_area_sum = 0;
	double partial_area_sum = 0;
	uint32_t used = 0;
	while (used < header[2])
		{
		const unsigned extra = record[position++];
		uint32_t index = overlap_s2::read_u32_le(record + position);
		position += 4;
		partial_f_x_area_sum += f_x_area[index];
		partial_area_sum += spatial_area[index];
		for (unsigned i=0; i<extra; ++i)
			{
			index += record[position++];
			partial_f_x_area_sum += f_x_area[index];
			partial_area_sum += spatial_area[index];
			}
		used += extra + 1;
		}
	used = 0;
	while (used < header[3])
		{
		const unsigned extra = record[position++];
		uint32_t index = overlap_s2::read_u32_le(record + position);
		position += 4;
		partial_f_x_area_sum -= f_x_area[index];
		partial_area_sum -= spatial_area[index];
		for (unsigned i=0; i<extra; ++i)
			{
			index += record[position++];
			partial_f_x_area_sum -= f_x_area[index];
			partial_area_sum -= spatial_area[index];
			}
		used += extra + 1;
		}
	f_x_area_sum[sequence_index] = partial_f_x_area_sum;
	area_sum[sequence_index] = partial_area_sum;
	}

// Subsequent fields reuse the fully propagated area denominator from field one.
static SMOOTHING_NOINLINE void accumulate_compressed_overlap_record_numerator(
	const uint8_t * const record, const double * const weighted_field,
	const size_t sequence_index, double * const field_sums)
	{
	uint32_t header[4];
	std::memcpy(header, record, sizeof(header));
	size_t position = 16;
	double sum = 0;
	uint32_t used = 0;
	while (used < header[2])
		{
		const unsigned extra = record[position++];
		uint32_t index = overlap_s2::read_u32_le(record + position);
		position += 4;
		sum += weighted_field[index];
		for (unsigned i=0; i<extra; ++i)
			{
			index += record[position++];
			sum += weighted_field[index];
			}
		used += extra + 1;
		}
	used = 0;
	while (used < header[3])
		{
		const unsigned extra = record[position++];
		uint32_t index = overlap_s2::read_u32_le(record + position);
		position += 4;
		sum -= weighted_field[index];
		for (unsigned i=0; i<extra; ++i)
			{
			index += record[position++];
			sum -= weighted_field[index];
			}
		used += extra + 1;
		}
	field_sums[sequence_index] = sum;
	}

#undef SMOOTHING_NOINLINE

void smooth_field_using_overlap_detection(const double * const area_size, const double * const f,
	const size_t number_of_points, const OverlapCacheData * const overlap_cache_data_pointer,
	double * const f_smoothed, const int number_of_threads, const bool show_output_log = true)
	{
    smoothing_kdtree::validate_number_of_threads(number_of_threads);

	smoothing_kdtree::validate_point_count(number_of_points);
	smoothing_kdtree::require_buffer(area_size, number_of_points);
	smoothing_kdtree::require_buffer(f, number_of_points);
	smoothing_kdtree::require_buffer(f_smoothed, number_of_points);
	smoothing_kdtree::require_buffer(overlap_cache_data_pointer, 1);
	if (overlap_cache_data_pointer->records.size() != number_of_points ||
		overlap_cache_data_pointer->sequence.size() != number_of_points ||
		overlap_cache_data_pointer->spatial_to_original.size() != number_of_points)
		throw invalid_argument("Overlap-cache point count mismatch.");
	if (number_of_points == 0) return;
	// Preparation and accumulation overwrite every element before it is read.
	std::unique_ptr<double[]> f_x_area(new double[number_of_points]);
	std::unique_ptr<double[]> spatial_area(new double[number_of_points]);
	std::unique_ptr<double[]> f_x_area_sum(new double[number_of_points]);
	std::unique_ptr<double[]> area_sum(new double[number_of_points]);
	#pragma omp parallel for num_threads(number_of_threads) schedule(static)
	for (size_t spatial_index=0; spatial_index<number_of_points; ++spatial_index)
		{
		const uint32_t original_index = overlap_cache_data_pointer->spatial_to_original[spatial_index];
		spatial_area[spatial_index] = area_size[original_index];
		f_x_area[spatial_index] = f[original_index] * spatial_area[spatial_index];
		}

	const auto begin4 = std::chrono::steady_clock::now();

	// Decode independent added/removed differences in parallel.
	smoothing_kdtree::overlap_accumulation_for(number_of_points, number_of_threads,
		[&](const size_t sequence_index)
		{
		accumulate_compressed_overlap_record(
			overlap_cache_data_pointer->records[sequence_index].get(),
			f_x_area.get(), spatial_area.get(),
			sequence_index,
			f_x_area_sum.get(), area_sum.get());
		});

	// Apply sequence dependencies serially.
	for (size_t sequence_index=0; sequence_index<number_of_points; ++sequence_index)
		{
		const overlap_sequence_link link = overlap_cache_data_pointer->sequence[sequence_index];
		const uint32_t previous_sequence_index = link.previous_sequence_index;
		if (sequence_index != previous_sequence_index)
			{
			f_x_area_sum[sequence_index] += f_x_area_sum[previous_sequence_index];
			area_sum[sequence_index] += area_sum[previous_sequence_index];
			}
		}

	// Division, scatter and missing-value handling are independent after propagation.
	#pragma omp parallel for num_threads(number_of_threads) schedule(static)
	for (size_t sequence_index=0; sequence_index<number_of_points; ++sequence_index)
		{
		const uint32_t original_index =
			overlap_cache_data_pointer->sequence[sequence_index].current_original_index;
		if (area_size[original_index] != 0 && area_sum[sequence_index] > 0)
			f_smoothed[original_index] =
				f_x_area_sum[sequence_index]/area_sum[sequence_index];
		else
			f_smoothed[original_index] = 0;
		}

	const auto output_end = std::chrono::steady_clock::now();
	if (show_output_log)
		cout << "----- Smoothing "
			<< std::chrono::duration<double>(output_end - begin4).count() << " s" << endl;
	}

void smooth_multiple_fields_simultaneously_using_overlap_detection(
	const double * const area_size, const double * const fields,
	const size_t number_of_points, const size_t number_of_fields,
	const OverlapCacheData * const overlap_cache_data_pointer,
	double * const smoothed_fields, const int number_of_threads, const bool show_output_log = true)
	{
	smoothing_kdtree::validate_number_of_threads(number_of_threads);
	smoothing_kdtree::validate_point_count(number_of_points);
	if (number_of_fields == 0)
		throw invalid_argument("At least one field is required.");
	smoothing_kdtree::require_buffer(overlap_cache_data_pointer, 1);
	if (overlap_cache_data_pointer->records.size() != number_of_points ||
		overlap_cache_data_pointer->sequence.size() != number_of_points ||
		overlap_cache_data_pointer->spatial_to_original.size() != number_of_points)
		throw invalid_argument("Overlap-cache point count mismatch.");
	if (number_of_points == 0) return;
	if (number_of_fields > numeric_limits<size_t>::max() / number_of_points)
		throw length_error("Multiple-field array dimensions are too large.");
	const size_t value_count = number_of_points * number_of_fields;
	smoothing_kdtree::require_buffer(area_size, number_of_points);
	smoothing_kdtree::require_buffer(fields, value_count);
	smoothing_kdtree::require_buffer(smoothed_fields, value_count);

	const auto begin = std::chrono::steady_clock::now();
	// Every element is written before it is read; no zero-fill or per-field allocation.
	std::unique_ptr<double[]> spatial_area(new double[number_of_points]);
	std::unique_ptr<double[]> weighted_field(new double[number_of_points]);
	std::unique_ptr<double[]> field_sums(new double[number_of_points]);
	std::unique_ptr<double[]> denominator(new double[number_of_points]);
	#pragma omp parallel for num_threads(number_of_threads) schedule(static)
	for (size_t spatial_index=0; spatial_index<number_of_points; ++spatial_index)
		spatial_area[spatial_index] =
			area_size[overlap_cache_data_pointer->spatial_to_original[spatial_index]];

	for (size_t field_index=0; field_index<number_of_fields; ++field_index)
		{
		const double * const input = fields + field_index * number_of_points;
		double * const output = smoothed_fields + field_index * number_of_points;
		#pragma omp parallel for num_threads(number_of_threads) schedule(static)
		for (size_t spatial_index=0; spatial_index<number_of_points; ++spatial_index)
			weighted_field[spatial_index] =
				input[overlap_cache_data_pointer->spatial_to_original[spatial_index]] * spatial_area[spatial_index];
		if (field_index == 0)
			{
			smoothing_kdtree::overlap_accumulation_for(number_of_points, number_of_threads,
				[&](const size_t sequence_index)
				{
				accumulate_compressed_overlap_record(overlap_cache_data_pointer->records[sequence_index].get(),
					weighted_field.get(), spatial_area.get(), sequence_index,
					field_sums.get(), denominator.get());
				});
			}
		else
			{
			smoothing_kdtree::overlap_accumulation_for(number_of_points, number_of_threads,
				[&](const size_t sequence_index)
				{
				accumulate_compressed_overlap_record_numerator(
					overlap_cache_data_pointer->records[sequence_index].get(), weighted_field.get(),
					sequence_index, field_sums.get());
				});
			}
		if (field_index == 0)
			{
			for (size_t sequence_index=0; sequence_index<number_of_points; ++sequence_index)
				{
				const uint32_t previous = overlap_cache_data_pointer->sequence[sequence_index].previous_sequence_index;
				if (sequence_index != previous)
					{
					field_sums[sequence_index] += field_sums[previous];
					denominator[sequence_index] += denominator[previous];
					}
				}
			}
		else
			{
			for (size_t sequence_index=0; sequence_index<number_of_points; ++sequence_index)
				{
				const uint32_t previous = overlap_cache_data_pointer->sequence[sequence_index].previous_sequence_index;
				if (sequence_index != previous)
					field_sums[sequence_index] += field_sums[previous];
				}
			}
		#pragma omp parallel for num_threads(number_of_threads) schedule(static)
		for (size_t sequence_index=0; sequence_index<number_of_points; ++sequence_index)
			{
			const uint32_t original = overlap_cache_data_pointer->sequence[sequence_index].current_original_index;
			output[original] = area_size[original] != 0 && denominator[sequence_index] > 0 ?
				field_sums[sequence_index] / denominator[sequence_index] : 0;
			}
		}
	const auto end = std::chrono::steady_clock::now();
	if (show_output_log)
		cout << "----- Smoothing " << number_of_fields << " fields simultaneously "
			<< std::chrono::duration<double>(end - begin).count() << " s" << endl;
	}

// S2 ordering with collapsed subtree and individual-point reference channels.
// Blocks store a 32-bit base index and up to 255 subsequent 8-bit deltas.
// Outputs retain input order; subtree grouping can change floating-point rounding.
struct KdTreeCacheData {
    kdtree::KdTree<3> kdtree;
    vector<uint32_t> spatial_to_original;
    vector<unique_ptr<uint8_t[]>> records;
    vector<size_t> record_sizes;
    double radius = 0;
};

namespace saved_s2 {
// Skilling's public-domain 2D axes-to-transpose transform.
// Programming the Hilbert curve (2004), doi:10.1063/1.1751381.
// S2-style cube-face projection and spatial ordering.
// https://s2geometry.io/devguide/s2cell_hierarchy
uint64_t hilbert2d_key(uint32_t x, uint32_t y, unsigned bits = 21) {
    if (!bits || bits > 21) throw invalid_argument("Invalid 2D Hilbert precision.");
    uint32_t a[2] = {x, y};
    const uint32_t top = uint32_t(1) << (bits - 1);
    for (uint32_t q = top; q > 1; q >>= 1) {
        for (unsigned axis = 0; axis < 2; ++axis) {
            if (a[axis] & q) a[0] ^= q - 1;
            else {
                const uint32_t exchange = (a[0] ^ a[axis]) & (q - 1);
                a[0] ^= exchange;
                a[axis] ^= exchange;
            }
        }
    }
    a[1] ^= a[0];
    uint32_t correction = 0;
    for (uint32_t q = top; q > 1; q >>= 1)
        if (a[1] & q) correction ^= q - 1;
    a[0] ^= correction;
    a[1] ^= correction;
    uint64_t key = 0;
    for (int bit = static_cast<int>(bits) - 1; bit >= 0; --bit) {
        key = (key << 1) | ((a[0] >> bit) & 1u);
        key = (key << 1) | ((a[1] >> bit) & 1u);
    }
    return key;
}
uint64_t spherical_s2_style_key(double x, double y, double z) {
    const double p[3] = {x, y, z};
    unsigned axis = 0;
    for (unsigned k = 0; k < 3; ++k) {
        if (!std::isfinite(p[k])) throw invalid_argument("Non-finite spherical ordering coordinate.");
        if (std::abs(p[k]) > std::abs(p[axis])) axis = k;
    }
    if (p[axis] == 0) throw invalid_argument("Zero direction in spherical ordering.");
    const unsigned face = axis + (p[axis] < 0 ? 3 : 0);
    // Local face basis vectors expressed as Cartesian axis/sign pairs.
    const unsigned u_axis[6] = {1, 0, 0, 2, 2, 1};
    const unsigned v_axis[6] = {2, 2, 1, 1, 0, 0};
    const int u_sign[6] = {1, -1, -1, -1, -1, 1};
    const int v_sign[6] = {1, 1, -1, -1, 1, 1};
    const double scale = std::abs(p[axis]);
    const double uv[2] = {u_sign[face] * p[u_axis[face]] / scale,
                          v_sign[face] * p[v_axis[face]] / scale};
    uint32_t ij[2];
    for (unsigned k = 0; k < 2; ++k) {
        const double u = std::max(-1.0, std::min(1.0, uv[k]));
        const double st = u >= 0 ? 0.5 * std::sqrt(1 + 3 * u)
                                 : 1 - 0.5 * std::sqrt(1 - 3 * u);
        ij[k] = static_cast<uint32_t>(std::min(2097151.0,
                   std::max(0.0, st) * 2097152.0));
    }
    if (face & 1u) std::swap(ij[0], ij[1]);
    return (uint64_t(face) << 42) | hilbert2d_key(ij[0], ij[1]);
}

unique_ptr<uint8_t[]> encode_record(const vector<uint32_t> &ids, const uint32_t *r,
                                   size_t &size) {
    unique_ptr<uint8_t[]> record;
    // First pass sizes; second writes. Every target owns a separate allocation.
    for (int pass = 0; pass < 2; ++pass) {
        size_t pos = 8;
        uint8_t *out = pass ? record.get() : nullptr;
        if (out) std::memcpy(out, r, 8);
        for (size_t channel = 0, begin = 0; channel < 2; ++channel) {
            const size_t end = begin + r[channel];
            for (size_t j = begin; j < end;) {
                const uint32_t base = ids[j];
                size_t stop = j + 1;
                // Count byte stores length minus one: 0..255 means 1..256 references.
                while (stop < end && stop - j < 256 && ids[stop] - ids[stop - 1] <= 255) ++stop;
                if (out) {
                    out[pos] = static_cast<uint8_t>(stop - j - 1);
                    for (unsigned b = 0; b < 4; ++b) out[pos + 1 + b] = static_cast<uint8_t>(base >> (8 * b));
                    for (size_t k = j + 1; k < stop; ++k) out[pos + 5 + k - j - 1] = static_cast<uint8_t>(ids[k] - ids[k - 1]);
                }
                pos += 5 + stop - j - 1;
                j = stop;
            }
            begin = end;
        }
        if (!pass) { size = pos; record.reset(new uint8_t[pos]); }
        else if (pos != size) throw runtime_error("Spatial record size mismatch.");
    }
    // Validate all decoded indices and both channel boundaries, including duplicates.
    size_t pos = 8, j = 0;
    for (size_t channel = 0; channel < 2; ++channel) {
        const size_t end = j + r[channel];
        while (j < end) {
            if (size - pos < 5) throw runtime_error("Truncated spatial block.");
            const unsigned extra = record[pos++];
            uint32_t base = 0;
            for (unsigned b = 0; b < 4; ++b) base |= uint32_t(record[pos++]) << (8 * b);
            if (extra >= end - j || extra > size - pos || ids[j++] != base)
                throw runtime_error("Spatial block validation failed.");
            uint32_t previous = base;
            for (unsigned k = 0; k < extra; ++k) {
                const uint32_t decoded = previous + record[pos++];
                if (ids[j++] != decoded) throw runtime_error("Spatial decoded index mismatch.");
                previous = decoded;
            }
        }

    }
    if (pos != size) throw runtime_error("Trailing spatial bytes.");
    return record;
}
bool collect_collapsed(const smoothing_kdtree::Node *node,
    const smoothing_kdtree::Point &point, float radius_squared,
    vector<uint32_t> &branches, vector<uint32_t> &singles) {
    if (!node) return true;
    if (!smoothing_kdtree::bounds_intersect(point.coords, radius_squared,
                                           node->coords_min, node->coords_max)) return false;
    const bool inside = smoothing_kdtree::distance_squared(point, node->val) <= radius_squared;
    if (inside && smoothing_kdtree::bounds_inside(point.coords, radius_squared,
                                                 node->coords_min, node->coords_max)) {
        branches.push_back(static_cast<uint32_t>(node->val.index));
        return true;
    }
    const size_t branch_start = branches.size(), single_start = singles.size();
    if (inside) {
        singles.push_back(static_cast<uint32_t>(node->val.index));
    }
    // Evaluate BOTH children, even when the point or the first child is outside.
    const bool right = collect_collapsed(node->right_child(), point, radius_squared, branches, singles);
    const bool left = collect_collapsed(node->left_child(), point, radius_squared, branches, singles);
    if (inside && right && left) {
        branches.resize(branch_start);
        singles.resize(single_start);
        branches.push_back(static_cast<uint32_t>(node->val.index));
        return true;
    }
    return false;
}
void native_subtree_sums(const smoothing_kdtree::Node *node,
    const double *products, const double *area,
    double *sums, double *areas) {
    if (!node) return;
    const size_t id = node->val.index;
    double value = products[id], weight = area[id];
    if (node->left_child()) {
        native_subtree_sums(node->left_child(), products, area, sums, areas);
        value += sums[node->left_child()->val.index];
        weight += areas[node->left_child()->val.index];
    }
    if (node->right_child()) {
        native_subtree_sums(node->right_child(), products, area, sums, areas);
        value += sums[node->right_child()->val.index];
        weight += areas[node->right_child()->val.index];
    }
    sums[id] = value;
    areas[id] = weight;
}

// Independent chunks followed by out-of-place merge rounds. Pair ordering uses
// original point ID to break equal spatial keys, regardless of thread count.
void sort_keys(vector<pair<uint64_t, uint32_t>> &keys, int threads) {
    const size_t n=keys.size();
    if (threads==1 || n<8192) { std::sort(keys.begin(), keys.end()); return; }
    const size_t chunks=std::min<size_t>(threads, (n+4095)/4096);
    const size_t width=(n+chunks-1)/chunks;
    smoothing_kdtree::parallel_for(chunks, 1, threads, [&](size_t c) {
        const size_t begin=c*width, end=std::min(n, begin+width);
        if (begin<end) std::sort(keys.begin()+begin, keys.begin()+end);
    });
    vector<pair<uint64_t, uint32_t>> temporary(n);
    for (size_t span=width; span<n; span*=2) {
        const size_t runs=(n+2*span-1)/(2*span);
        const size_t pieces=std::min<size_t>(threads,(2*span+4095)/4096);
        smoothing_kdtree::parallel_for(runs*pieces, 1, threads, [&](size_t task) {
            const size_t r=task/pieces, part=task%pieces;
            const size_t begin=r*2*span, mid=std::min(n,begin+span), end=std::min(n,begin+2*span);
            const size_t na=mid-begin, nb=end-mid, length=end-begin;
            const size_t lo=length*part/pieces, hi=length*(part+1)/pieces;
            // Co-rank an output boundary in the two sorted runs. This splits
            // even the final merge across workers, with disjoint output ranges.
            const auto cut=[&](size_t k) {
                size_t low=k>nb ? k-nb : 0, high=std::min(k,na);
                while (low<=high) {
                    const size_t a=low+(high-low)/2, b=k-a;
                    if (a && b<nb && keys[begin+a-1]>keys[mid+b]) high=a-1;
                    else if (b && a<na && keys[mid+b-1]>keys[begin+a]) low=a+1;
                    else return a;
                }
                throw runtime_error("Invalid S2 merge boundary.");
            };
            const size_t a0=cut(lo), a1=cut(hi);
            std::merge(keys.begin()+begin+a0,keys.begin()+begin+a1,
                       keys.begin()+mid+lo-a0,keys.begin()+mid+hi-a1,
                       temporary.begin()+begin+lo);
        });
        keys.swap(temporary);
    }
}

// Split only upper levels. Each worker computes a disjoint subtree; combining
// upper nodes afterward keeps precisely the serial point/left/right sum order.
void parallel_sums(const smoothing_kdtree::Node *root, size_t n,
    const double *products, const double *weights,
    double *sums, double *areas, int threads) {
    unsigned depth=0;
    while ((size_t(1)<<depth)<size_t(threads)*4 && (n>>(depth+1))>=4096) ++depth;
    if (threads==1 || !depth) { native_subtree_sums(root,products,weights,sums,areas); return; }
    vector<const smoothing_kdtree::Node*> frontier, upper;
    vector<pair<const smoothing_kdtree::Node*,unsigned>> pending;
    if (root) pending.emplace_back(root,0);
    while (!pending.empty()) {
        const auto entry=pending.back(); pending.pop_back();
        const auto *node=entry.first;
        if (entry.second==depth) { frontier.push_back(node); continue; }
        upper.push_back(node);
        if (node->left_child()) pending.emplace_back(node->left_child(),entry.second+1);
        if (node->right_child()) pending.emplace_back(node->right_child(),entry.second+1);
    }
    smoothing_kdtree::parallel_for(frontier.size(),1,threads,[&](size_t i) {
        native_subtree_sums(frontier[i],products,weights,sums,areas);
    });
    for (auto it=upper.rbegin(); it!=upper.rend(); ++it) {
        const auto *node=*it;
        const size_t id=node->val.index;
        double value=products[id], weight=weights[id];
        if (node->left_child()) { value+=sums[node->left_child()->val.index]; weight+=areas[node->left_child()->val.index]; }
        if (node->right_child()) { value+=sums[node->right_child()->val.index]; weight+=areas[node->right_child()->val.index]; }
        sums[id]=value; areas[id]=weight;
    }
}

// Prepare the independent subtrees once for a batch of fields. The traversal and
// upper-node combination preserve the single-field point/left/right sum order.
struct BatchSubtreeSumPlan {
    const smoothing_kdtree::Node *root;
    vector<const smoothing_kdtree::Node*> frontier, upper;
};

BatchSubtreeSumPlan prepare_batch_subtree_sum_plan(
    const smoothing_kdtree::Node *root,size_t n,int threads) {
    BatchSubtreeSumPlan plan;
    plan.root=root;
    unsigned depth=0;
    while ((size_t(1)<<depth)<size_t(threads)*4 && (n>>(depth+1))>=4096) ++depth;
    if (threads==1 || !depth) return plan;
    vector<pair<const smoothing_kdtree::Node*,unsigned>> pending;
    if (root) pending.emplace_back(root,0);
    while (!pending.empty()) {
        const auto entry=pending.back(); pending.pop_back();
        const auto *node=entry.first;
        if (entry.second==depth) { plan.frontier.push_back(node); continue; }
        plan.upper.push_back(node);
        if (node->left_child()) pending.emplace_back(node->left_child(),entry.second+1);
        if (node->right_child()) pending.emplace_back(node->right_child(),entry.second+1);
    }
    return plan;
}

void batch_subtree_sums(const BatchSubtreeSumPlan &plan,
    const double *products,const double *weights,double *sums,double *areas,int threads) {
    if (plan.upper.empty()) { native_subtree_sums(plan.root,products,weights,sums,areas); return; }
    smoothing_kdtree::parallel_for(plan.frontier.size(),1,threads,[&](size_t i) {
        native_subtree_sums(plan.frontier[i],products,weights,sums,areas);
    });
    for (auto it=plan.upper.rbegin();it!=plan.upper.rend();++it) {
        const auto *node=*it;
        const size_t id=node->val.index;
        double value=products[id],weight=weights[id];
        if (node->left_child()) { value+=sums[node->left_child()->val.index]; weight+=areas[node->left_child()->val.index]; }
        if (node->right_child()) { value+=sums[node->right_child()->val.index]; weight+=areas[node->right_child()->val.index]; }
        sums[id]=value; areas[id]=weight;
    }
}

void native_subtree_field_sums(const smoothing_kdtree::Node *node,
    const double *products,double *sums) {
    if (!node) return;
    const size_t id=node->val.index;
    double value=products[id];
    if (node->left_child()) {
        native_subtree_field_sums(node->left_child(),products,sums);
        value+=sums[node->left_child()->val.index];
    }
    if (node->right_child()) {
        native_subtree_field_sums(node->right_child(),products,sums);
        value+=sums[node->right_child()->val.index];
    }
    sums[id]=value;
}

void batch_subtree_field_sums(const BatchSubtreeSumPlan &plan,
    const double *products,double *sums,int threads) {
    if (plan.upper.empty()) { native_subtree_field_sums(plan.root,products,sums); return; }
    smoothing_kdtree::parallel_for(plan.frontier.size(),1,threads,[&](size_t i) {
        native_subtree_field_sums(plan.frontier[i],products,sums);
    });
    for (auto it=plan.upper.rbegin();it!=plan.upper.rend();++it) {
        const auto *node=*it;
        const size_t id=node->val.index;
        double value=products[id];
        if (node->left_child()) value+=sums[node->left_child()->val.index];
        if (node->right_child()) value+=sums[node->right_child()->val.index];
        sums[id]=value;
    }
}

unique_ptr<KdTreeCacheData> generate_cache_data_in_memory(
    const double *lat,const double *lon,size_t n,double radius,int threads, const bool show_output_log = true) {
    const auto logging_begin = std::chrono::steady_clock::now();
    smoothing_kdtree::validate_number_of_threads(threads);
    radius=cap_smoothing_radius(radius);
    auto points=generate_vector_of_kdtree_points_from_lat_lon_points_provided_as_arrays(lat,lon,n,threads);
    unique_ptr<KdTreeCacheData> kdtree_cache_data(new KdTreeCacheData());
    kdtree_cache_data->radius=radius;
    vector<pair<uint64_t,uint32_t>> keys(n);
    smoothing_kdtree::parallel_for(n,200,threads,[&](size_t i) {
        keys[i]=std::make_pair(spherical_s2_style_key(points[i].coords[0],points[i].coords[1],points[i].coords[2]),
                               static_cast<uint32_t>(i));
    });
    sort_keys(keys,threads);
    kdtree_cache_data->spatial_to_original.resize(n);
    {
        vector<smoothing_kdtree::Point> reordered(n);
        smoothing_kdtree::parallel_for(n,200,threads,[&](size_t i) {
            const uint32_t old=keys[i].second;
            kdtree_cache_data->spatial_to_original[i]=old;
            reordered[i]=points[old];
            reordered[i].index=static_cast<smoothing_kdtree::Index>(i);
        });
        points.swap(reordered);
    }
    vector<pair<uint64_t,uint32_t>>().swap(keys);
    if (n && !kdtree_cache_data->kdtree.buildKdTree_and_do_not_change_the_kdtree_points_vector(points))
        throw runtime_error("Cannot build KD-tree cache tree.");
    kdtree_cache_data->records.resize(n); kdtree_cache_data->record_sizes.resize(n);
    const float distance=static_cast<float>(great_circle_distance_to_euclidian_distance(radius));
    const float squared=distance*distance;
    smoothing_kdtree::parallel_for(n,200,threads,[&](size_t i) {
        vector<uint32_t> branches,singles;
        collect_collapsed(kdtree_cache_data->kdtree.root,points[i],squared,branches,singles);
        const uint32_t counts[2]={static_cast<uint32_t>(branches.size()),static_cast<uint32_t>(singles.size())};
        std::sort(branches.begin(),branches.end());
        std::sort(singles.begin(),singles.end());
        branches.insert(branches.end(),singles.begin(),singles.end());
        kdtree_cache_data->records[i]=encode_record(branches,counts,kdtree_cache_data->record_sizes[i]);
    });
    if (show_output_log)
        cout << "----- KD-tree cache generation " << std::chrono::duration<double>(
            std::chrono::steady_clock::now()-logging_begin).count() << " s" << endl;
    return kdtree_cache_data;
}

void validate_record(const uint8_t *p,size_t size,size_t n) {
    if (size<8) throw runtime_error("Truncated KD-tree cache record header.");
    uint32_t counts[2]; std::memcpy(counts,p,8);
    if (counts[0]>n || counts[1]>n-counts[0]) throw runtime_error("Invalid KD-tree cache record counts.");
    size_t pos=8;
    for (unsigned c=0;c<2;++c) {
        uint32_t used=0, previous=0;
        bool first=true;
        while (used<counts[c]) {
            if (size-pos<5) throw runtime_error("Truncated KD-tree cache delta block.");
            const unsigned extra=p[pos++];
            uint32_t id=0;
            for (unsigned b=0;b<4;++b) id|=uint32_t(p[pos++])<<(8*b);
            if (extra>=counts[c]-used || extra>size-pos || id>=n || (!first && id<=previous))
                throw runtime_error("Invalid KD-tree cache delta block.");
            first=false;
            for (unsigned k=0;k<extra;++k) {
                const unsigned delta=p[pos++];
                if (!delta || uint64_t(id)+delta>=n) throw runtime_error("Invalid KD-tree cache 8-bit delta.");
                id+=delta;
            }
            previous=id; used+=extra+1;
        }
    }
    if (pos!=size) throw runtime_error("Trailing KD-tree cache record bytes.");
}

string make_cache_filename(const string &kdtree_cache_data_folder,double radius) {
    radius=cap_smoothing_radius(radius);
    if (radius==0) radius=0; // Normalize negative zero.
    static_assert(sizeof(double)==8,"Cache requires 64-bit double.");
    std::ostringstream cache_filename_stream;
    cache_filename_stream.imbue(std::locale::classic());
    // Decimal metres with leading zeros to at least eight characters.
    // Retain fractional precision instead of rounding distinct radii together.
    cache_filename_stream.precision(std::numeric_limits<double>::max_digits10);
    cache_filename_stream<<kdtree_cache_data_folder<<"kdtree_cache_data_r_"
        <<std::setfill('0')<<std::setw(8)<<radius<<"_m.bin";
    return cache_filename_stream.str();
}

void encode_record_reusable(const vector<uint32_t> &ids, const uint32_t *r,
                            vector<uint8_t> &record) {
    size_t size = 0;
    // First pass sizes; second writes. Keep the output vector's capacity across batches.
    for (int pass = 0; pass < 2; ++pass) {
        size_t pos = 8;
        uint8_t *out = pass ? record.data() : nullptr;
        if (out) std::memcpy(out, r, 8);
        for (size_t channel = 0, begin = 0; channel < 2; ++channel) {
            const size_t end = begin + r[channel];
            for (size_t j = begin; j < end;) {
                const uint32_t base = ids[j];
                size_t stop = j + 1;
                // Count byte stores length minus one: 0..255 means 1..256 references.
                while (stop < end && stop - j < 256 && ids[stop] - ids[stop - 1] <= 255) ++stop;
                if (out) {
                    out[pos] = static_cast<uint8_t>(stop - j - 1);
                    for (unsigned b = 0; b < 4; ++b) out[pos + 1 + b] = static_cast<uint8_t>(base >> (8 * b));
                    for (size_t k = j + 1; k < stop; ++k) out[pos + 5 + k - j - 1] = static_cast<uint8_t>(ids[k] - ids[k - 1]);
                }
                pos += 5 + stop - j - 1;
                j = stop;
            }
            begin = end;
        }
        if (!pass) { size = pos; record.resize(pos); }
        else if (pos != size) throw runtime_error("Spatial record size mismatch.");
    }
    // Validate all decoded indices and both channel boundaries, including duplicates.
    size_t pos = 8, j = 0;
    for (size_t channel = 0; channel < 2; ++channel) {
        const size_t end = j + r[channel];
        while (j < end) {
            if (size - pos < 5) throw runtime_error("Truncated spatial block.");
            const unsigned extra = record[pos++];
            uint32_t base = 0;
            for (unsigned b = 0; b < 4; ++b) base |= uint32_t(record[pos++]) << (8 * b);
            if (extra >= end - j || extra > size - pos || ids[j++] != base)
                throw runtime_error("Spatial block validation failed.");
            uint32_t previous = base;
            for (unsigned k = 0; k < extra; ++k) {
                const uint32_t decoded = previous + record[pos++];
                if (ids[j++] != decoded) throw runtime_error("Spatial decoded index mismatch.");
                previous = decoded;
            }
        }

    }
    if (pos != size) throw runtime_error("Trailing spatial bytes.");
}

unique_ptr<KdTreeCacheData> prepare_cache_geometry(
    const double *lat,const double *lon,size_t n,double radius,int threads,
    vector<smoothing_kdtree::Point> &points) {
    smoothing_kdtree::validate_number_of_threads(threads);
    radius=cap_smoothing_radius(radius);
    points=generate_vector_of_kdtree_points_from_lat_lon_points_provided_as_arrays(lat,lon,n,threads);
    unique_ptr<KdTreeCacheData> kdtree_cache_data(new KdTreeCacheData());
    kdtree_cache_data->radius=radius;
    vector<pair<uint64_t,uint32_t>> keys(n);
    smoothing_kdtree::parallel_for(n,200,threads,[&](size_t i) {
        keys[i]=std::make_pair(spherical_s2_style_key(points[i].coords[0],points[i].coords[1],points[i].coords[2]),
                               static_cast<uint32_t>(i));
    });
    sort_keys(keys,threads);
    kdtree_cache_data->spatial_to_original.resize(n);
    {
        vector<smoothing_kdtree::Point> reordered(n);
        smoothing_kdtree::parallel_for(n,200,threads,[&](size_t i) {
            const uint32_t old=keys[i].second;
            kdtree_cache_data->spatial_to_original[i]=old;
            reordered[i]=points[old];
            reordered[i].index=static_cast<smoothing_kdtree::Index>(i);
        });
        points.swap(reordered);
    }
    vector<pair<uint64_t,uint32_t>>().swap(keys);
    if (n && !kdtree_cache_data->kdtree.buildKdTree_and_do_not_change_the_kdtree_points_vector(points))
        throw runtime_error("Cannot build KD-tree cache tree.");
    return kdtree_cache_data;
}

struct generation_worker_buffers {
    vector<uint32_t> branches, singles;
};

void collect_sorted_record(const smoothing_kdtree::Node *root,
    const smoothing_kdtree::Point &point, float squared,
    generation_worker_buffers &scratch, uint32_t (&counts)[2]) {
    scratch.branches.clear();
    scratch.singles.clear();
    collect_collapsed(root,point,squared,scratch.branches,scratch.singles);
    counts[0]=static_cast<uint32_t>(scratch.branches.size());
    counts[1]=static_cast<uint32_t>(scratch.singles.size());
    std::sort(scratch.branches.begin(),scratch.branches.end());
    std::sort(scratch.singles.begin(),scratch.singles.end());
    scratch.branches.insert(scratch.branches.end(),scratch.singles.begin(),scratch.singles.end());
}

// Write one radius's KD-tree cache data using shared geometry and reusable batch slots.
// A persistent worker team keeps search buffers across every batch of this radius.
void write_cache_data_from_prepared_geometry_batched(
    const KdTreeCacheData &kdtree_cache_data,
    const vector<smoothing_kdtree::Point> &points,double radius,const string &kdtree_cache_data_folder,int threads,
    const vector<uint8_t> &tree,vector<vector<uint8_t>> &records,vector<uint8_t> &write_buffer,
    const bool show_output_log = true) {
    radius=cap_smoothing_radius(radius);
    const size_t n=points.size();
    if (show_output_log)
        cout << "----- Generating kd-tree cache (radius " << radius << " m)" << endl;
    smoothing_io::File file(nullptr,&fclose);
    const char magic[8]={'S','2','C','A','C','H','E','1'};
    const uint32_t endian=0x01020304, point_count=static_cast<uint32_t>(n);
    const uint64_t tree_size=tree.size();
    file=smoothing_io::open_cache(make_cache_filename(kdtree_cache_data_folder,radius),"wb",show_output_log);
    smoothing_io::write(file.get(),magic,1,8);
    smoothing_io::write(file.get(),&endian,4,1);
    smoothing_io::write(file.get(),&point_count,4,1);
    smoothing_io::write(file.get(),&radius,8,1);
    smoothing_io::write(file.get(),&tree_size,8,1);
    smoothing_io::write(file.get(),tree.data(),1,tree.size());
    smoothing_io::write(file.get(),kdtree_cache_data.spatial_to_original.data(),4,n);

    unsigned next_progress_percent=10;
    if (show_output_log) cout << "---0 %" << endl;
    const size_t batch_size=8192;
    const float distance=static_cast<float>(great_circle_distance_to_euclidian_distance(radius));
    const float squared=distance*distance;
    std::exception_ptr failure;
    bool failed_batch=false;
    #pragma omp parallel num_threads(threads)
    {
        generation_worker_buffers scratch;
        for (size_t begin=0;begin<n;begin+=batch_size) {
            const size_t count=std::min(batch_size,n-begin);
            #pragma omp barrier
            #pragma omp for schedule(static)
            for (size_t offset=0;offset<count;++offset) {
                try {
                    uint32_t counts[2];
                    collect_sorted_record(kdtree_cache_data.kdtree.root,points[begin+offset],squared,scratch,counts);
                    encode_record_reusable(scratch.branches,counts,records[offset]);
                } catch (...) {
                    #pragma omp critical (saved_s2_generation_failure)
                    { if (!failure) failure=std::current_exception(); }
                }
            }
            // Both implicit barriers are required: writing sees completed records,
            // and no worker reuses a slot until writing the previous batch finishes.
            #pragma omp single
            {
                if (!failure) {
                    try {
                        size_t bytes=0;
                        for (size_t offset=0;offset<count;++offset) {
                            const size_t size=records[offset].size();
                            if (bytes>numeric_limits<size_t>::max()-sizeof(uint64_t) ||
                                size>numeric_limits<size_t>::max()-sizeof(uint64_t)-bytes)
                                throw length_error("Saved KD-tree write batch is too large.");
                            bytes+=sizeof(uint64_t)+size;
                        }
                        write_buffer.resize(bytes);
                        size_t position=0;
                        for (size_t offset=0;offset<count;++offset) {
                            const uint64_t size=records[offset].size();
                            std::memcpy(write_buffer.data()+position,&size,sizeof(size));
                            position+=sizeof(size);
                            std::memcpy(write_buffer.data()+position,records[offset].data(),
                                        records[offset].size());
                            position+=records[offset].size();
                        }
                        smoothing_io::write(file.get(),write_buffer.data(),1,write_buffer.size());
                        if (show_output_log) {
                            while (next_progress_percent < 100 &&
                                   uint64_t(begin+count)*100 >= uint64_t(next_progress_percent)*n) {
                                cout << "---" << next_progress_percent << " %" << endl;
                                next_progress_percent+=10;
                            }
                        }
                    } catch (...) { failure=std::current_exception(); }
                }
                failed_batch=static_cast<bool>(failure);
            }
            if (failed_batch) break;
        }
    }
    if (failure) std::rethrow_exception(failure);
    smoothing_io::close(file);
    if (show_output_log) cout << "---100 %" << endl;

}

// Geometry and serialized tree bytes are independent of smoothing radius.
// Process radii sequentially so encoded memory remains one bounded batch,
// rather than scaling with the number of complete caches or requested radii.
void generate_cache_data_and_write_batched(
    const double *lat,const double *lon,size_t n,const vector<double> &radii,
    const string &kdtree_cache_data_folder,int threads, const bool show_output_log = true) {
    const auto logging_begin = std::chrono::steady_clock::now();
    smoothing_kdtree::validate_number_of_threads(threads);
    if (radii.empty()) throw invalid_argument("At least one smoothing radius is required.");
    vector<double> capped_radii;
    capped_radii.reserve(radii.size());
    vector<string> filenames;
    filenames.reserve(radii.size());
    for (double radius:radii) {
        capped_radii.push_back(cap_smoothing_radius(radius));
        filenames.push_back(make_cache_filename(kdtree_cache_data_folder,capped_radii.back()));
    }
    std::sort(filenames.begin(),filenames.end());
    if (std::adjacent_find(filenames.begin(),filenames.end())!=filenames.end())
        throw invalid_argument("Smoothing radii must have distinct cache filenames after capping.");

    vector<smoothing_kdtree::Point> points;
    auto kdtree_cache_data=prepare_cache_geometry(lat,lon,n,capped_radii.front(),threads,points);
    const auto tree=kdtree_cache_data->kdtree.serialize();
    const size_t capacity=std::min<size_t>(8192,n);
    vector<vector<uint8_t>> records(capacity);
    vector<uint8_t> write_buffer;
    for (size_t i=0;i<radii.size();++i)
        write_cache_data_from_prepared_geometry_batched(*kdtree_cache_data,points,capped_radii[i],kdtree_cache_data_folder,threads,tree,records,write_buffer,show_output_log);
    if (show_output_log)
        cout << "----- KD-tree cache generation " << std::chrono::duration<double>(
            std::chrono::steady_clock::now()-logging_begin).count() << " s" << endl;
}

void generate_cache_data_and_write_batched(
    const double *lat,const double *lon,size_t n,double radius,const string &kdtree_cache_data_folder,int threads, const bool show_output_log = true) {
    generate_cache_data_and_write_batched(lat,lon,n,vector<double>(1,radius),kdtree_cache_data_folder,threads, show_output_log);
}

unique_ptr<KdTreeCacheData> read_cache_data_from_binary_file(const string &cache_filename, const bool show_output_log = true) {
    const auto logging_begin = std::chrono::steady_clock::now();
    auto file=smoothing_io::open_cache(cache_filename,"rb",show_output_log);
    char magic[8]; uint32_t endian,n;
    smoothing_io::read(file.get(),magic,1,8);
    if (std::memcmp(magic,"S2CACHE1",8)) throw runtime_error("Not a version-1 collapsed S2 KD-tree cache; regenerate it.");
    smoothing_io::read(file.get(),&endian,4,1);
    if (endian!=0x01020304) throw runtime_error("Incompatible KD-tree cache byte order.");
    smoothing_io::read(file.get(),&n,4,1);
    smoothing_kdtree::validate_point_count(n);
    unique_ptr<KdTreeCacheData> kdtree_cache_data(new KdTreeCacheData());
    smoothing_io::read(file.get(),&kdtree_cache_data->radius,8,1);
    validate_cache_radius(kdtree_cache_data->radius);
    uint64_t tree_size;
    smoothing_io::read(file.get(),&tree_size,8,1);
    const uint64_t expected_tree_size=40+uint64_t(n)*(3*sizeof(kdtree::PointType)+sizeof(kdtree::IndexType)+16);
    if (tree_size!=expected_tree_size || tree_size>std::numeric_limits<size_t>::max())
        throw runtime_error("Invalid KD-tree cache tree size.");
    {
        vector<uint8_t> tree(static_cast<size_t>(tree_size));
        smoothing_io::read(file.get(),tree.data(),1,tree.size());
        kdtree_cache_data->kdtree.deserialize(tree);
    }
    vector<uint8_t> seen(n,0);
    vector<const smoothing_kdtree::Node*> pending;
    if (kdtree_cache_data->kdtree.root) pending.push_back(kdtree_cache_data->kdtree.root);
    size_t visited=0;
    while (!pending.empty()) {
        const auto *node=pending.back(); pending.pop_back();
        const size_t id=node->val.index;
        if (id>=n || seen[id]) throw runtime_error("Invalid KD-tree cache tree indices.");
        seen[id]=1; ++visited;
        if (node->left_child()) pending.push_back(node->left_child());
        if (node->right_child()) pending.push_back(node->right_child());
    }
    if (visited!=n) throw runtime_error("KD-tree cache tree count mismatch.");
    std::fill(seen.begin(),seen.end(),0);
    kdtree_cache_data->spatial_to_original.resize(n);
    smoothing_io::read(file.get(),kdtree_cache_data->spatial_to_original.data(),4,n);
    for (uint32_t id:kdtree_cache_data->spatial_to_original) {
        if (id>=n || seen[id]) throw runtime_error("Invalid KD-tree cache S2 point permutation.");
        seen[id]=1;
    }
    kdtree_cache_data->records.resize(n); kdtree_cache_data->record_sizes.resize(n);
    for (size_t i=0;i<n;++i) {
        uint64_t size; smoothing_io::read(file.get(),&size,8,1);
        if (size<8 || size>8+uint64_t(n)*5 || size>std::numeric_limits<size_t>::max())
            throw runtime_error("Invalid KD-tree cache record length.");
        kdtree_cache_data->record_sizes[i]=static_cast<size_t>(size);
        kdtree_cache_data->records[i].reset(new uint8_t[kdtree_cache_data->record_sizes[i]]);
        smoothing_io::read(file.get(),kdtree_cache_data->records[i].get(),1,kdtree_cache_data->record_sizes[i]);
        validate_record(kdtree_cache_data->records[i].get(),kdtree_cache_data->record_sizes[i],n);
    }
    smoothing_io::require_end(file.get()); smoothing_io::close(file);
    if (show_output_log)
        cout << "----- KD-tree cache loading " << std::chrono::duration<double>(
            std::chrono::steady_clock::now()-logging_begin).count() << " s" << endl;
    return kdtree_cache_data;
}
} // namespace saved_s2

void generate_kdtree_cache_data_and_write_to_disk(
    const double *lat,const double *lon,size_t n,double radius,const string kdtree_cache_data_folder,int threads, const bool show_output_log = true) {
    saved_s2::generate_cache_data_and_write_batched(lat,lon,n,radius,kdtree_cache_data_folder,threads, show_output_log);
}

void free_kdtree_cache_data(char *&kdtree_cache_data_pointer) {
    delete reinterpret_cast<KdTreeCacheData*>(kdtree_cache_data_pointer);
    kdtree_cache_data_pointer=nullptr;
}

void smooth_field_using_kdtree_cache_data(
    const double *area,const double *field,size_t n,const char *kdtree_cache_data_pointer,double *output,int threads, const bool show_output_log = true) {
    const auto logging_begin = std::chrono::steady_clock::now();
    smoothing_kdtree::validate_number_of_threads(threads);
    smoothing_kdtree::validate_point_count(n);
    smoothing_kdtree::require_buffer(kdtree_cache_data_pointer,1);
    smoothing_kdtree::require_buffer(area,n);
    smoothing_kdtree::require_buffer(field,n);
    smoothing_kdtree::require_buffer(output,n);
    const auto &kdtree_cache_data=*reinterpret_cast<const KdTreeCacheData*>(kdtree_cache_data_pointer);
    if (kdtree_cache_data.records.size()!=n) throw invalid_argument("KD-tree cache point count mismatch.");
    unique_ptr<double[]> raw_products,raw_weights,raw_sums,raw_areas;
    double *products,*weights,*sums,*areas;
    // Preparation writes every product/weight. Each subtree writes its
    // sums before its parent reads them; decoding follows the full pass.
    raw_products.reset(new double[n]); raw_weights.reset(new double[n]);
    raw_sums.reset(new double[n]); raw_areas.reset(new double[n]);
    products=raw_products.get(); weights=raw_weights.get();
    sums=raw_sums.get(); areas=raw_areas.get();
    const auto prepare_point=[&](const size_t j) {
        const uint32_t original=kdtree_cache_data.spatial_to_original[j];
        weights[j]=area[original]; products[j]=field[original]*weights[j];
    };
    #pragma omp parallel for num_threads(threads) schedule(static)
    for (size_t j=0;j<n;++j) prepare_point(j);
    saved_s2::parallel_sums(kdtree_cache_data.kdtree.root,n,products,weights,sums,areas,threads);
    const auto decode_record=[&](const uint32_t i) {
        const uint8_t *p=kdtree_cache_data.records[i].get();
        uint32_t counts[2]; std::memcpy(counts,p,8); p+=8;
        double numerator=0,denominator=0;
        for (unsigned channel=0;channel<2;++channel) {
            const double *v=channel ? products : sums;
            const double *a=channel ? weights : areas;
            uint32_t consumed=0;
            while (consumed<counts[channel]) {
                const unsigned extra=*p++;
                uint32_t id=uint32_t(p[0]) | (uint32_t(p[1])<<8) | (uint32_t(p[2])<<16) | (uint32_t(p[3])<<24);
                p+=4; numerator+=v[id]; denominator+=a[id];
                for (unsigned k=0;k<extra;++k) { id+=*p++; numerator+=v[id]; denominator+=a[id]; }
                consumed+=extra+1;
            }
        }
        output[kdtree_cache_data.spatial_to_original[i]]=weights[i]>0 && denominator>0 ? numerator/denominator : 0;
    };
    #pragma omp parallel for num_threads(threads) schedule(static)
    for (uint32_t i=0;i<n;++i) decode_record(i);
    if (show_output_log)
        cout << "----- Smoothing " << std::chrono::duration<double>(
            std::chrono::steady_clock::now()-logging_begin).count() << " s" << endl;
}

void smooth_multiple_fields_simultaneously_using_kdtree_cache_data(
    const double *area_size,const double *fields,size_t number_of_points,size_t number_of_fields,
    const KdTreeCacheData *kdtree_cache_data_pointer,double *smoothed_fields,int number_of_threads, const bool show_output_log = true) {
    const auto logging_begin = std::chrono::steady_clock::now();
    smoothing_kdtree::validate_number_of_threads(number_of_threads);
    smoothing_kdtree::validate_point_count(number_of_points);
    if (!number_of_fields) throw invalid_argument("At least one field is required.");
    smoothing_kdtree::require_buffer(kdtree_cache_data_pointer,1);
    const auto &kdtree_cache_data=*kdtree_cache_data_pointer;
    if (kdtree_cache_data.records.size()!=number_of_points ||
        kdtree_cache_data.spatial_to_original.size()!=number_of_points)
        throw invalid_argument("KD-tree cache point count mismatch.");
    if (!number_of_points) return;
    if (number_of_fields>static_cast<size_t>(numeric_limits<ptrdiff_t>::max())/
        sizeof(double)/number_of_points)
        throw length_error("Multiple-field array dimensions are too large.");
    const size_t value_count=number_of_points*number_of_fields;
    smoothing_kdtree::require_buffer(area_size,number_of_points);
    smoothing_kdtree::require_buffer(fields,value_count);
    smoothing_kdtree::require_buffer(smoothed_fields,value_count);

    // A fixed set of uninitialized buffers is reused for every field. Shared
    // areas, subtree area sums and final denominators are calculated only once.
    unique_ptr<double[]> weights(new double[number_of_points]);
    unique_ptr<double[]> products(new double[number_of_points]);
    unique_ptr<double[]> sums(new double[number_of_points]);
    unique_ptr<double[]> areas(new double[number_of_points]);
    unique_ptr<double[]> denominators(new double[number_of_points]);
    const auto plan=saved_s2::prepare_batch_subtree_sum_plan(
        kdtree_cache_data.kdtree.root,number_of_points,number_of_threads);
    #pragma omp parallel for num_threads(number_of_threads) schedule(static)
    for (size_t i=0;i<number_of_points;++i)
        weights[i]=area_size[kdtree_cache_data.spatial_to_original[i]];

    for (size_t field_index=0;field_index<number_of_fields;++field_index) {
        const double *input=fields+field_index*number_of_points;
        double *output=smoothed_fields+field_index*number_of_points;
        #pragma omp parallel for num_threads(number_of_threads) schedule(static)
        for (size_t i=0;i<number_of_points;++i)
            products[i]=input[kdtree_cache_data.spatial_to_original[i]]*weights[i];
        if (field_index==0)
            saved_s2::batch_subtree_sums(plan,products.get(),weights.get(),sums.get(),areas.get(),number_of_threads);
        else
            saved_s2::batch_subtree_field_sums(plan,products.get(),sums.get(),number_of_threads);
        if (field_index==0) {
            #pragma omp parallel for num_threads(number_of_threads) schedule(static)
            for (uint32_t i=0;i<number_of_points;++i) {
                const uint8_t *p=kdtree_cache_data.records[i].get();
                uint32_t counts[2]; std::memcpy(counts,p,8); p+=8;
                double numerator=0,denominator=0;
                for (unsigned channel=0;channel<2;++channel) {
                    const double *values=channel ? products.get() : sums.get();
                    const double *area_values=channel ? weights.get() : areas.get();
                    uint32_t consumed=0;
                    while (consumed<counts[channel]) {
                        const unsigned extra=*p++;
                        uint32_t id=uint32_t(p[0]) | (uint32_t(p[1])<<8) | (uint32_t(p[2])<<16) | (uint32_t(p[3])<<24);
                        p+=4; numerator+=values[id]; denominator+=area_values[id];
                        for (unsigned k=0;k<extra;++k) { id+=*p++; numerator+=values[id]; denominator+=area_values[id]; }
                        consumed+=extra+1;
                    }
                }
                denominators[i]=denominator;
                output[kdtree_cache_data.spatial_to_original[i]]=weights[i]>0 && denominator>0 ? numerator/denominator : 0;
            }
        } else {
            #pragma omp parallel for num_threads(number_of_threads) schedule(static)
            for (uint32_t i=0;i<number_of_points;++i) {
                const uint8_t *p=kdtree_cache_data.records[i].get();
                uint32_t counts[2]; std::memcpy(counts,p,8); p+=8;
                double numerator=0;
                for (unsigned channel=0;channel<2;++channel) {
                    const double *values=channel ? products.get() : sums.get();
                    uint32_t consumed=0;
                    while (consumed<counts[channel]) {
                        const unsigned extra=*p++;
                        uint32_t id=uint32_t(p[0]) | (uint32_t(p[1])<<8) | (uint32_t(p[2])<<16) | (uint32_t(p[3])<<24);
                        p+=4; numerator+=values[id];
                        for (unsigned k=0;k<extra;++k) { id+=*p++; numerator+=values[id]; }
                        consumed+=extra+1;
                    }
                }
                const double denominator=denominators[i];
                output[kdtree_cache_data.spatial_to_original[i]]=weights[i]>0 && denominator>0 ? numerator/denominator : 0;
            }
        }
    }
    if (show_output_log)
        cout << "----- Smoothing " << number_of_fields << " fields simultaneously " << std::chrono::duration<double>(
            std::chrono::steady_clock::now()-logging_begin).count() << " s" << endl;
}




