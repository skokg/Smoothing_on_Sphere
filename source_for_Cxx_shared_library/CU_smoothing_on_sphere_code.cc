#include "CU_kdtree_with_index.cc"

double Earth_radius= 6371*1000;

double great_circle_distance_to_euclidian_distance(const double great_circle_distance)
	{
	return(2*Earth_radius*sin(great_circle_distance/(2*Earth_radius)));
	}

double euclidian_distance_to_great_circle_distance(const double euclidian_distance)
	{
	return(2*Earth_radius*asin(euclidian_distance/(2*Earth_radius)));
	}

void generate_fxareasize_and_areasize_Bounding_Box_data(const kdtree::KdTreeNode * const node, const double * const f, const double * const area_size,  vector <double> &fxareasize_BB_data, vector <double> &area_size_BB_data)
	{
	size_t index = node->val.index;
	double fxareasize_BB = f[index]*area_size[index];
	double areasize_BB = area_size[index];
	if(node->leftNode != nullptr)
		{
		generate_fxareasize_and_areasize_Bounding_Box_data(node->leftNode, f, area_size, fxareasize_BB_data, area_size_BB_data);
		fxareasize_BB+= fxareasize_BB_data[node->leftNode->val.index];
		areasize_BB+= area_size_BB_data[node->leftNode->val.index];
		}
	if(node->rightNode != nullptr)
		{
		generate_fxareasize_and_areasize_Bounding_Box_data(node->rightNode, f, area_size, fxareasize_BB_data, area_size_BB_data);
		fxareasize_BB+= fxareasize_BB_data[node->rightNode->val.index];
		areasize_BB+= area_size_BB_data[node->rightNode->val.index];
		}

	fxareasize_BB_data[index] =fxareasize_BB;
	area_size_BB_data[index] =areasize_BB;
	}


void generate_fxareasize_and_areasize_Bounding_Box_data_multiple_fields_simultaneously(const kdtree::KdTreeNode * const node, const double * const * const f, const double * const * const area_size,  const size_t number_of_fields, vector <vector <double>> &fxareasize_BB_data, vector <vector <double>> &area_size_BB_data  )
	{
	size_t index = node->val.index;
	vector <double> fxareasize_BB (number_of_fields,0);
	vector <double> areasize_BB (number_of_fields,0);
	for (unsigned long iff=0; iff < number_of_fields; iff++)
		{
		fxareasize_BB[iff] = f[iff][index]*area_size[iff][index];
		areasize_BB[iff] = area_size[iff][index];
		}
	if(node->leftNode != nullptr)
		{
		generate_fxareasize_and_areasize_Bounding_Box_data_multiple_fields_simultaneously(node->leftNode, f, area_size, number_of_fields, fxareasize_BB_data, area_size_BB_data);
		for (unsigned long iff=0; iff < number_of_fields; iff++)
			{
			fxareasize_BB[iff] += fxareasize_BB_data[iff][node->leftNode->val.index];
			areasize_BB[iff] += area_size_BB_data[iff][node->leftNode->val.index];
			}
		}
	if(node->rightNode != nullptr)
		{
		generate_fxareasize_and_areasize_Bounding_Box_data_multiple_fields_simultaneously(node->rightNode, f, area_size, number_of_fields, fxareasize_BB_data, area_size_BB_data);
		for (unsigned long iff=0; iff < number_of_fields; iff++)
			{
			fxareasize_BB[iff] += fxareasize_BB_data[iff][node->rightNode->val.index];
			areasize_BB[iff] += area_size_BB_data[iff][node->rightNode->val.index];
			}
		}

	for (unsigned long iff=0; iff < number_of_fields; iff++)
		{
		fxareasize_BB_data[iff][index] =fxareasize_BB[iff];
		area_size_BB_data[iff][index] =areasize_BB[iff];
		}
	}



void get_fxareasize_and_areasize_sums_of_points_in_radius( const kdtree::KdTreeNode * const node, const kdtree::KdTree &kdtree, const kdtree::Point_str& point, const float distance, const float distance_sqaured, const vector <double> &fxareasize, const double * const areasize, const vector <double> &fxareasize_BB_data, const vector <double> &area_size_BB_data, double &fxareasize_sum, double &areasize_sum)
	{

	if(node == nullptr) return;

	if (kdtree.test_if_the_hypersphere_and_hyperrectangle_intersect_using_sqr_radius( point.coords, distance_sqaured, node->coords_min, node->coords_max ))
		{
		float dist_sqr = kdtree.calDist_sqr(point, node->val);
		bool inside_sphere = false;
		if(dist_sqr < distance_sqaured)
			inside_sphere = true;

		if (inside_sphere && kdtree.test_if_the_hyperrectangle_is_fully_inside_the_sphere_using_sqr_radius( point.coords, distance_sqaured, node->coords_min, node->coords_max))
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

			get_fxareasize_and_areasize_sums_of_points_in_radius(node->rightNode, kdtree, point, distance, distance_sqaured, fxareasize, areasize, fxareasize_BB_data, area_size_BB_data,fxareasize_sum, areasize_sum );
			get_fxareasize_and_areasize_sums_of_points_in_radius(node->leftNode, kdtree, point, distance, distance_sqaured, fxareasize, areasize, fxareasize_BB_data, area_size_BB_data,fxareasize_sum, areasize_sum );
			}
		}
	}

void get_fxareasize_and_areasize_sums_of_points_in_radius_multiple_fields_simultaneously( const kdtree::KdTreeNode * const node, const kdtree::KdTree &kdtree, const kdtree::Point_str& point, const float distance, const float distance_sqaured, const vector <vector <double>> &fxareasize, const double * const * const areasize, const vector <vector <double>> &fxareasize_BB_data, const vector <vector <double>> &area_size_BB_data, vector <double> &fxareasize_sum, vector <double> &areasize_sum, const size_t number_of_fields)
	{

	if(node == nullptr) return;

	if (kdtree.test_if_the_hypersphere_and_hyperrectangle_intersect_using_sqr_radius( point.coords, distance_sqaured, node->coords_min, node->coords_max ))
		{
		float dist_sqr = kdtree.calDist_sqr(point, node->val);
		bool inside_sphere = false;
		if(dist_sqr < distance_sqaured)
			inside_sphere = true;

		if (inside_sphere && kdtree.test_if_the_hyperrectangle_is_fully_inside_the_sphere_using_sqr_radius( point.coords, distance_sqaured, node->coords_min, node->coords_max))
			{
			for (unsigned long iff=0; iff < number_of_fields; iff++)
				{
				fxareasize_sum[iff]+=fxareasize_BB_data[iff][node->val.index];
				areasize_sum[iff]+=area_size_BB_data[iff][node->val.index];
				}
			}

		else
			{
			if (inside_sphere)
				{
				for (unsigned long iff=0; iff < number_of_fields; iff++)
					{
					fxareasize_sum[iff]+=fxareasize[iff][node->val.index];
					areasize_sum[iff]+=areasize[iff][node->val.index];
					}
				}

			get_fxareasize_and_areasize_sums_of_points_in_radius_multiple_fields_simultaneously(node->rightNode, kdtree, point, distance, distance_sqaured, fxareasize, areasize, fxareasize_BB_data, area_size_BB_data,fxareasize_sum, areasize_sum, number_of_fields );
			get_fxareasize_and_areasize_sums_of_points_in_radius_multiple_fields_simultaneously(node->leftNode, kdtree, point, distance, distance_sqaured, fxareasize, areasize, fxareasize_BB_data, area_size_BB_data,fxareasize_sum, areasize_sum, number_of_fields );
			}
		}
	}

vector <kdtree::Point_str> generate_vector_of_kdtree_points_from_lat_lon_points(const vector <double> &lat, const vector <double> &lon)
	{
	ERRORIF(lat.size() != lon.size());

	// create kdtree points
	vector <kdtree::Point_str> kdtree_points (lat.size(),{{(float)0,(float)0,(float)0},(size_t)-1});
	#pragma omp parallel for num_threads(NUMBER_OF_THREADS)
	for (size_t il=0; il < lon.size(); il++)
		{
		double x,y,z;
		spherical_to_cartesian_coordinates(deg2rad(lat[il]), deg2rad(lon[il]), Earth_radius, x, y, z);
		kdtree_points[il].coords[0] = (float)x;
		kdtree_points[il].coords[1] = (float)y;
		kdtree_points[il].coords[2] = (float)z;
		kdtree_points[il].index = il;
		}
	return(kdtree_points);
	}

vector <kdtree::Point_str> generate_vector_of_kdtree_points_from_lat_lon_points_provided_as_arrays(const double * const lat, const double * const lon, const size_t size)
	{
	// create kdtree points
	vector <kdtree::Point_str> kdtree_points (size,{{(float)0,(float)0,(float)0},(size_t)-1});
	#pragma omp parallel for num_threads(NUMBER_OF_THREADS)
	for (size_t il=0; il < size; il++)
		{
		double x,y,z;
		spherical_to_cartesian_coordinates(deg2rad(lat[il]), deg2rad(lon[il]), Earth_radius, x, y, z);
		kdtree_points[il].coords[0] = (float)x;
		kdtree_points[il].coords[1] = (float)y;
		kdtree_points[il].coords[2] = (float)z;
		kdtree_points[il].index = il;
		}
	return(kdtree_points);
	}

void smooth_field_using_kd_tree(const double r_kernel_in_metres, const kdtree::KdTree &kdtree, const double * const lat, const double * const lon, const double * const area_size, const double * const f, const size_t number_of_points, double * const f_smoothed)
	{

	auto begin3 = std::chrono::high_resolution_clock::now();
	vector <double> fxareasize (number_of_points,0);
	for (unsigned long il=0; il < number_of_points; il++)
		fxareasize[il] = f[il]*area_size[il];
	//cout << sum_vector(fxareasize) << endl;
	//cout << sum_vector(area_size) << endl;

	vector <double> fxareasize_BB_data (number_of_points,0);
	vector <double> area_size_BB_data (number_of_points,0);

	generate_fxareasize_and_areasize_Bounding_Box_data(kdtree.root, f, area_size,  fxareasize_BB_data, area_size_BB_data);

	/*cout << " ------------- " << endl;
	cout << kdtree.root->val.index << endl;
	cout << fxareasize_BB_data[kdtree.root->val.index] << endl;
	cout << area_size_BB_data[kdtree.root->val.index] << endl;
	cout << fxareasize_BB_data[kdtree.root->rightNode->val.index] << endl;
	cout << fxareasize_BB_data[kdtree.root->leftNode->val.index] << endl;
	*/

	float r_Tunel_Distance = great_circle_distance_to_euclidian_distance(r_kernel_in_metres);
	//cout << "r_kernel: " << r_kernel << endl;
	//cout << "r_Tunel_Distance: " << r_Tunel_Distance << endl;
	float r_Tunel_Distance_squared = r_Tunel_Distance * r_Tunel_Distance;


	//vector <double> f_smoothed (number_of_points,0);

	long couter_threshold = number_of_points/100;

	//#pragma omp parallel for num_threads(NUMBER_OF_THREADS)
	//#pragma omp parallel for num_threads(NUMBER_OF_THREADS) schedule(static, 10)
	#pragma omp parallel for num_threads(NUMBER_OF_THREADS) schedule(dynamic,200)
	for (unsigned long il=0; il < number_of_points; il++)
	//for (unsigned long il=0; il < lat.size(); il++)
		{
		if (il % couter_threshold ==0)
		 cout << "-- " << round_to_digits ( (double)il / (double)number_of_points * 100.0 , 1) << " % " << endl;

		double fxareasize_sum = 0;
		double areasize_sum = 0;

		double x,y,z;
		spherical_to_cartesian_coordinates(deg2rad(lat[il]), deg2rad(lon[il]), Earth_radius, x, y, z);

		get_fxareasize_and_areasize_sums_of_points_in_radius(kdtree.root, kdtree, {{(float)x, (float)y,(float)z},(size_t)il}, r_Tunel_Distance, r_Tunel_Distance_squared, fxareasize, area_size, fxareasize_BB_data, area_size_BB_data, fxareasize_sum, areasize_sum);

		if (area_size[il] > 0)
			f_smoothed[il] = fxareasize_sum/areasize_sum;
		else
			f_smoothed[il] = 0;

		//f_smoothed.set(ix,iy,sum/(double)node_list.size());
		//cout << il << " " << f_smoothed[il]<< endl;
		}
	cout << "----- Smoothing " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin3).count() * 1e-9 << " s" << endl;

	//return(f_smoothed);
	}

void smooth_field_using_kd_tree_multiple_fields_simultaneously(const double r_kernel_in_metres, const kdtree::KdTree &kdtree, const double * const lat, const double * const lon, const double * const * const area_size, const double * const * const f, const size_t number_of_points, const size_t number_of_fields, double * const * const f_smoothed)
	{

	auto begin3 = std::chrono::high_resolution_clock::now();
	vector <vector <double>> fxareasize (number_of_fields, vector <double> (number_of_points,0));
	//vector <double> fxareasize (number_of_points,0);
	for (unsigned long iff=0; iff < number_of_fields; iff++)
		for (unsigned long il=0; il < number_of_points; il++)
			fxareasize[iff][il] = f[iff][il]*area_size[iff][il];
	//cout << sum_vector(fxareasize) << endl;
	//cout << sum_vector(area_size) << endl;

	vector <vector <double>> fxareasize_BB_data (number_of_fields, vector <double> (number_of_points,0));
	vector <vector <double>> area_size_BB_data (number_of_fields, vector <double> (number_of_points,0));

	generate_fxareasize_and_areasize_Bounding_Box_data_multiple_fields_simultaneously(kdtree.root, f, area_size, number_of_fields, fxareasize_BB_data, area_size_BB_data);

	//cout << " ------------- " << endl;
	//cout << kdtree.root->val.index << endl;
	//cout << fxareasize_BB_data[kdtree.root->val.index] << endl;
	//cout << area_size_BB_data[kdtree.root->val.index] << endl;
	//cout << fxareasize_BB_data[kdtree.root->rightNode->val.index] << endl;
	//cout << fxareasize_BB_data[kdtree.root->leftNode->val.index] << endl;


	float r_Tunel_Distance = great_circle_distance_to_euclidian_distance(r_kernel_in_metres);
	//cout << "r_kernel: " << r_kernel << endl;
	//cout << "r_Tunel_Distance: " << r_Tunel_Distance << endl;
	float r_Tunel_Distance_squared = r_Tunel_Distance * r_Tunel_Distance;


	//vector <double> f_smoothed (number_of_points,0);

	long couter_threshold = number_of_points/100;

	//#pragma omp parallel for num_threads(NUMBER_OF_THREADS)
	//#pragma omp parallel for num_threads(NUMBER_OF_THREADS) schedule(static, 10)
	#pragma omp parallel for num_threads(NUMBER_OF_THREADS) schedule(dynamic,1000)
	for (unsigned long il=0; il < number_of_points; il++)
	//for (unsigned long il=0; il < lat.size(); il++)
		{
		if (il % couter_threshold ==0)
		 cout << "-- " << round_to_digits ( (double)il / (double)number_of_points * 100.0 , 1) << " % " << endl;

		vector <double> fxareasize_sum (number_of_fields,0);
		vector <double> areasize_sum (number_of_fields,0);

		double x,y,z;
		spherical_to_cartesian_coordinates(deg2rad(lat[il]), deg2rad(lon[il]), Earth_radius, x, y, z);

		get_fxareasize_and_areasize_sums_of_points_in_radius_multiple_fields_simultaneously(kdtree.root, kdtree, {{(float)x, (float)y,(float)z},(size_t)il}, r_Tunel_Distance, r_Tunel_Distance_squared, fxareasize, area_size, fxareasize_BB_data, area_size_BB_data, fxareasize_sum, areasize_sum, number_of_fields);

		for (unsigned long iff=0; iff < number_of_fields; iff++)
			{
			if (area_size[iff][il] > 0)
				f_smoothed[iff][il] = fxareasize_sum[iff]/areasize_sum[iff];
			else
				f_smoothed[iff][il] = 0;
			}

		//f_smoothed.set(ix,iy,sum/(double)node_list.size());
		//cout << il << " " << f_smoothed[il]<< endl;
		}
	cout << "----- Smoothing " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin3).count() * 1e-9 << " s" << endl;

	//return(f_smoothed);

	}

void Read_smoothing_data_from_binary_file(const string &fname, uint32_t **&data_pointer, const uint32_t &number_of_points)
	{
	cout << "Reading file: " << fname << endl;
	ERRORIF(data_pointer != nullptr);

	data_pointer = new uint32_t*[number_of_points];

	FILE * pFile;
	pFile = fopen ( fname.c_str() , "rb" );
	if (pFile==NULL)
		error(AT,FU, "Problems opening file: " + fname);

	vector <uint32_t> temp_buffer (number_of_points,0);

	uint32_t current_smoothing_point;
	uint32_t previous_smoothing_point;
	uint32_t added_smoothing_points_size;
	uint32_t removed_smoothing_points_size;

	uint32_t counter = 0;
	size_t result;
	while (counter < number_of_points)
		{

		result = fread ((char*) &current_smoothing_point,1,sizeof(uint32_t),pFile);
		if (result != sizeof(uint32_t))
			error(AT,FU, "Problems reading file: " + fname);

		//cout << current_smoothing_point << endl;

		result = fread ((char*) &previous_smoothing_point,1,sizeof(uint32_t),pFile);
		if (result != sizeof(uint32_t))
			error(AT,FU, "Problems reading file: " + fname);

		//cout << previous_smoothing_point << endl;

		result = fread ((char*) &added_smoothing_points_size,1,sizeof(uint32_t),pFile);
		if (result != sizeof(uint32_t))
			error(AT,FU, "Problems reading file: " + fname);

		//cout << added_smoothing_points_size << endl;

		result = fread ((char*) &removed_smoothing_points_size,1,sizeof(uint32_t),pFile);
		if (result != sizeof(uint32_t))
			error(AT,FU, "Problems reading file: " + fname);

		//cout << removed_smoothing_points_size << endl;

		data_pointer[counter] = new uint32_t[added_smoothing_points_size + removed_smoothing_points_size + 4];

		data_pointer[counter][0] = current_smoothing_point;
		data_pointer[counter][1] = previous_smoothing_point;
		data_pointer[counter][2] = added_smoothing_points_size;
		data_pointer[counter][3] = removed_smoothing_points_size;

		result = fread ((char*)  &data_pointer[counter][4],sizeof(uint32_t),added_smoothing_points_size,pFile);
		if (result != added_smoothing_points_size)
			error(AT,FU, "Problems reading file: " + fname);

		result = fread ((char*)  &data_pointer[counter][4 + added_smoothing_points_size],sizeof(uint32_t),removed_smoothing_points_size,pFile);
		if (result != removed_smoothing_points_size)
			error(AT,FU, "Problems reading file: " + fname);

		counter++;

		}

	if (counter != number_of_points)
		error(AT,FU, "Problems reading file: " + fname);

	fclose (pFile);
	}


void free_smoothing_data_memory(uint32_t **&data_pointer, const uint32_t number_of_points)
	{
	for (uint32_t il = 0; il < number_of_points; il++)
		 delete [] data_pointer[il];

	delete [] data_pointer;

	data_pointer=nullptr;
	}

void calculate_sqr_distance_to_all_points_in_a_vector(const float &x_p, const float &y_p, const float &z_p, const vector <float> &x, const vector <float> &y, const vector <float> &z, vector <float> &sqr_distance)
	{
	#pragma omp parallel for num_threads(NUMBER_OF_THREADS)
	for (unsigned long il=0; il < x.size(); il++)
		{
		float dx = x_p - x[il];
		float dy = y_p - y[il];
		float dz = z_p - z[il];
		sqr_distance[il] = dx*dx + dy*dy + dz*dz;
		}
	}


void generate_smoothing_data_for_the_overlap_detection_and_write_it_to_disk(const double * const lat, const double * const lon, const size_t number_of_points, const vector <double> &smoothing_kernel_radius_in_metres, const string output_folder, const size_t starting_point)
	{

	cout << "----- Generating smoothing sequence for the points " << endl;

	// ------------------------------------------------
	// Generate a smoothing sequence for the points
	// ------------------------------------------------
	// create kdtree points
	vector <kdtree::Point_str> kdtree_points = generate_vector_of_kdtree_points_from_lat_lon_points_provided_as_arrays(lat, lon, number_of_points);

	// initialization of list and other parmeters
	//vector <uint32_t> list_previous_smoothing_point (number_of_points,0);
	//vector <uint32_t> list_total_steps (number_of_points,0);
	//vector <double> list_distance (number_of_points,-1);

	vector <uint32_t> iterative_sequence_current_smoothing_point;
	vector <uint32_t> iterative_sequence_previous_smoothing_point;

	uint32_t current_point_index=starting_point;
	uint32_t rebalance_kdtree_threshold = 100;

	// create kdtree for unassigned points
	kdtree::KdTree kdtree_unassigned_points;
	//auto begin2 = std::chrono::high_resolution_clock::now();
	kdtree_unassigned_points.buildKdTree_and_do_not_change_the_kdtree_points_vector(kdtree_points);
	//cout << "----- kdtree construction " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin2).count() * 1e-9 << endl;

	// create empty kdtree for assigned points
	kdtree::KdTree kdtree_assigned_points;

	long output_stepX = number_of_points / 10;

	uint32_t number_of_unassigned_points=number_of_points;
	uint32_t rebalance_kdtree_threshold_count=0;
	while (number_of_unassigned_points > 0)
		{

		if (number_of_unassigned_points % output_stepX == 0)
			{
			cout << "---" << round_to_digits( (double)(number_of_points - number_of_unassigned_points)/(double)number_of_points * 100, 1) << " %" << endl;
			//cout << "----- time " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin4).count() * 1e-9 << "s" << endl;
			}

		long ind1;
		long ind2;

		// special case for the first point
		if (number_of_unassigned_points == number_of_points)
			{
			ind1 = 	current_point_index;
			ind2 = 	current_point_index;
			}

		else
			{
			// the closest unassigned point to the current_point_index - for the first point this is going to be the same point
			auto node1 = kdtree_unassigned_points.findNearestNode(kdtree_points[current_point_index]);
			//temp.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() * 1e-9);
			ERRORIF(node1 == nullptr);
			ind1=node1->val.index;

			auto node2 = kdtree_assigned_points.findNearestNode(kdtree_points[ind1]);
			//temp.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() * 1e-9);
			ERRORIF(node2 == nullptr);
			ind2=node2->val.index;
			}

		// move ind1 point from the one kdtree to another
		kdtree_unassigned_points.deleteNode(kdtree_points[ind1]);
		kdtree_assigned_points.addNode(kdtree_points[ind1]);

		// add point data to sequences and lists
		iterative_sequence_current_smoothing_point.push_back(ind1);
		iterative_sequence_previous_smoothing_point.push_back(ind2);

		//list_previous_smoothing_point[ind1] = ind2;
		//list_total_steps[ind1] = list_total_steps[ind2] + 1;
		//list_distance[ind1] = euclidian_distance_to_great_circle_distance(sqrt(squared_euclidian_distance_multidimensional(kdtree_points[ind1].coords, kdtree_points[ind2].coords)));


		number_of_unassigned_points--;

		rebalance_kdtree_threshold_count++;
		if (rebalance_kdtree_threshold_count == rebalance_kdtree_threshold)
			{
			//auto begin3 = std::chrono::high_resolution_clock::now();
			kdtree_assigned_points.rebuild_a_balanced_KdTree();
			//cout << "----- kdtree rebalancing " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin3).count() * 1e-9 << endl;
			//cout << rebalance_kdtree_threshold << " " << number_of_points - number_of_unassigned_points << " " << number_of_unassigned_points << " " << current_point_index << " " << ind1 << " " << ind2 << endl;
			rebalance_kdtree_threshold_count=0;
			rebalance_kdtree_threshold*=2;
			}
		current_point_index = ind1;
		}

	// free some memory
	free_vector_memory(kdtree_points);
	kdtree_unassigned_points.free_memory();
	kdtree_assigned_points.free_memory();

	//Write_vector_to_netcdf_file_compressed_multiple_types_supported(output_folder + "out_iterative_sequence_current_smoothing_point.nc",iterative_sequence_current_smoothing_point,1);
	//Write_vector_to_netcdf_file_compressed_multiple_types_supported(output_folder + "out_iterative_sequence_previous_smoothing_point.nc",iterative_sequence_previous_smoothing_point,1);
	//Write_vector_to_netcdf_file_compressed_multiple_types_supported(output_folder + "out_list_previous_smoothing_point.nc",list_previous_smoothing_point,1);
	//Write_vector_to_netcdf_file_compressed_multiple_types_supported(output_folder + "out_list_total_steps.nc",list_total_steps,1);
	//Write_vector_to_netcdf_file_compressed_multiple_types_supported(output_folder + "out_list_distance.nc",list_distance,1);

	// -------------------------------------------------------------------------------------------
	// Generate the smoothing information based on the previusly determined smoothing sequence
	// -------------------------------------------------------------------------------------------


	// read iterative_sequence
	//vector <uint32_t> iterative_sequence_current_smoothing_point;
	//Read_vector_from_netcdf_file_multiple_types_supported(output_folder + "out_iterative_sequence_current_smoothing_point.nc", "var", iterative_sequence_current_smoothing_point);
	//vector <uint32_t> iterative_sequence_previous_smoothing_point;
	//Read_vector_from_netcdf_file_multiple_types_supported(output_folder + "out_iterative_sequence_previous_smoothing_point.nc", "var", iterative_sequence_previous_smoothing_point);*/

	cout << "----- Overlap detection and writing the smoothing data to the disk " << endl;

	vector <float> squared_kernel_tunnel_distance_radius_in_metres;
	for (unsigned long il=0; il < smoothing_kernel_radius_in_metres.size(); il++)
		squared_kernel_tunnel_distance_radius_in_metres.push_back( pow ( great_circle_distance_to_euclidian_distance( (double)smoothing_kernel_radius_in_metres[il] ),2));

	// create x,y,z vectors
	vector <float> x (number_of_points,0);
	vector <float> y (number_of_points,0);
	vector <float> z (number_of_points,0);
	for (unsigned long il=0; il < number_of_points; il++)
		{
		double x_temp,y_temp,z_temp;
		spherical_to_cartesian_coordinates(deg2rad(lat[il]), deg2rad(lon[il]), Earth_radius, x_temp, y_temp, z_temp);
		x[il] = (float)x_temp;
		y[il] = (float)y_temp;
		z[il] = (float)z_temp;
		}

	vector < vector <uint32_t> > smoothing_information (smoothing_kernel_radius_in_metres.size(), vector <uint32_t> (0,0));

	vector <FILE *> pFile;

	for (unsigned long ir=0; ir < smoothing_kernel_radius_in_metres.size(); ir++)
		{
		string fname = output_folder + "out_smoothing_information_r_"+output_leading_zero_string(round(smoothing_kernel_radius_in_metres[ir]),8)+"_m.bin";
		cout << "Opening file for write: " << fname << endl;
		pFile.push_back(fopen(fname.c_str(),"wb"));  // w for write, b for binary
		}


	vector <float> sqr_distance_current (number_of_points,0);
	vector <float> sqr_distance_previous (number_of_points,0);

	uint32_t current_smoothing_point_in_previous_step = 0;

	vector <vector <uint32_t>> added_smoothing_points  (smoothing_kernel_radius_in_metres.size(), vector <uint32_t> (0,0));
	vector <vector <uint32_t>> removed_smoothing_points  (smoothing_kernel_radius_in_metres.size(), vector <uint32_t> (0,0));

	//auto begin4 = std::chrono::high_resolution_clock::now();
	long output_step = iterative_sequence_current_smoothing_point.size() / 1000;
	for (unsigned long il=0; il < iterative_sequence_current_smoothing_point.size(); il++)
		{
		if (il % output_step == 0)
			{
			cout << "---" << round_to_digits( (double)il/(double)iterative_sequence_current_smoothing_point.size() * 100, 1) << " %" << endl;
			//cout << "----- time " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin4).count() * 1e-9 << "s" << endl;
			}

		uint32_t current_smoothing_point = iterative_sequence_current_smoothing_point[il];
		uint32_t previous_smoothing_point = iterative_sequence_previous_smoothing_point[il];

		// reset vectors
		for (uint32_t ir=0; ir < smoothing_kernel_radius_in_metres.size(); ir++)
			{
			added_smoothing_points[ir].clear();
			removed_smoothing_points[ir].clear();
			}

		// spacial case for the first point
		if (il == 0)
			{
			calculate_sqr_distance_to_all_points_in_a_vector(x[current_smoothing_point], y[current_smoothing_point], z[current_smoothing_point], x, y, z, sqr_distance_current);
			for (uint32_t ir=0; ir < smoothing_kernel_radius_in_metres.size(); ir++)
				for (uint32_t ip=0; ip < number_of_points; ip++)
					if (sqr_distance_current[ip] < squared_kernel_tunnel_distance_radius_in_metres[ir])
						added_smoothing_points[ir].push_back(ip);
			}

		else
			{
			calculate_sqr_distance_to_all_points_in_a_vector(x[current_smoothing_point], y[current_smoothing_point], z[current_smoothing_point], x, y, z, sqr_distance_current);

			// only calculate distances if previous_current_smoothing_point != previous_smoothing_point
			if (current_smoothing_point_in_previous_step != previous_smoothing_point)
				calculate_sqr_distance_to_all_points_in_a_vector(x[previous_smoothing_point], y[previous_smoothing_point], z[previous_smoothing_point], x, y, z, sqr_distance_previous);

			//auto begin3 = std::chrono::high_resolution_clock::now();
			for (uint32_t ir=0; ir < smoothing_kernel_radius_in_metres.size(); ir++)
				#pragma omp parallel for num_threads(NUMBER_OF_THREADS)
				for (uint32_t ip=0; ip < number_of_points; ip++)
					{
					bool inside_current = false;
					if (sqr_distance_current[ip] < squared_kernel_tunnel_distance_radius_in_metres[ir])
						inside_current = true;

					bool inside_previous = false;
					if (sqr_distance_previous[ip] < squared_kernel_tunnel_distance_radius_in_metres[ir])
						inside_previous = true;

					if (inside_current && !inside_previous)
						{
						#pragma omp critical (added)
						{added_smoothing_points[ir].push_back(ip);}
						}

					if (!inside_current && inside_previous)
						{
						#pragma omp critical (removed)
						{removed_smoothing_points[ir].push_back(ip);}
						}
					}
			//cout << "----- time " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin3).count() * 1e-9 << endl;
			}

		for (unsigned long ir=0; ir < smoothing_kernel_radius_in_metres.size(); ir++)
			{
			fwrite((char*)&current_smoothing_point,1,sizeof(uint32_t),pFile[ir]);
			fwrite((char*)&previous_smoothing_point,1,sizeof(uint32_t),pFile[ir]);

			uint32_t added_smoothing_points_size = added_smoothing_points[ir].size();
			uint32_t removed_smoothing_points_size = removed_smoothing_points[ir].size();
			fwrite((char*)&added_smoothing_points_size,1,sizeof(uint32_t),pFile[ir]);
			fwrite((char*)&removed_smoothing_points_size,1,sizeof(uint32_t),pFile[ir]);

			fwrite((char*)&added_smoothing_points[ir][0],added_smoothing_points_size,sizeof(uint32_t),pFile[ir]);
			fwrite((char*)&removed_smoothing_points[ir][0],removed_smoothing_points_size,sizeof(uint32_t),pFile[ir]);
			}

		sqr_distance_previous.swap(sqr_distance_current);
		current_smoothing_point_in_previous_step = current_smoothing_point;
		}

	for (unsigned long ir=0; ir < smoothing_kernel_radius_in_metres.size(); ir++)
		{
		fclose (pFile[ir]);
		}

	//cout << "----- time4 " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin4).count() * 1e-9 << endl;
	}



void smooth_field_using_overlap_detection(const double * const area_size, const double * const f, const size_t number_of_points, const uint32_t * const * const data_pointer,  double * const f_smoothed)
	{

	vector <double> f_x_area (number_of_points,0);
	#pragma omp parallel for num_threads(NUMBER_OF_THREADS)
	for (uint32_t il=0; il < number_of_points; il++)
		f_x_area[il] = f[il]*area_size[il];

	vector <double> f_x_area_sum (number_of_points,0);
	vector <double> area_sum (number_of_points,0);

	auto begin4 = std::chrono::high_resolution_clock::now();

	// first loop to precalculate the partial sums of all points - can be parallelized
	#pragma omp parallel for num_threads(NUMBER_OF_THREADS)
	for (uint32_t il = 0; il < number_of_points; il++)
		{
		double partial_f_x_area_sum=0;
		double partial_area_sum=0;

		const uint32_t * const data_pointer_for_point = data_pointer[il];

		// terms from added points
		uint32_t start_index = 4;
		uint32_t end_index = start_index + data_pointer_for_point[2];
		for (uint32_t ip = start_index; ip < end_index; ip++)
			{
			uint32_t index = data_pointer_for_point[ip];
			partial_f_x_area_sum+=f_x_area[index];
			partial_area_sum+=area_size[index];
			}

		// terms from removed points
		start_index = end_index ;
		end_index = start_index + data_pointer_for_point[3];
		for (uint32_t ip = start_index; ip < end_index; ip++)
			{
			uint32_t index = data_pointer_for_point[ip];
			partial_f_x_area_sum-=f_x_area[index];
			partial_area_sum-=area_size[index];
			}

		f_x_area_sum[data_pointer_for_point[0]] = partial_f_x_area_sum;
		area_sum[data_pointer_for_point[0]] = partial_area_sum;
		}

	// second loop to calculate the smoothed values - cannot be paralelized since it is iterative
	for (uint32_t il = 0; il < number_of_points; il++)
		{

		const uint32_t * const data_pointer_for_point = data_pointer[il];

		uint32_t current_smoothing_point = data_pointer_for_point[0];
		uint32_t previous_smoothing_point = data_pointer_for_point[1];

		// check if this point has a reference point
		if (current_smoothing_point != previous_smoothing_point)
			{
			f_x_area_sum[current_smoothing_point] += f_x_area_sum[previous_smoothing_point];
			area_sum[current_smoothing_point] += area_sum[previous_smoothing_point];
			}

		if (area_sum[current_smoothing_point] > 0)
			f_smoothed[current_smoothing_point] = f_x_area_sum[current_smoothing_point]/area_sum[current_smoothing_point];
		else
			f_smoothed[current_smoothing_point]=0;
		}

	// a third loop to set the missing data points to zero value - can be parallelized
	#pragma omp parallel for num_threads(NUMBER_OF_THREADS)
	for (uint32_t il = 0; il < number_of_points; il++)
		{
		if (area_size[il] == 0)
			f_smoothed[il]=0;
		}

	cout << "----- Smoothing " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin4).count() * 1e-9 << " s" << endl;
	}


void identify_branches_and_single_nodes_inside_the_radius( const kdtree::KdTreeNode * const node, const kdtree::KdTree &kdtree, const kdtree::Point_str& point, const float distance, const float distance_sqaured, vector <uint32_t> &branches_fully_inside_the_radius, vector <uint32_t> &single_nodes_inside_the_radius)
	{

	if(node == nullptr) return;

	if (kdtree.test_if_the_hypersphere_and_hyperrectangle_intersect_using_sqr_radius( point.coords, distance_sqaured, node->coords_min, node->coords_max ))
		{
		float dist_sqr = kdtree.calDist_sqr(point, node->val);
		bool inside_sphere = false;
		if(dist_sqr < distance_sqaured)
			inside_sphere = true;

		if (inside_sphere && kdtree.test_if_the_hyperrectangle_is_fully_inside_the_sphere_using_sqr_radius( point.coords, distance_sqaured, node->coords_min, node->coords_max))
			{
			branches_fully_inside_the_radius.push_back((uint32_t)node->val.index);
			}

		else
			{
			if (inside_sphere)
				{
				single_nodes_inside_the_radius.push_back((uint32_t)node->val.index);
				}

			identify_branches_and_single_nodes_inside_the_radius(node->rightNode, kdtree, point, distance, distance_sqaured, branches_fully_inside_the_radius, single_nodes_inside_the_radius );
			identify_branches_and_single_nodes_inside_the_radius(node->leftNode, kdtree, point, distance, distance_sqaured, branches_fully_inside_the_radius, single_nodes_inside_the_radius );
			}
		}
	}

void generate_smoothing_data_for_the_kdtree_based_approach_and_write_it_to_the_disk(const double * const lat, const double * const lon, const size_t number_of_points, const double smoothing_kernel_radius_in_metres, const string output_folder)
	{
	auto begin = std::chrono::high_resolution_clock::now();

	// construct kd-tree
	vector <kdtree::Point_str> kdtree_points = generate_vector_of_kdtree_points_from_lat_lon_points_provided_as_arrays(lat, lon, number_of_points);
	kdtree::KdTree kdtree;
	kdtree.buildKdTree_and_do_not_change_the_kdtree_points_vector(kdtree_points);

	cout << "----- Construction of kd-tree " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() * 1e-9 << " s" << endl;

	// calculate the smoothing data
	vector <vector <uint32_t>> branches_fully_inside_the_radius (number_of_points, vector <uint32_t> (0,0));
	vector <vector <uint32_t>> single_nodes_inside_the_radius (number_of_points, vector <uint32_t> (0,0));

	float tunnel_distance_radius_in_metres = (float) great_circle_distance_to_euclidian_distance(smoothing_kernel_radius_in_metres);
	float squared_tunnel_distance_radius_in_metres = tunnel_distance_radius_in_metres*tunnel_distance_radius_in_metres;

	begin = std::chrono::high_resolution_clock::now();

	#pragma omp parallel for num_threads(NUMBER_OF_THREADS)
	for (uint32_t il = 0; il < number_of_points; il++)
		{
		identify_branches_and_single_nodes_inside_the_radius( kdtree.root, kdtree, kdtree_points[il], tunnel_distance_radius_in_metres, squared_tunnel_distance_radius_in_metres, branches_fully_inside_the_radius[il], single_nodes_inside_the_radius[il]);
		}
	cout << "----- Generating smoothing data " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() * 1e-9 << " s" << endl;

	//for (uint32_t il = 0; il < number_of_points; il++)
	//	cout << il << " " << branches_fully_inside_the_radius[il].size() << " " << single_nodes_inside_the_radius[il].size() << endl;

	begin = std::chrono::high_resolution_clock::now();

	string fname = output_folder + "out_kdtree_smoothing_information_r_"+output_leading_zero_string(round(smoothing_kernel_radius_in_metres),8)+"_m.bin";

	// kd-tree data
	char *data_stream_pointer = nullptr;
	uint64_t size = 0;
	kdtree.write_KdTree_data_into_a_binary_data_stream(data_stream_pointer, size);

	//cout << "bb " << size << endl;

	FILE * pFile;
	pFile = fopen ( fname.c_str() , "wb" );
		if (pFile==NULL)
	error(AT,FU, "Problems opening file: " + fname);

	// write kd-tree data size
	size_t result = fwrite(&size,sizeof(uint64_t),1,pFile);
	if (result != 1)
		error(AT,FU, "Problems writing kdtree to file: " + fname);

	// write kd-tree data
	result = fwrite(data_stream_pointer,sizeof(char),size,pFile);
	if (result != size)
		error(AT,FU, "Problems writing kdtree to file: " + fname);

	delete [] data_stream_pointer;

	// write number of points
	uint32_t count = number_of_points;
	result = fwrite(&count,sizeof(uint32_t),1,pFile);
	if (result != 1)
		error(AT,FU, "Problems writing to file: " + fname);

	// write smoothing data
	for (uint32_t il = 0; il < number_of_points; il++)
		{
		uint32_t count1 = branches_fully_inside_the_radius[il].size();
		result = fwrite(&count1,sizeof(uint32_t),1,pFile);
		if (result != 1)
			error(AT,FU, "Problems writing to file: " + fname);

		uint32_t count2 = single_nodes_inside_the_radius[il].size();
		result = fwrite(&count2,sizeof(uint32_t),1,pFile);
		if (result != 1)
			error(AT,FU, "Problems writing to file: " + fname);

		result = fwrite(&branches_fully_inside_the_radius[il][0],sizeof(uint32_t),count1,pFile);
		if (result != count1)
			error(AT,FU, "Problems writing to file: " + fname);

		result = fwrite(&single_nodes_inside_the_radius[il][0],sizeof(uint32_t),count2,pFile);
		if (result != count2)
			error(AT,FU, "Problems writing to file: " + fname);
		}

	fclose (pFile);

	cout << "----- Writing smoothing data to disk " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() * 1e-9 << " s" << endl;

	}


struct smoothing_data_for_the_kdtree_based_approach
	{
	kdtree::KdTree kdtree;
	vector <uint32_t *> smoothing_data_for_point;

	void free_memory()
		{
		kdtree.free_memory();
		for (uint32_t il=0; il<smoothing_data_for_point.size(); il++)
			{
			if (smoothing_data_for_point[il] != nullptr)
				delete[] smoothing_data_for_point[il];
			}
		}

	};


void free_smoothing_data_memory_kdtree(char *&data_pointer)
	{
	smoothing_data_for_the_kdtree_based_approach * const smoothing_data = (smoothing_data_for_the_kdtree_based_approach * const) data_pointer;
	smoothing_data->free_memory();
	data_pointer=nullptr;
	}



void Read_smoothing_data_for_the_kdtree_based_approach_from_binary_file(const string &fname, char *&data_pointer)
	{
	cout << "Reading file: " << fname << endl;
	ERRORIF(data_pointer != nullptr);

	FILE * pFile;
	pFile = fopen ( fname.c_str() , "rb" );
	if (pFile==NULL)
		error(AT,FU, "Problems opening file: " + fname);


	// read kd-tree data size
	uint64_t kd_tree_data_size;
	size_t result = fread (&kd_tree_data_size,sizeof(uint64_t),1,pFile);
	if (result != 1)
		error(AT,FU, "Problems reading smoothig data from file: " + fname);

	//cout << kd_tree_data_size << endl;

	char * kd_tree_data = new char[kd_tree_data_size];

	// read kd-tree data
	result = fread (kd_tree_data,sizeof(char),kd_tree_data_size,pFile);
	if (result != kd_tree_data_size)
		error(AT,FU, "Problems reading smoothig data from file: " + fname);

	smoothing_data_for_the_kdtree_based_approach * smoothing_data = new smoothing_data_for_the_kdtree_based_approach();
	// reconstruct kd-tree
	smoothing_data->kdtree.reconstruct_KdTree_data_from_a_binary_data_stream(kd_tree_data, kd_tree_data_size);
	delete[] kd_tree_data;

	//cout << smoothing_data->kdtree.count_number_of_nonnull_nodes() << endl;

	// read number of points
	uint32_t number_of_points;
	result = fread (&number_of_points,sizeof(uint32_t),1,pFile);
	if (result != 1)
		error(AT,FU, "Problems reading smoothig data from file: " + fname);


	// read smoothing data
	for (uint32_t il = 0; il < number_of_points; il++)
		{
		uint32_t count1;
		result = fread (&count1,sizeof(uint32_t),1,pFile);
		if (result != 1)
			error(AT,FU, "Problems reading smoothig data from file: " + fname);

		uint32_t count2;
		result = fread (&count2,sizeof(uint32_t),1,pFile);
		if (result != 1)
			error(AT,FU, "Problems reading smoothig data from file: " + fname);

		uint32_t * temp = new uint32_t[2+count1+count2];
		temp[0] = count1;
		temp[1] = count2;

		result = fread (&temp[2],sizeof(uint32_t),count1+count2,pFile);
		if (result != count1+count2)
			error(AT,FU, "Problems reading smoothig data from file: " + fname);

		smoothing_data->smoothing_data_for_point.push_back(temp);
		}

	fclose (pFile);

	data_pointer = (char *) smoothing_data;
	}

void smooth_field_using_smoothing_data_for_the_kdtree_approach(const double * const area_size, const double * const f, const size_t number_of_points, const char * const data_pointer,  double * const f_smoothed)
	{
	auto begin = std::chrono::high_resolution_clock::now();

	const smoothing_data_for_the_kdtree_based_approach * const smoothing_data = (const smoothing_data_for_the_kdtree_based_approach * const) data_pointer;

	ERRORIF(number_of_points != smoothing_data->smoothing_data_for_point.size());

	vector <double> fxareasize (number_of_points,0);
	for (unsigned long il=0; il < number_of_points; il++)
		fxareasize[il] = f[il]*area_size[il];
	//cout << sum_vector(fxareasize) << endl;
	//cout << sum_vector(area_size) << endl;

	vector <double> fxareasize_BB_data (number_of_points,0);
	vector <double> area_size_BB_data (number_of_points,0);

	generate_fxareasize_and_areasize_Bounding_Box_data(smoothing_data->kdtree.root, f, area_size,  fxareasize_BB_data, area_size_BB_data);

	#pragma omp parallel for num_threads(NUMBER_OF_THREADS)
	for (uint32_t il=0; il < number_of_points; il++)
		{
		double fxareasize_temp=0;
		double area_size_temp=0;

		const uint32_t * const smoothing_data_for_single_point =  smoothing_data->smoothing_data_for_point[il];
		uint32_t start_index = 2;
		uint32_t end_index = start_index + smoothing_data_for_single_point[0];
		for (uint32_t ip=start_index; ip < end_index; ip++)
			{
			fxareasize_temp+=fxareasize_BB_data[smoothing_data_for_single_point[ip]];
			area_size_temp+=area_size_BB_data[smoothing_data_for_single_point[ip]];
			}
		start_index = end_index;
		end_index = start_index + smoothing_data_for_single_point[1];
		for (uint32_t ip=start_index; ip < end_index; ip++)
			{
			fxareasize_temp+=fxareasize[smoothing_data_for_single_point[ip]];
			area_size_temp+=area_size[smoothing_data_for_single_point[ip]];
			}

		if (area_size[il] > 0)
			f_smoothed[il] = fxareasize_temp / area_size_temp;
		else
			f_smoothed[il] = 0;

		}

	cout << "----- Smoothing " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin).count() * 1e-9 << " s" << endl;

	}




