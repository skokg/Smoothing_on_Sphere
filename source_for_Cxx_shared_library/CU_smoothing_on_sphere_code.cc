#include "CU_kdtree_with_index.cc"

double Earth_radius= 6371*1000;

double great_circle_distance_to_euclidian_distance(double great_circle_distance)
	{
	return(2*Earth_radius*sin(great_circle_distance/(2*Earth_radius)));
	}

double euclidian_distance_to_great_circle_distance(double euclidian_distance)
	{
	return(2*Earth_radius*asin(euclidian_distance/(2*Earth_radius)));
	}

void generate_fxareasize_and_areasize_Bounding_Box_data(const std::shared_ptr<kdtree::KdTreeNode>& node, const double *&f, const double *&area_size,  vector <double> &fxareasize_BB_data, vector <double> &area_size_BB_data)
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


void get_fxareasize_and_areasize_sums_of_points_in_radius( const std::shared_ptr<kdtree::KdTreeNode>& node, const kdtree::KdTree &kdtree, const kdtree::Point_str& point, const float &distance, const float &distance_sqaured, const vector <double> &fxareasize, const double *&areasize, const vector <double> &fxareasize_BB_data, const vector <double> &area_size_BB_data, double &fxareasize_sum, double &areasize_sum)
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

vector <kdtree::Point_str> generate_vector_of_kdtree_points_from_lat_lon_points_provided_as_arrays(const double *&lat, const double *&lon, size_t size)
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

void smooth_field_using_kd_tree(double r_kernel_in_metres, const kdtree::KdTree &kdtree, const double *lat, const double *lon, const double *area_size, const double *f, const size_t number_of_points, double *f_smoothed)
	{

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

	auto begin3 = std::chrono::high_resolution_clock::now();
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
	cout << "----- smoothing " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin3).count() * 1e-9 << " s" << endl;

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


void generate_smoothing_data_for_the_overlap_detection_and_write_it_to_disk(const double *lat, const double *lon, size_t number_of_points, vector <double> smoothing_kernel_radius_in_metres, string output_folder, size_t starting_point)
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



void smooth_field_using_overlap_detection(const double *area_size, const double *f, const size_t number_of_points, uint32_t **data_pointer,  double *f_smoothed)
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

		uint32_t *data_pointer_for_point = data_pointer[il];

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

		uint32_t *data_pointer_for_point = data_pointer[il];

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

	cout << "----- smoothing " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - begin4).count() * 1e-9 << " s" << endl;
	}





