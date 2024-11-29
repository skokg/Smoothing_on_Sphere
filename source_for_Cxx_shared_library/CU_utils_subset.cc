
// error function
void error(const char *location,const char *function, string txt)
	{
	cout << "ERROR: " << location << " " << function << ": " <<  txt << endl;
	exit(1);
	}


// this foreces the vector to free its memory - only using .clear() does not free the memory !
template<typename T>
void free_vector_memory(std::vector<T>& vec)
	{
	std::vector<T>().swap(vec);
	}


template<typename T>
string output_vector_as_string(const std::vector<T>& vec, string separator)
	{
	ostringstream s1;
	for (long il=0; il < (long)vec.size(); il++)
		{
		s1 << vec[il];
		if (il< (long)vec.size() - 1)
			s1 << separator;
		}
	return(s1.str());
	}


double deg2rad(double kot)
	{
	return (kot*M_PI/180.0);
	}

double rad2deg(double kot)
	{
	return (180.0*kot/M_PI);
	}

void spherical_to_cartesian_coordinates(double lat, double lon, double r, double &x, double &y, double &z)
	{
	x=r*cos(lon) * cos(lat);
	y=r*sin(lon) * cos(lat);
	z=r*sin(lat);
	}


// for long number il it adds 0 in front so the number in mest chaacters long
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

// rounds a double to digits
double round_to_digits(double x, int digits)
	{
	return(round(pow(10,(double)digits)*x)/pow(10,(double)digits));
	}
