
import numpy as np
import os
import ctypes

# search for the PAD C++ shared library file (PAD_on_sphere_Cxx_shared_library.so) in the same folder
libc = ctypes.CDLL(os.path.abspath(os.path.expanduser(os.path.dirname(__file__)))+ os.path.sep + "smoothing_on_sphere_Cxx_shared_library.so") 

# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

ND_POINTER_1D = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C")

# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

libc.free_mem_double_array.argtypes = [ctypes.POINTER(ctypes.c_double)]
libc.free_mem_double_array.restype = None

def check_array(array, array_name):
	
	if type(array) is not np.ndarray:
		print("ERROR: the "+array_name+" is not a numpy array of type numpy.ndarray - this is not permitted! Perhaps it is a masked arrays, which is also not permitted. Returning \"None\" as result!")
		return(False)
	
	# check dimensions of fields
	if array.ndim != 1:
		print("ERROR: the "+array_name+" array nedd to be one-dimensional. Returning \"None\" as result!")
		return(False)
	
	# compare the array has some elements
	if array.size == 0:
		print("ERROR: the dimension of the "+array_name+" arrays is zero. Returning \"None\" as result!")
		return(False)
	
	# detect non-numeric values
	result=np.where(np.isfinite(array) == False)
	if len(result[0]) > 0:
		print("ERROR: the "+array_name+" array contain some non-numeric values. Returning \"None\" as result!")
		return(False)
	
	# detect masked array
	if isinstance(array, np.ma.MaskedArray):
		print("ERROR: the "+array_name+" array is a masked array which is not allowed. Returning \"None\" as result!")
		return(False)	
	
	return(True)

# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------
libc.construct_KdTree_ctypes.argtypes = [ND_POINTER_1D, ND_POINTER_1D, ctypes.c_size_t]
libc.construct_KdTree_ctypes.restype = ctypes.c_void_p

def construct_KdTree(lat, lon):
	
	if check_array(lat, "lat") != True:
		return(None)
	
	if check_array(lon, "lon") != True:
		return(None)
	
	# compare lat lon dimensions
	if lat.shape != lon.shape:
		print("ERROR: the lat and lon arrays do not have the same shape. Returning \"None\" as result!")
		return(None)
	
	
	# convert array to np.float64 and continousarray if necessary - this format is required for the interaction with the C++ code
	lat = np.ascontiguousarray(lat, dtype = np.float64)
	lon = np.ascontiguousarray(lon, dtype = np.float64)
	
	kdtree_pointer = libc.construct_KdTree_ctypes(lat, lon, len(lat))
	
	return(kdtree_pointer)


libc.free_KdTree_memory_ctypes.argtypes = [ctypes.c_void_p]
libc.free_KdTree_memory_ctypes.restype = None

def free_KdTree_memory(kdtree_pointer):
	libc.free_KdTree_memory_ctypes(kdtree_pointer)

# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

libc.smooth_field_using_KdTree_ctypes.argtypes = [ctypes.c_double, ctypes.c_void_p, ND_POINTER_1D, ND_POINTER_1D, ND_POINTER_1D, ND_POINTER_1D, ctypes.c_size_t, ND_POINTER_1D]
#libc.smooth_field_using_KdTree_ctypes.restype = ctypes.POINTER(ctypes.c_double)
libc.smooth_field_using_KdTree_ctypes.restype = None

def smooth_field_using_KdTree(smoothing_kernel_radius_in_metres, kdtree_pointer, lat, lon, area_size, f):
	
	
	if check_array(lat, "lat") != True:
		return(None)
	
	if check_array(lon, "lon") != True:
		return(None)
	
	if check_array(area_size, "area_size") != True:
		return(None)
	
	if check_array(f, "f") != True:
		return(None)
	
	# compare lat lon dimensions
	if lat.shape != lon.shape or lat.shape != area_size.shape or lat.shape != f.shape:
		print("ERROR: the lat, lon, area_size and f arrays do not have the same shape. Returning \"None\" as result!")
		return(None)
	
	# convert array to np.float64 and contiguousarray  - this format is required for the interation with the C++ code
	lat = np.ascontiguousarray(lat, dtype = np.float64)
	lon = np.ascontiguousarray(lon, dtype = np.float64)
	area_size = np.ascontiguousarray(area_size, dtype = np.float64)
	f = np.ascontiguousarray(f, dtype = np.float64)
	
	if np.ndim(smoothing_kernel_radius_in_metres) != 0:
		print("ERROR: the smoothing_kernel_radius_in_metres should be a single value! Returning \"None\" as result!")
		return(None)
	
	# cast smoothing_kernel_radius_in_metres to np.float64 just in case - - this format is required for the interaction with the C++ code
	smoothing_kernel_radius_in_metres = np.float64(smoothing_kernel_radius_in_metres)
	
	# reserve memory for the outputed smoothed field
	f_smoothed = np.ascontiguousarray(np.zeros(len(f),dtype=np.float64))
	
	# calculate the smoothed values
	libc.smooth_field_using_KdTree_ctypes(smoothing_kernel_radius_in_metres, kdtree_pointer, lat, lon, area_size, f, len(f), f_smoothed)
	
	return(f_smoothed)

# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

libc.generate_smoothing_data_for_the_overlap_detection_and_write_it_to_disk_ctypes.argtypes = [ND_POINTER_1D, ND_POINTER_1D, ctypes.c_size_t, ND_POINTER_1D, ctypes.c_size_t, ctypes.c_char_p, ctypes.c_size_t]
#libc.smooth_field_using_KdTree_ctypes.restype = ctypes.POINTER(ctypes.c_double)
libc.generate_smoothing_data_for_the_overlap_detection_and_write_it_to_disk_ctypes.restype = None

def generate_smoothing_data_for_the_overlap_detection_based_approach_and_write_it_to_the_disk(lat, lon, smoothing_kernel_radius_in_metres, smoothing_data_folder):
	
	if check_array(lat, "lat") != True:
		return(None)
	
	if check_array(lon, "lon") != True:
		return(None)
	
	# compare lat lon dimensions
	if lat.shape != lon.shape:
		print("ERROR: the lat, lon, area_size and f arrays do not have the same shape. Returning \"None\" as result!")
		return(None)
	
	# convert array to np.float64 and contiguousarray  - this format is required for the interation with the C++ code
	lat = np.ascontiguousarray(lat, dtype = np.float64)
	lon = np.ascontiguousarray(lon, dtype = np.float64)
	
	
	if np.ndim(smoothing_kernel_radius_in_metres) > 1:
		print("ERROR: the smoothing_kernel_radius_in_metres needs to bo a single value or a 1D array. Returning \"None\" as result!")
		return(None)
	
	# if smoothing_kernel_radius_in_metres is a single value convert it to an array 
	if np.ndim(smoothing_kernel_radius_in_metres) == 0:
		smoothing_kernel_radius_in_metres = np.ascontiguousarray([smoothing_kernel_radius_in_metres], dtype = np.float64)
	
	# convert array to np.float64 and contiguousarray  - this format is required for the interaction with the C++ code
	smoothing_kernel_radius_in_metres = np.ascontiguousarray(smoothing_kernel_radius_in_metres, dtype = np.float64)
	
	# calculate the smoothed values
	libc.generate_smoothing_data_for_the_overlap_detection_and_write_it_to_disk_ctypes(lat,lon,len(lat),smoothing_kernel_radius_in_metres, len(smoothing_kernel_radius_in_metres), smoothing_data_folder, 0)	


# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

libc.read_smoothing_data_from_binary_file_ctypes.argtypes = [ctypes.c_char_p, ctypes.c_double, ctypes.c_size_t]
libc.read_smoothing_data_from_binary_file_ctypes.restype = ctypes.c_void_p


def read_smoothing_data_from_binary_file(smoothing_data_folder, smoothing_kernel_radius_in_metres, number_of_points):
	
	if np.ndim(smoothing_kernel_radius_in_metres) != 0:
		print("ERROR: the smoothing_kernel_radius_in_metres needs to bo a single value in the read_smoothing_data_from_binary_file function. Returning \"None\" as result!")
		return(None)
	
	smoothing_data_pointer = libc.read_smoothing_data_from_binary_file_ctypes(smoothing_data_folder, smoothing_kernel_radius_in_metres, number_of_points)
	return(smoothing_data_pointer)


# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

libc.free_smoothing_data_memory_ctypes.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
libc.free_smoothing_data_memory_ctypes.restype = None

def free_smoothing_data_memory(smoothing_data_pointer, number_of_points):
	libc.free_smoothing_data_memory_ctypes(smoothing_data_pointer, number_of_points)
	smoothing_data_pointer = None


# -----------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------

libc.smooth_field_using_overlap_detection_ctypes.argtypes = [ND_POINTER_1D, ND_POINTER_1D, ctypes.c_size_t, ctypes.c_void_p, ND_POINTER_1D]
libc.smooth_field_using_overlap_detection_ctypes.restype = None

def smooth_field_using_overlap_detection(area_size, f, smoothing_data_pointer):
	
	if check_array(area_size, "area_size") != True:
		return(None)
	
	if check_array(f, "f") != True:
		return(None)
	
	# compare dimensions
	if area_size.shape != f.shape:
		print("ERROR: the area_size and f arrays do not have the same shape. Returning \"None\" as result!")
		return(None)
	
	# convert array to np.float64 and contiguousarray  - this format is required for the interation with the C++ code
	area_size = np.ascontiguousarray(area_size, dtype = np.float64)
	f = np.ascontiguousarray(f, dtype = np.float64)
	
	# reserve memory for the outputed smoothed field
	f_smoothed = np.ascontiguousarray(np.zeros(len(f),dtype=np.float64))
	
	# calculate the smoothed values
	libc.smooth_field_using_overlap_detection_ctypes(area_size, f, len(f), smoothing_data_pointer, f_smoothed)
	
	return(f_smoothed)


