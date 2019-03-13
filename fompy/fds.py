
# -*- coding: utf-8 -*-
"""
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
fds.py 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This module includes the routines used to import data into a FoMpy Dataset. Some useful
examples showing how to import the data using these functions can be seen below:

Example
-------

If the parser is not defined, the default parser is used, where a simple voltage-current file
is assumed to be defined as input::

	import fompy
	path_file = './data/default/'
	fds = fompy.dataset(path_file)

or::

	import fompy
	path_file = './data/default/'
	fds = fompy.dataset(path_file,'data*',  parser=fompy.default)

Moreover, another types of parsers can be selected::

	import fompy
	path_file_JCJB = './data/sim_FinFET_vd_high/'
	fds = fompy.dataset(path_file_JCJB, parser=fompy.JCJB)

	path_file_mc = './data/mc_data/'
	fds = fompy.dataset(path_file_mc, parser=fompy.MC)


If the user wishes to remove several IV curves included in the parent folder, two options called
'interval' and 'exclude' can be passed, so the indexes defined are removed from the FoMpy Dataset::

	import fompy
	path_file_JCJB = './data/sim_FinFET_vd_high/'
	fds = fompy.dataset(path_file_JCJB, parser=fompy.JCJB, exclude=[5,6])
	fds = fompy.dataset(path_file_JCJB, parser=fompy.JCJB, interval=[0,8])

"""

from abc import ABCMeta, abstractmethod
import os, glob
from enum import Enum
import numpy as np
from scipy.stats import mode

#--------------------------------------------------------------------------------------------------------------
class _FompyDatasetStrategy(metaclass=ABCMeta):
	"""
	Abstract class containing the most important attributes of a semiconductor's IV curve.

	"""
	def __init__(self, **kwargs):
		self.dataset = []
		self.n_sims = 0
		self.sanity_array = []
		self.norm = 1
		self.ext_method = 'SD'
		self.drain_bias_label = None
		self.drain_bias_value = None
		self.interpolation = 'cubic_spline'
		self.filter = None


	def print_parameters(self):
		"""
		Prints the details of a FompyDataset.
		"""		
		print("\nNumber of simulations: ", self.n_sims)
		print("Normalization value: ", self.norm)
		print("Extraction method: ", self.ext_method)
		print("Drain Bias label: ", self.drain_bias_label)
		print("Drain Bias value: ", self.drain_bias_value)
		print("Interpolation: ", self.interpolation)
		print("Filtering: ", self.filter, "\n")



#--------------------------------------------------------------------------------------------------------------
class FompyDataset(_FompyDatasetStrategy):
	"""
	Class containing the simulated IV curves and their parameters.

	...

	Attributes
	----------
	dataset : array: double[]
		List containing all IV curves
	n_sims : int
		Number of simulated IV curves
	sanity_array : array: double[]
		Array of ones by default. If a simulation has failed, its index is converted to zero.
	norm : double
		Normalization value applied to the IV curve
	ext_method : str
		FompyDataset default method used to extrad the figures of merit
	drain_bias_label : str
		Either high or low.
	drain_bias_value : double
		Drain voltage value used to simulate the IV curves.
	interpolation : str
		IV curve interpolation method. Can either be 'cubic_spline',
		'akima','pchip' or 'linear'.
	filter : str
		IV curve filter method. So far only the polar filter has been implemented ('polar_filter').

	Methods
	-------
	print_parameters
		Prints all the non-array attributes.

	"""
	def __init__(self,**kwargs):
		_FompyDatasetStrategy.__init__(self, **kwargs)
		for key in kwargs:
			setattr(self, key, kwargs[key])

#--------------------------------------------------------------------------------------------------------------
class dataDAO(metaclass=ABCMeta):
	"""
	Data Acces Object Interface used to import the simulated IV curves to a FompyDataset.
	"""	

	@abstractmethod
	def load(self):
		pass
	@abstractmethod
	def save(self):
		pass

class daoFile(dataDAO, FompyDataset):
	"""
	Data Acces Object used to import from a set of files the simulated IV curves to a FompyDataset.

	...

	Attributes
	----------
	parser : str
		Type of parser, user-defined, that imports the simulated IV curves in a specific format.

	"""	

	def __init__(self, parser = None):
		self.parser = parser

	def load(self, path, parser = None, interval = None, exclude = None, skiprows = None, comments = None):
		"""
		Methods
		-------
		load(path,filename_user = None, parser = None, interval = None, exclude = None, skiprows = None, comments = None)
			Class method that extracts :math:`V_{TH}` of a semiconductor's IV curve.

		Parameters
		----------
		path : array_like, shape (n,)
			1-d array containing values of the independent variable.
		parser : function
			Function that implements how the data is imported to a Fompy Dataset. The list of available functions includes:
			'JCJB' and 'MC'.
		interval : array_like
			List of two int values: start(index of the first simulation to load into the Fompy Dataset)
			and end(index of the last simulation to load into the FompyDataset)
		exclude : array_like
			Index values of simulations to exclude.
		skiprows : int
			Number of rows to skip at the begining of a file. 0 rows are skipped by default.
		comments : str
			All the lines starting with this character are considered comments.
			'#' is used by default.

		"""
		fds = FompyDataset()
		path_subdirs = []
		path_filenames = []
		for dirname, dirnames, filenames in os.walk(path):
			# print path to all subdirectories.
			for subdirname in dirnames:
				path_subdirs.append(os.path.join(dirname, subdirname))
			# print path to all filenames.
			for filename in filenames:
				path_filenames.append(os.path.join(dirname, filename))

		path_subdirs = sorted(path_subdirs)
		path_filenames = sorted(path_filenames)
		# print(path_filenames)
		# print(path_subdirs)

		for i in range(len(path_subdirs)):
			try:
				if(path_subdirs[i] in str(path_filenames)):
					fds.sanity_array.append(1)
				else:
					fds.sanity_array.append(-1)
			except IndexError:
				fds.sanity_array.append(-1)
		
		# print(fds.sanity_array)

		if(parser==None):
			print('No parsers have been defined')
		else:
			parser(fds, path, path_subdirs, path_filenames, interval, exclude)
		
		# print(fds.dataset)
		fds = exclude_indexes(fds,interval, exclude)

		return fds

	def save(self, fds, path):
		length = len(fds.dataset)
		f = open(str(path), "w+")
		for i in range(length):
			if(fds.dataset[i] is not np.nan):
				f.write("#Simulation id {0}, Vd = [{1}]\n".format(i+1, fds.drain_bias_value))
				f.write("#Vg[V]\tId[V]\n")
				arr = np.array2string(fds.dataset[i][:], separator="\t")
				f.write(arr+"\n")
			else:
				f.write("Simulation id {0}\n".format(i+1))
				f.write("FAILURE TO LOAD THE DATA\n")
		f.close()


#--------------------------------------------------------------------------------------------------------------
# def default(fds, path, path_subdirs, path_filenames, interval, exclude, skiprows, comments):
# 	"""Function that imports the simulated data from a given file and
# 	stores it into a FoMpy Dataset. The format of this file is asummed to have comments
# 	starting with '#', without a header and two columns, one for the voltage and one for the currents, separated
# 	by '\t' delimiters.

# 	Parameters
# 	----------
# 	fds : FoMpy Dataset
# 		Structure of data containing the most important parameters of a semiconductor's IV curve.
# 	path : str
# 		Parent path where the simulations are stored
# 	path_subdirs : str
# 		List of alphabetically sorted subdirectories found inside the parent directory.
# 	path_filenames : str
# 		List of alphabetically sorted files found inside the parent directory.
# 	interval : array_like
# 		List of two int values: start(index of the first simulation to load into the Fompy Dataset)
# 		and end(index of the last simulation to load into the FompyDataset)
# 	exclude : array_like
# 		Index values of simulations to exclude.
# 	skiprows : int
# 		Number of rows to skip at the begining of a file. 0 rows are skipped by default.
# 	comments : str
# 		All the lines starting with this character are considered comments.
# 		'#' is used by default.
		
# 	Returns
# 	-------
# 	fds : FoMpy Dataset
# 		Structure of data containing the FoMpy Dataset.

# 	"""
# 	extension = [".txt", ".dat", ".out", ".csv"]
	
# 	print(path_filenames)
# 	for item in path_filenames:
# 		# print(item)
# 		file_name = os.path.basename(item)
# 		for ext in extension:
# 			if(file_name.endswith(ext) is True):
# 				print('Reading data from %s file -> %s' % (ext ,file_name))
# 		# print(file_name)
# 		if(filename_user is not None):
# 			if(str(filename_user) in file_name):
# 				if(skiprows != None) & (comments != None):
# 					data_temp = np.loadtxt(item,comments=comments, skiprows=skiprows, unpack=True, delimiter='\t')
# 				else:
# 					data_temp = np.loadtxt(item)
# 				a, b = data_temp[0], data_temp[1]
# 				arr = np.column_stack((a, b))
# 				fds.dataset.append(arr)
# 		else:
# 			if(skiprows != None) & (comments != None):
# 				data_temp = np.loadtxt(item,comments=comments, skiprows=skiprows, unpack=True, delimiter='\t')
# 			else:
# 				data_temp = np.loadtxt(item)
# 			a, b = data_temp[0], data_temp[1]
# 			arr = np.column_stack((a, b))
# 			fds.dataset.append(arr)
			
# 		# print(fds.dataset)

# 	return fds

#--------------------------------------------------------------------------------------------------------------
def JCJB(fds, path, path_subdirs, path_filenames, interval, exclude):
	"""Function that imports the simulated data from a JCJB file and
	stores it into a FoMpy Dataset.

	Parameters
	----------
	fds : FoMpy Dataset
		Structure of data containing the most important parameters of a semiconductor's IV curve.
	path : str
		Parent path where the simulations are stored
	path_subdirs : str
		List of alphabetically sorted subdirectories found inside the parent directory.
	path_filenames : str
		List of alphabetically sorted files found inside the parent directory.
	interval : array_like
		List of two int values: start(index of the first simulation to load into the Fompy Dataset)
		and end(index of the last simulation to load into the FompyDataset)
	exclude : array_like
		Index values of simulations to exclude.

	Returns
	-------
	fds : FoMpy Dataset
		Structure of data containing the FoMpy Dataset.

	"""

	try:
		for path in path_filenames:
			if('JCJB.dat.0.4' in path):
				data_temp = np.loadtxt(path)
				a, b, c = data_temp[:,0], data_temp[:,1], data_temp[:,2]
				condition = np.mod(a, np.amax(a))==0
				vd_temp = np.copy(np.extract(condition, a))
				vg_temp = np.copy(np.extract(condition, b))
				id_temp = np.copy(np.extract(condition, c))
				Vd = np.copy(vd_temp[np.argmin(vg_temp):len(vg_temp)])
				Vg = np.copy(vg_temp[np.argmin(vg_temp):len(vg_temp)])
				Id = np.copy(id_temp[np.argmin(vg_temp):len(vg_temp)])
				arr = np.column_stack((Vd, Vg, Id))
				arr_f = np.column_stack((arr[:,1], arr[:,2]))
				vd = np.unique(Vd)
				fds.dataset.append(arr_f)
				fds.drain_bias_value = vd
	except (ValueError, IndexError):
		print('An error has ocurred during the JCJB data load')

	return fds
#--------------------------------------------------------------------------------------------------------------
def MC(fds, path, path_subdirs, path_filenames, interval, exclude):
	"""Function that imports the simulated data from a MC set of files and
	stores it into a FoMpy Dataset.

	Parameters
	----------
	fds : FoMpy Dataset
		Structure of data containing the most important parameters of a semiconductor's IV curve.
	path : str
		Parent path where the simulations are stored
	path_subdirs : str
		List of alphabetically sorted subdirectories found inside the parent directory.
	path_filenames : str
		List of alphabetically sorted files found inside the parent directory.
	interval : array_like
		List of two int values: start(index of the first simulation to load into the Fompy Dataset)
		and end(index of the last simulation to load into the FompyDataset)
	exclude : array_like
		Index values of simulations to exclude.
	skiprows : int
		Number of rows to skip at the begining of a file. 0 rows are skipped by default.
	comments : str
		All the lines starting with this character are considered comments.
		'#' is used by default.

	Returns
	-------
	fds : FoMpy Dataset
		Structure of data containing the FoMpy Dataset.

	"""

	index=0
	for i in path_subdirs:
		try:
			if('Graphs' in i):
				path_subdirs.pop(index)
		except IndexError:
			pass
		index = index + 1

	i = 0
	time_arr = []
	avg_curr_arr = []
	vg_arr = []
	if (interval is not None):
		ran = range(interval[0], interval[1])
	for path in path_filenames:
		if('fichero_particula' in path):
			try:
				time, avg_curr = np.genfromtxt(path, unpack=True, usecols = (0,11))
				time_arr.append(time)
				avg_curr_arr.append(np.abs(avg_curr))
			except ValueError:
				print('Non-valid data found in the file.')				
			i = i+1
		elif('voltajes' in path):
			vg, vd = np.genfromtxt(path, delimiter='\t',  skip_header = 1, usecols = (1,2), unpack=True)
			fds._drain_bias_value = vd
			vg_arr.append(vg)

	time_f = []
	avg_curr_f = []
	index_time = []
	for i in range(len(avg_curr_arr)):
		try:
			if(i==0):
				index_1 = np.argwhere(np.in1d(time_arr[i],np.intersect1d(time_arr[i], time_arr[i+1])) == True)
			index_2 = np.argwhere(np.in1d(time_arr[i],np.intersect1d(time_arr[i], time_arr[i+1])) == True)
			if(index_1[-1] <= index_2[-1]):
				index_time.append(index_1[-1])
			else:
				index_time.append(index_2[-1])
			min = np.min(index_time[-1])
		except IndexError:
			pass
	for i in range(len(avg_curr_arr)):
		try:
				time_f.append(time_arr[i][min])
				avg_curr_f.append(avg_curr_arr[i][min])										
		except IndexError:
			pass

	arr_f = np.column_stack((vg_arr, np.abs(avg_curr_f)))
	fds.dataset.append(arr_f)

	return fds

#--------------------------------------------------------------------------------------------------------------
def exclude_indexes(fds,interval, exclude):
	"""Function that generates an array of zeros and ones. If the simulation
	number 5 has failed, either because the folder is empty, or not enough voltages
	have been simulated, then fds.sanity_array[4] is set to zero.

	Parameters
	----------
	fds : FoMpy Dataset
		Structure of data containing the most important parameters of a semiconductor's IV curve.
	interval : array_like
		List of two int values: start(index of the first simulation to load into the Fompy Dataset)
		and end(index of the last simulation to load into the FompyDataset)
	exclude : array_like
		Index values of simulations to exclude.

	Returns
	-------
	fds : FoMpy Dataset
		Structure of data containing the sanity array.

	"""
	try:
		line_length = []
		index = 0
		for i in fds.sanity_array:
			try:
				if(i!=-1):
					line_length.append(len(fds.dataset[index][:,1][:]))
					index = index + 1
				else:
					line_length.append(0)
			except IndexError:
				line_length.append(0)
		index_crashed = line_length.index(0)
		fds.dataset.insert(index_crashed, np.nan)

		mode_length = mode(line_length)[0]
		for i in range(len(line_length)):
			if(line_length[i] != mode_length):
				fds.sanity_array.pop(i)
				fds.sanity_array.insert(i, -1)
	except ValueError:
		pass

	if (interval is not None):
		temp_dataset= fds.dataset[interval[0]:interval[1]]
		fds.dataset=[]
		fds.dataset=temp_dataset

	if (exclude is not None):
		if(len(exclude)!=1):
			index = 0
			for i in exclude:
				fds.dataset.pop(i-index)
				index = index+1
		else:
			print(exclude)
			fds.dataset.pop(exclude)

	fds.n_sims =len(fds.dataset)

	return fds

