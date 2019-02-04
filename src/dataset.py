
from abc import ABCMeta, abstractmethod
from glob import glob
from src.fom_methods import *
import numpy as np

class _FompyDataset(metaclass=ABCMeta):
	'Factory interface containing the methods to import data as a Fompy Dataset'
	def __init__(self, **kwargs):
		self.dataset = []
		self.n_sims = None
		self.norm = None
		self.ext_method = None
		self.drain_bias = None
		self.drain_bias_value = None

	def print_parameters(self):
		"""
		Prints the details of a FompyDataset.
		"""		
		print("\nNumber of simulations: ", self.n_sims)
		print("Normalization value: ", self.norm)
		print("Extraction method: ", self.ext_method)
		print("Drain Bias: ", self.drain_bias)
		print("Drain Bias value: ", self.drain_bias_value, "\n")

	def extraction(self,fom = None):

		if(str(fom.name) == 'ioff') or (fom == None):
			temp = ioff()
			ioff_fds3 = temp.extract(self,method=ioffMethod.default)
			print('The {0} array has been extracted'.format(fom.name))
		return ioff_fds3


class FompyDataset(_FompyDataset):
	'Class containing the methods and attributes to instantiate a Fompy Dataset '

	def __init__(self,**kwargs):
		_FompyDataset.__init__(self, **kwargs)
		for key in kwargs:
			if(kwargs[key] == drain_bias.Low):
				raise Exception('FompyDataset_HDB cannot have attribute drain bias set to Low!')
			else:
				setattr(self, key, kwargs[key])


class dataDAO(metaclass=ABCMeta):

	@abstractmethod
	def load(self):
		pass
	@abstractmethod
	def save(self):
		pass

class daoJCJB(dataDAO, FompyDataset):
	'Create or save a dataset from a set of files'

	def load(self, globstr):
		fds = FompyDataset()
		filenames = sorted(glob(globstr))
		fds.n_sims = len(filenames)
		for path in filenames:
			data_temp = np.loadtxt(path)
			data = extractJCJB(data_temp[:,0], data_temp[:,1], data_temp[:,2])
			arr = np.column_stack((data[:,1], data[:,2]))
			fds.dataset.append(arr)
		return fds

	def save(self, fds, path):
		length = len(fds.dataset)
		f = open(str(path), "w+")
		for i in range(length):
			f.write("Simulation id {0}\n".format(i))
			f.write("Vd[V]\tVg[V]\tIs[V]\t?[V]\tId[V]\n")
			arr = np.array2string(fds.dataset[i][:], separator="\t")
			f.write(arr+"\n")
		f.close()

class daoArray(dataDAO, FompyDataset):

	def load(self, arraydata):
		pass
	def save(self):
		pass

class daoFile(dataDAO, FompyDataset):

	def load(self, arraydata):
		pass
	def save(self):
		pass


#===========================================================================#
 #    ///\    |||   |||   \\\ ///       ||||||  |||   |||  |||\\  |||  |||||||
 #   /// \\   |||   |||    \\|//        |||     |||   |||  ||| \\ |||  |||
 #  ///---\\  |||   |||    //|\\        |||||   |||   |||  |||  \\|||  |||
 # ///     \\ |||||||||   /// \\\       |||     |||||||||  |||   \|||  |||||||
#===========================================================================#

# This is a wrapper that prevents from constructing the DAO explicitly
# when using the library

def JCJBtoDataset(globstr):

	dao_dataset = daoJCJB()
	fds = dao_dataset.load(globstr)
	print('The data has been loaded from file "%s"' % (globstr))
	return fds

def datasettoFile(fds,globstr):

	dao_dataset = daoJCJB()
	dao_dataset.save(fds, globstr)
	print('The data has been save to file in "%s"' % (globstr))
	return fds

def extractJCJB(a, b, c):

	condition = np.mod(a, np.amax(a))==0
	vd_temp = np.copy(np.extract(condition, a))
	vg_temp = np.copy(np.extract(condition, b))
	id_temp = np.copy(np.extract(condition, c))
	Vd = np.copy(vd_temp[np.argmin(vg_temp):len(vg_temp)])
	Vg = np.copy(vg_temp[np.argmin(vg_temp):len(vg_temp)])
	Id = np.copy(id_temp[np.argmin(vg_temp):len(vg_temp)])
	arr = np.column_stack((Vd, Vg, Id))
	return arr		

class drain_bias(Enum):
	'This class implements the drain bias used in the simulation'
	High, Low = range(2)

