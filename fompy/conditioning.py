# -*- coding: utf-8 -*-
"""
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
conditioning.py 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This module includes the routines used to precondition the input data
before any extraction is done like interpolating the
data points, normalizing the data and filtering points that are too noisy.

Example
-------
In order to normalize data the following steps have to be executed either in a script or
python3 terminal::

	import fompy
	path_file_JCJB = './data/sim_FinFET_vd_high/'
	fds = fompy.dataset(path_file_JCJB, parser=fompy.JCJB)
	norm_value = 35.8/10**9
	fompy.normalize(fds, norm_value)

In order to filter data the following steps have to be executed either in a script or
python3 terminal::

	import fompy
	path_file_JCJB = './data/sim_FinFET_vd_high/'
	fds = fompy.dataset(path_file_JCJB, parser=fompy.JCJB)
	fompy.filter(fds, theta_crit = 1.52)

Additionally, by default FoMpy uses the cubic spline interpolation in order to 
extract accurately the FoMs. If the user wishes to change the interpolation methods, 
there are several implemented and tested methods that can be selected. 
The list of available methods includes: 'akima','pchip' or 'linear'::

	import fompy
	path_file_JCJB = './data/sim_FinFET_vd_high/'
	fds = fompy.dataset(path_file_JCJB, parser=fompy.JCJB)
	fds.interpolation = 'linear'

"""


from abc import ABCMeta, abstractmethod
import numpy as np
#----------------------------------------------------------------------------------------------------------------
class _normalizerStrategy(metaclass=ABCMeta):
	"""
	Abstract class for normalizing a semiconductor's IV curve.

	"""

	@abstractmethod
	def normalize(self):
		pass

#----------------------------------------------------------------------------------------------------------------
class normalizer(_normalizerStrategy):
	"""
	Normalizer class.

	...

	Methods
	-------
	normalize(fds, norm)
		Class method used for normalizing the currents a semiconductor's IV curve.

	Parameters
	----------
	fds : FoMpy Dataset
		Structure of data containing the most important parameters of a semiconductor's IV curve.
	norm : float
		float value used to normalize the currents contained in the FoMpy dataset IV curves.
	"""

	def normalize(self, fds, norm):
		try:
			for i in range(len(fds.dataset)):
				fds.dataset[i][:,1][:] = np.divide(fds.dataset[i][:,1][:],norm)
		except TypeError:
			pass
#----------------------------------------------------------------------------------------------------------------
class _interpolationStrategy(metaclass=ABCMeta):
	"""
	Abstract class containing several tested interpolation strategies that can be used in a semiconductor's IV curve.

	"""

	@abstractmethod
	def spline_interpol(self, x,y,s):
		pass

	@abstractmethod
	def akima_interpol(self,  x,y,s):
		pass        

	@abstractmethod
	def pchip_interpol(self,  x,y,s):
		pass  

	@abstractmethod
	def lin_interpol(self,x,y):
		pass  

#----------------------------------------------------------------------------------------------------------------
class interpolator(_interpolationStrategy):
	"""
	Class containing several tested interpolation strategies that can be used in a semiconductor's IV curve.
	For further documentation go to  https://docs.scipy.org/doc/scipy/reference/interpolate.html

	...

	Methods
	-------
	spline_interpol(x ,y , n , d , s)
		Interpolate data with a piecewise cubic polynomial which is twice
		continuously differentiable [1]. The result is represented as a PPoly
		instance with breakpoints matching the given data. The natural boundary
		condition is selected by default.

	Parameters
	----------
	x : array_like, shape (n,)
		1-d array containing values of the independent variable.
	y : array_like
		Array containing values of the dependent variable.
		It can have arbitrary number of dimensions, but the length along axis
		must match the length of x. Values must be finite.
	n : int
		Length of the output interpolated array.
	d : int
		Degree of the smoothing spline. Must be <= 5.
		Default is k=3, a cubic spline.
	s : float
		Positive smoothing factor used to choose the number of knots. 

	Methods
	-------
	akima_interpol(x,y,n)
		Fit piecewise cubic polynomials, given vectors x and y.
		The interpolation method by Akima uses a continuously differentiable
		sub-spline built from piecewise cubic polynomials. The resultant curve passes
		through the given data points and will appear smooth and natural.

	Parameters
	----------
	x : array_like, shape (n,)
		1-D array of monotonically increasing real values.
	y : array_like
		N-D array of real values. The length of y
		along the first axis must be equal to the length of x.
	n : int
		Length of the output interpolated array.

	Methods
	-------
	pchip_interpol(x, y, n)
		Convenience function for pchip interpolation. xi and yi are arrays of
		values used to approximate some function f, with yi = f(xi).
		The interpolant uses monotonic cubic splines to find the value of
		new points x and the derivatives there.

	Parameters
	----------
	x : array_like, shape (n,)
		A 1-D array of monotonically increasing real values.
		x cannot include duplicate values (otherwise f is overspecified)
	y : array_like
		A 1-D array of real values. y’s length along
		the interpolation axis must be equal to the length
		of x. If N-D array, use axis parameter to select correct axis.
	n : int
		Length of the output interpolated array.

	Methods
	-------
	lin_interpol(x,y, n)
		Interpolate a 1-D function.x and y are arrays of values used to approximate
		some function f: y = f(x). This class returns a function
		whose call method uses interpolation to find the value of new points.

	Parameters
	----------
	x : array_like, shape (n,)
		A 1-D array of real values.
	y : array_like
		A N-D array of real values. The length of y along
		the interpolation axis must be equal to the length of x.
	n : int
		Length of the output interpolated array.

	"""
	def spline_interpol(self, x = None ,y = None, n = None, d = None, s = 0):

		if(d is 3) or (d is None):
			from scipy.interpolate import CubicSpline
			cs_natural = CubicSpline(x, y, bc_type='natural')
			xvals = np.linspace(x[0], x[-1], n)
			return (xvals, cs_natural(xvals))
		elif(d is 5) :
			from scipy import interpolate
			spline = interpolate.splrep(x, y, k=5)
			xvals = np.linspace(x[0], x[-1], n)
			yinterp = interpolate.splev(xvals, spline, der=0)
			return (xvals, yinterp)
		else:
			pass

	def akima_interpol(self, x,y,n):
		from scipy.interpolate import Akima1DInterpolator
		cs_natural = Akima1DInterpolator(x, y, axis=0)
		xvals = np.linspace(x[0], x[-1], n)
		return (xvals, cs_natural(xvals))

	def pchip_interpol(self, x, y, n):
		from scipy.interpolate import PchipInterpolator
		cs_natural = PchipInterpolator(x, y)
		xvals = np.linspace(x[0], x[-1], n)		
		return (xvals, cs_natural(xvals)) 

	def lin_interpol(self,x,y, n):
		from scipy.interpolate import interp1d
		cs_natural = interp1d(x, y)
		xvals = np.linspace(x[0], x[-1], n)		
		return (xvals, cs_natural(xvals))

#----------------------------------------------------------------------------------------------------------------
class _filterStrategy(metaclass=ABCMeta):
	"""
	Abstract class for filtering data of a semiconductor's IV curve.

	"""

	@abstractmethod
	def polar_filter(self):
		pass

	# @abstractmethod
	# def gauss_filter(self):
	# 	pass

#--------------------------------------------------------------------------------------------------------------
class filter_tool(_filterStrategy):
	"""
	Class for filtering noisy data from a semiconductor's IV curve.

	...

	Methods
	-------
	polar_filter(fds,theta_crit, show_theta = False)
		Class method that filters data from a semiconductor's IV curve checking the increase in the angle
		of the points with respect to the origin of an IV curve and removes all the increments that go above
		defined a threshold.

	Parameters
	----------
	fds : FoMpy Dataset
		Structure of data containing the most important parameters of a semiconductor's IV curve.
	theta_crit : float
		Threshold value for the biggest allowed increase in the angle between to consecutive points
		of an IV curve.
	show_theta : bool
		If True it prints all the angle increments between two consecutive points so the user can
		choose a suitable theta crit.	
	"""

	def polar_filter(self,fds,theta_crit, show_theta = False):
		for i in range(len(fds.dataset)):
			try:
				theta = angle_wrt_0(fds.dataset[i][:,0],np.log10(fds.dataset[i][:,1])) #The angle wrt 0 is calculated
				step = np.gradient(angle_wrt_0(fds.dataset[i][:,0],np.log10(fds.dataset[i][:,1]))) #The difference between the angles is calculated
				
				j = 0
				step_lst = list(np.abs(step))
				if(show_theta is True):
					print(step_lst)
				indexes = []
				for index, elem in list(enumerate(step_lst)):
					if elem > theta_crit: #If the angle increment is above threshold the index is stored in 'indexes'
						indexes.append(index-j)
						step_lst.pop(index-j)
						j = j+1
			except (TypeError, ValueError):
				pass

			for item in indexes: #And then deleted
				fds.dataset[i] = np.delete(fds.dataset[i],item, axis=0)


#----------------------------------------------------------------------------------------------------------------
def angle_wrt_0(x, y):
	"""Function that calculates the polar angle of a point with respect to the origin of coordinates.
	It returns the angle in degrees and it is assummed that the angle is contained in the 1rst cuadrant.

	Parameters
	----------
	x: float
		The x-axis coodinate of the point.
	y: float
		The y-axis coodinate of the point.

	Returns
	-------
	theta: float
		Angle between the given point and the origin
		of coordinates in degrees. It is assummed that the angle is contained
		in the 1rst cuadrant.

	"""
	if(len(x) is not len(y)):
		raise('The sizes of x and y must be the same')
	if(len(x) is 1):
		if(x == 0):
			theta = 90
		else:
			theta = np.arctan(y/x)*180/np.pi
	else:
		theta = []
		for i in range(len(x)):
			if(x[i] == 0):
				theta.append(90)
			else:
				theta.append(np.arctan(y[i]/x[i])*180/np.pi)
	return np.abs(theta) #Theta is considered to be in the first cuadrant
