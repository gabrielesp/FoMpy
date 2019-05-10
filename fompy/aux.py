# -*- coding: utf-8 -*-
"""
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
aux.py 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This module includes several auxiliar routines used during the usage of FoMpy.

Example
-------
In order find the closest index between a given value and a array of values, the function 
find_closest is used. For example during the extraction of IOFF::

	voltage_index = find_closest(x_interp, vg_ext)
	parameter.append(y_interp[voltage_index])

In order calculate the discrete derivate of a given curve the function get_diff is used as follows (second derivative method
extraction)::

	d2 = get_diff(curve, order = 2)
	lower_bound=find_closest(d2[:,0],10*d2[1,0])
	higher_bound=find_closest(d2[2:,0],0.8*d2[-1,0])						
	warning_interval_limit_low = int(np.round(lower_bound*1.5))
	warning_interval_limit_high = int(np.round(higher_bound*0.8))
	vt_temp = np.argmax(d2[lower_bound:higher_bound,1])

where these the maximum value of the second derivative is found leaving out values outside a defined interval where the derivatives
can be too noisy.

Finally, the function checkPath is used in FoMpy to check if a directory exists and if it does not it creates it when saving a plot or results.

"""

import numpy as np
from abc import ABCMeta, abstractmethod
import os

def find_closest(A, target):
	"""
	Function that finds the closest item inside an array A to a target value.
	It returns the index of A with the closest value to the target (:math:`A[idx]~~{\\approx}~~target`).

	Parameters
	-------
	A : array_like
		Array containing all the values to be compared. It usually is
		an array of currents of a FoMpy dataset.
	target : float
		Target value to be found inside the array A.
	Returns
	-------
	idx : int
		Index of the item contained in the array A with a value closest to the target.
	"""
	idx = A.searchsorted(target) #A must be sorted
	idx = np.clip(idx, 1, len(A)-1)
	left = A[idx-1]
	right = A[idx]
	idx -= target - left < right - target
	return idx
#----------------------------------------------------------------------------------------------------------------
def get_diff(arr, order, type = 'central'):
	"""
	Function for calculating the n-order derivative of a discrete array

	Parameters
	-------
	order : int
		Order indicating the type of derivative to compute
	type: str
		The type of derivate can either be calculated in a forward way, 
		central or backward, depending on the evaluation coordinates
		of the derivative function (see https://en.wikipedia.org/wiki/Finite_difference)
	"""

	x = arr[:,0]
	y = arr[:,1]

	dx = np.gradient(x)
	dy = np.gradient(y)
	d1 = dy/dx
	d2y = np.gradient(d1)
	d2 = d2y/dx
	d3y = np.gradient(d2y)
	d3 = d3y/dx

	if order == 1:
		gm_c = np.column_stack((x, d1))
		return(gm_c)

	if order == 2:
		gm2_c = np.column_stack((x, d2))
		return(gm2_c)

	if order == 3:
		gm3_c = np.column_stack((x, d3))
		return(gm3_c)
#----------------------------------------------------------------------------------------------------------------
def checkPath(cwd):
	"""
	Function that checks if a given path exists, and if it doesn't, it creates it.

	Parameters
	-------
	cwd : str
		Path that the user wishes to check if it exists.
	"""
	try:
		os.makedirs(cwd)
	except OSError:
		if not os.path.isdir(cwd):
			raise
