# -*- coding: utf-8 -*-
"""
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
fom.py 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This module includes the routines used to extract several commonly used FoMs
in semiconductor simulations.

Example
-------

After the data is imported into a FoMpy Dataset, the user may use several routines implemented in FoMpy to extract the 'vth', 'ioff', 'ion', 'ss', 'ratio', 'power' or 'dibl'. Code examples explaining how to use them can be seen below.

In order to extrac the vth of a given IV curve, the following commands have to be used either in a script or 
in a python3 command line::

	import fompy
	path_file_high = './data/sim_FinFET_vd_high/'
	fds = fompy.dataset(path_file_JCJB, parser=fompy.JCJB)
	fds.drain_bias_label = 'High'
	vth_array = fompy.extract(fds, fom = 'vth')

If you are not using a JCJB file, the attribute of drain_bias_label has to be defined when using the data of a FoMpy dataset at high drain bias, otherwise low drain bias formulas will be used.
Also if a different FoM is needed, the user has to change the keyword fom, from 'vth' to another one contained in the list shown before.

Additionally, if the user wants to obtain the DIBL, two different FoMpy Datasets have to be imported::

	import fompy
	path_file_JCJB = './data/sim_FinFET_vd_high/'
	path_file_low = './data/sim_FinFET_vd_low/'

	fds_hdb = fompy.dataset(path_file_JCJB, parser='JCJB')
	fds_ldb = fompy.dataset(path_file_low, parser='JCJB')
	fds_hdb.drain_bias_value = 0.7
	fds_ldb.drain_bias_value = 0.05
	dibl_array = fompy.extract(fds_hdb, fds_ldb, fom = 'dibl')

"""

from abc import ABCMeta, abstractmethod
from enum import Enum
import numpy as np
from fompy.aux import find_closest, get_diff
from fompy.conditioning import interpolator
from fompy.plots import plotter
import os
from scipy.optimize import curve_fit

class _extractor(metaclass=ABCMeta):
	"""
	Abstract class containing the most common methods applied to an IV curve to extrat the figure of merit.

	"""
	@abstractmethod
	def extraction(self):
		pass
	@abstractmethod
	def save_results_to_file(self):
		pass
	@abstractmethod
	def plot(self):
		pass

#----------------------------------------------------------------------------------------------------------------	

class vth_ext(_extractor):
	"""
	Child class of _extractor that obtains the :math:`V_{TH}` figure of merit from a semiconductor's IV curve.

	"""
#----------------------------------------------------------------------------------------------------------------	
	def extraction(self,fds1,method=None, cc_criteria = None):
		"""
		Methods
		-------
		extraction(fds1,method=None, cc_criteria = None)
			Class method that extracts :math:`V_{TH}` of a semiconductor's IV curve.

		Parameters
		----------
		fds1 : FoMpy Dataset
			Structure of data containing the most important parameters of a semiconductor's IV curve.
			Needed for the extraction of any FoM.
		method : str
			Keyword indicating the desired method of extraction of the FoMs.The list of available methods includes:
			'SD', 'CC', 'TD' and 'LE'. If method is not defined the 'SD' is selected by default.
		cc_criteria : float, optional
			Float value used for the extraction of several FoMs using the constant current method.
		"""
		if(method == 'SD') or(method == None):
			length = len(fds1.dataset)
			parameter = []
			curves = []
			for i in range(length):
				try:
					curve = np.column_stack((fds1.dataset[i][:,0][:], fds1.dataset[i][:,1][:]))
					if(fds1.drain_bias_label=='High'):
						curve[:,1] = np.power( curve[:,1],0.5)
					d2_curve = get_diff(curve, order = 2)
					d2_interp_x, d2_interp_y = interpol(d2_curve[:,0], d2_curve[:,1],strategy = fds1.interpolation, n = 1000,s = 0)
					d2_interp = np.column_stack((d2_interp_x, d2_interp_y))

					lower_bound=find_closest(d2_interp[:,0],d2_interp[50,0])
					higher_bound=find_closest(d2_interp[:,0],d2_interp[-100,0])
					vt_temp = np.argmax(d2_interp[lower_bound:higher_bound,1])
					if (d2_interp[vt_temp+lower_bound, 0] <= d2_interp[lower_bound, 0]):
						print('Vt value outside of the confidence interval for simulation', i)
					if (d2_interp[vt_temp+lower_bound, 0] >= d2_interp[higher_bound, 0]):
						print('Vt value outside of the confidence interval for simulation', i)
					vt_SD = d2_interp[vt_temp+lower_bound,0]

					try:
						vt_temp = np.argmax(d2_interp[lower_bound:higher_bound+lower_bound,1])
						vt_SD = d2_interp[vt_temp+lower_bound,0]
					except ValueError:
						print('Multiple indexes in curve {}'.format(i))

					if(fds1.drain_bias_label=='High'):
						curve[:,1] = np.power( curve[:,1],2)
					parameter.append(float(("%0.3f"%vt_SD)))
					curve = np.column_stack((fds1.dataset[i][:,0][:], fds1.dataset[i][:,1][:]))		
					curves.append(curve)

				except (TypeError, ValueError):
					parameter.append(np.nan)		
					curves.append(np.nan)

			return parameter, curves

		elif(method == 'CC'):
			length = len(fds1.dataset)
			parameter = []
			curves = []
			for i in range(length):
				try:		
					x_interp,y_interp = interpol(fds1.dataset[i][:,0][:], fds1.dataset[i][:,1][:],strategy = fds1.interpolation, n = 1000,s = 0)
					curve = np.column_stack((x_interp, y_interp))
					try:
						vt_index_cc=find_closest(curve[:,1],float(cc_criteria))
						vt_CC_float =curve[vt_index_cc,0]
						parameter.append(float(("%0.3f"%vt_CC_float)))
						curve = np.column_stack((fds1.dataset[i][:,0][:], fds1.dataset[i][:,1][:]))		
						curves.append(curve)
					except (TypeError):
						print('No CC criteria defined')
						break
				except (ValueError, TypeError):
					parameter.append(np.nan)		
					curves.append(np.nan)

			return parameter, curves

		if(method == 'TD') or(method == None):
			length = len(fds1.dataset)
			parameter = []
			curves = []
			for i in range(length):
				try:
					curve = np.column_stack((fds1.dataset[i][:,0][:], fds1.dataset[i][:,1][:]))
					if(fds1.drain_bias_label=='High'):
						curve[:,1] = np.power( curve[:,1],0.5)
					d3_curve = get_diff(curve, order = 3)
					d3_interp_x, d3_interp_y = interpol(d3_curve[:,0], d3_curve[:,1],strategy = fds1.interpolation, n = 1000,s = 0)
					d3_interp = np.column_stack((d3_interp_x, d3_interp_y))

					lower_bound=find_closest(d3_interp[:,0],d3_interp[50,0])
					higher_bound=find_closest(d3_interp[:,0],d3_interp[-100,0])
					vt_temp = np.argmax(d3_interp[lower_bound:higher_bound,1])

					if (d3_interp[vt_temp+lower_bound, 0] <= d3_interp[lower_bound, 0]):
						print('Vt value outside of the confidence interval for simulation', i)
					if (d3_interp[vt_temp+lower_bound, 0] >= d3_interp[higher_bound, 0]):
						print('Vt value outside of the confidence interval for simulation', i)
					vt_TD = d3_interp[vt_temp+lower_bound,0]

					if(fds1.drain_bias_label=='High'):
						curve[:,1] = np.power( curve[:,1],2)						
					parameter.append(float(("%0.3f"%vt_TD)))
					curve = np.column_stack((fds1.dataset[i][:,0][:], fds1.dataset[i][:,1][:]))		
					curves.append(curve)

				except (TypeError, ValueError):
					parameter.append(np.nan)		
					curves.append(np.nan)

			return parameter, curves

		if(method == 'LE') or(method == None):

			def line(x, A, B):
				return A*x + B

			length = len(fds1.dataset)
			parameter = []
			curves = []
			A = []
			B = []
			for i in range(length):
				try:
					x_interp,y_interp = interpol(fds1.dataset[i][:,0][:], fds1.dataset[i][:,1][:],strategy = fds1.interpolation, n = 1000,s = 0)
					curve = np.column_stack((fds1.dataset[i][:,0][:], fds1.dataset[i][:,1][:]))
					
					if(fds1.drain_bias_label=='High'):
						curve[:,1] = np.power( curve[:,1],0.5)

					d1_curve = get_diff(curve, order = 1)
					d1_interp_x, d1_interp_y = interpol(d1_curve[:,0], d1_curve[:,1],strategy = fds1.interpolation, n = 1000,s = 0)
					d1 = np.column_stack((d1_interp_x, d1_interp_y))

					index_lower_bound=find_closest(d1[:,0],d1[50,0])
					index_higher_bound=find_closest(d1[:,0],d1[-100,0])

					x_interp, y_interp = interpol(fds1.dataset[i][:,0][:], fds1.dataset[i][:,1][:],strategy = fds1.interpolation, n = 1000,s = 0)
					curve = np.column_stack((x_interp, y_interp))

					try:
						d1_init = find_closest(x_interp,d1[0,0])
						d1_max_index = np.argmax(d1[index_lower_bound:index_higher_bound,1])+index_lower_bound+d1_init
						vt_temp_d1 = find_closest(x_interp,d1[d1_max_index,0])
					except ValueError:
						print('Multiple indexes')

					corriente_vt_temp =curve[vt_temp_d1,1]
					lim1 = int(vt_temp_d1*0.9)
					lim2 = int(vt_temp_d1*1.1)

					A_i,B_i = curve_fit(line, curve[lim1:lim2,0],curve[lim1:lim2,1])[0]
					fit = A_i*x_interp+B_i
					vt_index_le_fit = find_closest(fit,0)

					if(fds1.drain_bias_label=='High'):
						try:
							vt_LE = x_interp[vt_index_le_fit]
							curve[:,1] = np.power( curve[:,1],2)
						except NameError:
							print('The drain bias value has not been defined')
					else:
						vd_low = fds1.drain_bias_value
						vt_LE = curve[vt_index_le_fit,0]+vd_low/2
				
					parameter.append(float(("%0.3f"%vt_LE)))
					
					curve = np.column_stack((fds1.dataset[i][:,0][:], fds1.dataset[i][:,1][:]))
					curves.append(curve)
					A.append(A_i)
					B.append(B_i)
				except (TypeError, ValueError, IndexError):
					parameter.append(np.nan)		
					curves.append(np.nan)
					A.append(np.nan)
					B.append(np.nan)					
			return parameter, curves, A, B

#----------------------------------------------------------------------------------------------------------------	
	def save_results_to_file(self,path, parameter):
		"""
		save_results_to_file(path, parameter)
			Class method that saves the extracted :math:`V_{TH}` values.

		Parameters
		----------
		path : str
			Defines the path where the extracted results are saved to a file.
		parameter : array_like
			Array of extracted FoM values to be saved into the file.

		"""

		length = len(parameter)
		f = open(str(path), "w+")
		f.write("#Simulation_id\t#Parameter\n")
		for i in range(length):
			if (str(type(parameter)) is not '<class \'numpy.ndarray\'>'):
				parameter = np.array(parameter)
			arr = np.array2string(parameter[i], separator="\t")
			f.write("{0}\t{1}\n".format(i, arr))			
		f.close()
#----------------------------------------------------------------------------------------------------------------		
	def plot(self, fds1, parameter = None, method = None,cc_criteria = None, curves = None,backend= None, save_plot = None, A=None, B=None):
		"""
		plot(fds1, parameter = None, method = None,cc_crit = None, curves = None, save = None, A=None, B=None)
			Class method that plots the extracted :math:`V_{TH}` values.

		Parameters
		----------
		fds1 : FoMpy Dataset
			Structure of data containing the most important parameters of a semiconductor's IV curve.
			Needed for generating the plot of any FoM.
		parameter : array_like
			Array of extracted FoM values to be plotted.
		method : str
			Keyword indicating the desired method of extraction of the FoMs. The list of available methods includes: 'SD', 'CC', 'TD' and 'LE'. If method is not defined the 'SD' is selected by default.
		cc_criteria : float
				Current criteria used to extract vth with the CC criteria for the fomplot.
		curves : array_like
			Array of data containing the IV curves.
		backend : str
			String containing the name of the backend chosen to either plot or save the plots. The backends available are:
			'Agg', which only works whenever saving plots to files (non-GUI) and 'TkAgg' a GUI tools for visualizing the plots.
			'TkAgg' requires the package python3-tk installed in order to run.			
		save_plot : bool
			If True the generated plot is save to the defined path.
		A, B: float
			Parameters obtained during the vth LE extraction method used for the plots.

		"""

		temp_plot = plotter()
		for i in range(len(fds1.dataset)):
			try:
				if (method is not 'LE'):
					temp_plot.fomplot(i,fds1, fom = 'vth', voltages = curves[i][:,0], currents = curves[i][:,1],
				 	parameter = parameter[i], method = method,cc_criteria = cc_criteria,backend= backend, save_plot = save_plot)
				else:
					temp_plot.fomplot(i,fds1, fom = 'vth', voltages = curves[i][:,0], currents = curves[i][:,1],
				 	parameter = parameter[i], method = method,cc_criteria = cc_criteria,backend= backend, save_plot = save_plot, A = A[i], B=B[i])
			except (TypeError, ValueError):
				i = i+i


#----------------------------------------------------------------------------------------------------------------	

class ioff_ext(_extractor):
	"""
	Child class of _extractor that obtains the :math:`I_{OFF}` figure of merit from a semiconductor's IV curve.

	"""
#----------------------------------------------------------------------------------------------------------------	
	def extraction(self,fds1 = None,vg_ext = None):
		"""
		Methods
		-------
		extraction(fds1,method=None, cc_criteria = None)
			Class method that extracts :math:`I_{OFF}` of a semiconductor's IV curve.

		Parameters
		----------
		fds1 : FoMpy Dataset
			Structure of data containing the most important parameters of a semiconductor's IV curve.
			Needed for the extraction of any FoM.
		vg_ext : float
			Gate voltage value used to calculate IOFF at.
		"""
		parameter = []
		curves = []
		for i in range(len(fds1.dataset)):
			try:
				x_interp,y_interp = interpol(fds1.dataset[i][:,0][:], fds1.dataset[i][:,1][:],strategy=fds1.interpolation, n = 1000,s = 0)
				single_curve = np.column_stack((x_interp, y_interp))
				if(vg_ext is not 0) and (vg_ext is not None) :
					voltage_index = find_closest(x_interp, vg_ext)
					parameter.append(y_interp[voltage_index])
					curves.append(single_curve)

				else:
					parameter.append(y_interp[0])
					curves.append(single_curve)
			except (TypeError, ValueError):
				parameter.append(np.nan)		
				curves.append(np.nan)
		return parameter, curves
#----------------------------------------------------------------------------------------------------------------	
	def save_results_to_file(self,path, parameter):
		"""
		save_results_to_file(path, parameter)
			Class method that saves the extracted :math:`I_{OFF}` values.

		Parameters
		----------
		path : str
			Defines the path where the extracted results are saved to a file.
		parameter : array_like
			Array of extracted FoM values to be saved into the file.

		"""
		length = len(parameter)
		f = open(str(path), "w+")
		f.write("#Simulation id\t#Parameter\n")
		f.write("##############################\n###")
		for i in range(length):
			if (str(type(parameter)) is not '<class \'numpy.ndarray\'>'):
				parameter = np.array(parameter)
			arr = np.array2string(parameter[i], separator="\t")
			f.write("{0}\t{1}\n".format(i, arr))			
		f.close()
#----------------------------------------------------------------------------------------------------------------	
	def plot(self, fds1, parameter = None,vg_ext = None, curves = None,backend= None, save_plot = None):
		"""
		plot(fds1, parameter = None, method = None,cc_crit = None, curves = None, save = None, A=None, B=None)
			Class method that plots the extracted :math:`I_{OFF}` values.

		Parameters
		----------
		fds1 : FoMpy Dataset
			Structure of data containing the most important parameters of a semiconductor's IV curve.
			Needed for generating the plot of any FoM.
		parameter : array_like
			Array of extracted FoM values to be plotted.
		vg_ext : float
			Gate voltage value used to calculate IOFF at.
		curves : array_like
			Array of data containing the IV curves.
		backend : str
			String containing the name of the backend chosen to either plot or save the plots. The backends available are:
			'Agg', which only works whenever saving plots to files (non-GUI) and 'TkAgg' a GUI tools for visualizing the plots.
			'TkAgg' requires the package python3-tk installed in order to run.
		save_plot : bool
			If True the generated plot is save to the defined path.

		"""
		temp_plot = plotter()
		for i in range(len(fds1.dataset)):
			try:
				temp_plot.fomplot(i,fds1, fom = 'ioff', voltages = curves[i][:,0], currents = curves[i][:,1], parameter = parameter[i],vg_ext = vg_ext,backend= backend, save_plot = save_plot)
			except (TypeError, ValueError):
				pass

class ion_ext(_extractor):
	"""
	Child class of _extractor that obtains the :math:`I_{ON}` figure of merit from a semiconductor's IV curve.

	"""

	def extraction(self,fds1,vg_ext = None, vth = None):
		"""
		Methods
		-------
		extraction(fds1,vg_ext = None, vth = None)
			Class method that extracts :math:`I_{ON}` of a semiconductor's IV curve.

		Parameters
		----------
		fds1 : FoMpy Dataset
			Structure of data containing the most important parameters of a semiconductor's IV curve.
			Needed for the extraction of any FoM.
		vg_ext : float
			Gate voltage value used to calculate IOFF at.
		vth : array_like
			Array of vth extracted values used for obtaining ION (using the vth-dependant formula).

		"""

		if(vth is not None):
			length = len(fds1.dataset)
			parameter = []
			curves = []
			for i in range(length):
				try:
					x_interp,y_interp = interpol(fds1.dataset[i][:,0][:], fds1.dataset[i][:,1][:],strategy=fds1.interpolation, n = 1000,s = 0)
					single_curve = np.column_stack((x_interp, y_interp))				
					ion_index=find_closest(x_interp, vth[i]+fds1.drain_bias_value)	
					parameter.append(y_interp[ion_index][0])
					curves.append(single_curve)
				except (TypeError, ValueError):
					parameter.append(np.nan)		
					curves.append(np.nan)
		elif(vg_ext is not None) and (vth is None):
			length = len(fds1.dataset)
			parameter = []
			curves = []
			for i in range(length):
				try:
					x_interp,y_interp = interpol(fds1.dataset[i][:,0][:], fds1.dataset[i][:,1][:],strategy=fds1.interpolation, n = 1000,s = 0)
					single_curve = np.column_stack((x_interp, y_interp))				
					ion_index=find_closest(x_interp, vg_ext)	
					parameter.append(y_interp[ion_index])
					curves.append(single_curve)
				except (TypeError, ValueError):
					parameter.append(np.nan)		
					curves.append(np.nan)
		else:
			print('No criteria for the extracion of ION has been defined')			
				
		return parameter, curves
#----------------------------------------------------------------------------------------------------------------	
	def save_results_to_file(self,path, parameter):
		"""
		save_results_to_file(path, parameter)
			Class method that saves the extracted :math:`I_{ON}` values.

		Parameters
		----------
		path : str
			Defines the path where the extracted results are saved to a file.
		parameter : array_like
			Array of extracted FoM values to be saved into the file.

		"""
		length = len(parameter)
		f = open(str(path), "w+")
		f.write("#Simulation id\t#Parameter\n")
		f.write("##############################\n###")
		for i in range(length):
			if (str(type(parameter)) is not '<class \'numpy.ndarray\'>'):
				parameter = np.array(parameter)
			arr = np.array2string(parameter[i], separator="\t")
			f.write("{0}\t{1}\n".format(i, arr))			
		f.close()
#----------------------------------------------------------------------------------------------------------------	
	def plot(self, fds1 , parameter = None, curves = None, parameter_vth = None, vg_ext = None,backend= None, save_plot = None):
		"""
		plot(fds1, parameter = None, method = None,cc_crit = None, curves = None, save = None, A=None, B=None)
			Class method that plots the extracted :math:`I_{ON}` values.

		Parameters
		----------
		fds1 : FoMpy Dataset
			Structure of data containing the most important parameters of a semiconductor's IV curve.
			Needed for generating the plot of any FoM.
		parameter : array_like
			Array of extracted FoM values to be plotted.
		curves : array_like
			Array of data containing the IV curves.
		parameter_vth : array_like
			Array of extracted vth values, as a method to obtain ION depends on them.		
		vg_ext : float
			Gate voltage value used to calculate IOFF at.
		backend : str
			String containing the name of the backend chosen to either plot or save the plots. The backends available are:
			'Agg', which only works whenever saving plots to files (non-GUI) and 'TkAgg' a GUI tools for visualizing the plots.
			'TkAgg' requires the package python3-tk installed in order to run.			
		save_plot : bool
			If True the generated plot is save to the defined path.

		"""		
		temp_plot = plotter()
		for i in range(len(fds1.dataset)):
			try:
				if(vg_ext is not None) and (type(vg_ext) is float):
					temp_plot.fomplot(i,fds1, fom = 'ion', voltages = curves[i][:,0], currents = curves[i][:,1], parameter = parameter[i], vg_ext =vg_ext,backend= backend, save_plot = save_plot)
				else:
					temp_plot.fomplot(i,fds1, fom = 'ion', voltages = curves[i][:,0], currents = curves[i][:,1], parameter = parameter[i], parameter_vth = parameter_vth[i],backend= backend, save_plot = save_plot)
			except (TypeError, ValueError):
				pass

class ss_ext(_extractor):
#----------------------------------------------------------------------------------------------------------------	
	def extraction(self, fds1, vth = None, vg_start = None, vg_end = None):
		"""
		Methods
		-------
		extraction( fds1, vg_start = None, vg_end = None)
			Class method that extracts :math:`SS` of a semiconductor's IV curve.

		Parameters
		----------
		fds1 : FoMpy Dataset
			Structure of data containing the most important parameters of a semiconductor's IV curve.
			Needed for the extraction of any FoM.
		vth : array_like
			Array of vth extracted values used for obtaining SS (using the vth-dependant formula).			
		vg_start : float
			Gate voltage defining the start of the interval in which the Subthreshold Swing is extracted.
		vg_end : float
			Gate voltage defining the end of the interval in which the Subthreshold Swing is extracted.

		"""
		parameter = []
		curves = []
		vt_sd_medio = []
		for i in range(len(fds1.dataset)):
			try:
				curve = np.column_stack((fds1.dataset[i][:,0][:], fds1.dataset[i][:,1][:]))

				x_interp, y_interp = interpol(curve[:,0], curve[:,1],strategy = fds1.interpolation, n = 1000,s = 0)
				if(vg_start is None) and (type(vg_end) is float):
					index = find_closest(x_interp,vg_end)
					parameter.append((x_interp[0]-x_interp[index])*1000/(np.log10(y_interp[0])-np.log10(y_interp[index])))
				elif(vg_start is not None) and (vg_end is not None) :
					index_start = find_closest(x_interp,vg_start)
					index_end = find_closest(x_interp,vg_end)
					parameter.append((x_interp[index_start]-x_interp[index_end])*1000/(np.log10(y_interp[index_start])-np.log10(y_interp[index_end])))
				elif(vg_start is not None) and (type(vg_start) is float):
					index_start = find_closest(x_interp,vg_start)
					index = find_closest(x_interp,vth[i]/2)

					parameter.append((x_interp[index_start]-x_interp[index])*1000/(np.log10(y_interp[index_start])-np.log10(y_interp[index])))
				else:
					index = find_closest(x_interp,vth[i]/2)
					parameter.append((x_interp[0]-x_interp[index])*1000/(np.log10(y_interp[0])-np.log10(y_interp[index])))
				# print(parameter)
				curves.append(curve)
				vt_sd_medio.append(vth[i]/2)
			except (TypeError, ValueError):
				parameter.append(np.nan)
				curves.append(np.nan)
				vt_sd_medio.append(np.nan)

		return parameter, curves, vt_sd_medio
#----------------------------------------------------------------------------------------------------------------	
	def save_results_to_file(self,path, parameter):
		"""
		save_results_to_file(path, parameter)
			Class method that saves the extracted :math:`SS` values.

		Parameters
		----------
		path : str
			Defines the path where the extracted results are saved to a file.
		parameter : array_like
			Array of extracted FoM values to be saved into the file.

		"""		
		length = len(parameter)
		f = open(str(path), "w+")
		f.write("#Simulation id\t#Parameter\n")
		f.write("##############################\n###")
		for i in range(length):
			if (str(type(parameter)) is not '<class \'numpy.ndarray\'>'):
				parameter = np.array(parameter)
			arr = np.array2string(parameter[i], separator="\t")
			f.write("{0}\t{1}\n".format(i, arr))			
		f.close()
# ----------------------------------------------------------------------------------------------------------------	
	def plot(self, fds1, parameter = None, curves = None, parameter_ss = None, vg_start = None, vg_end = None, vt_sd_medio = None,backend= None, save_plot = None):
		"""
		plot(fds1, parameter = None, method = None,cc_crit = None, curves = None, save = None, A=None, B=None)
			Class method that plots the extracted :math:`SS` values.

		Parameters
		----------
		fds1 : FoMpy Dataset
			Structure of data containing the most important parameters of a semiconductor's IV curve.
			Needed for generating the plot of any FoM.
		parameter : array_like
			Array of extracted FoM values to be plotted.
		curves : array_like
			Array of data containing the IV curves.
		parameter_ss : array_like
			Array of extracted ss values.		
		vg_start : float
			Gate voltage defining the start of the interval in which the Subthreshold Swing is extracted.
		vg_end : float
			Gate voltage defining the end of the interval in which the Subthreshold Swing is extracted.
		vg_sd_medio : float
			Value in the middle of the interval between zero and vth extracted with the SD method.
			It is used only for defining a limit in the plot.
		backend : str
			String containing the name of the backend chosen to either plot or save the plots. The backends available are:
			'Agg', which only works whenever saving plots to files (non-GUI) and 'TkAgg' a GUI tools for visualizing the plots.
			'TkAgg' requires the package python3-tk installed in order to run.			
		save_plot : bool
			If True the generated plot is save to the defined path.

		"""
		temp_plot = plotter()
		for i in range(len(fds1.dataset)):
			try:
				temp_plot.fomplot(i,fds1, fom = 'ss', voltages = curves[i][:,0], currents = curves[i][:,1], parameter = parameter[i],vg_start = vg_start, vg_end = vg_end, vt_sd_medio = vt_sd_medio[i] ,backend= backend,save_plot = save_plot)
			except (TypeError, ValueError):
					pass

class dibl_ext(_extractor):

	def extraction(self,fds1,fds2,method = None, cc_criteria = None):
		"""
		Methods
		-------
		extraction(fds1,fds2,method = None)
			Class method that extracts the :math:`DIBL` of two given semiconductor's IV curves.

		Parameters
		----------
		fds1 : FoMpy Dataset
			Structure of data containing the most important parameters of a semiconductor's IV curve.
			Needed for the extraction of any FoM.
		fds2 : FoMpy Dataset
			additional structure of data containing the most important parameters of a semiconductor's IV curve.
			Needed for generating the plot of the calibration and the DIBL.
		method : str
			Keyword indicating the desired method of extraction of the FoMs. The list of available methods includes:
			'SD', 'CC', 'TD' and 'LE'. If method is not defined the 'SD' is selected by default.

		"""

		temp_vth = vth_ext()
		if(method == None):
			method = 'SD'
		if(method is not 'LE'):
			# fds1.drain_bias_label = 'High'
			vth_high, curve_high= temp_vth.extraction(fds1, method, cc_criteria)
			vth_low, curve_low= temp_vth.extraction(fds2, method, cc_criteria)
		else:
			# fds1.drain_bias_label = 'High'
			vth_high, curve_high, _, _= temp_vth.extraction(fds1, method, cc_criteria)
			vth_low, curve_low, _, _= temp_vth.extraction(fds2, method, cc_criteria)

		dibl = []
		corriente_low_arr = []
		for i in range(len(fds1.dataset)):
			try:
				vth_index_low = find_closest(curve_low[i][:,0],vth_low[i])
				corriente_low = curve_low[i][vth_index_low,1]

				# print(vth_high)
				# print(vth_low)
				# print(fds1.drain_bias_value)
				# print(fds2.drain_bias_value)
				# print(vth_high)
				# print(vth_low)
				dibl.append((-(vth_low[i] - vth_high[i])/(fds1.drain_bias_value-fds2.drain_bias_value)*1000))
				corriente_low_arr.append(corriente_low)
			except (TypeError, ValueError):
				dibl.append(np.nan)
				corriente_low_arr.append(np.nan)

		return dibl,curve_high, curve_low, vth_high, vth_low, corriente_low_arr

#----------------------------------------------------------------------------------------------------------------	
	def save_results_to_file(self,path, parameter):
		"""
		save_results_to_file(path, parameter)
			Class method that saves the extracted :math:`DIBL` values.

		Parameters
		----------
		path : str
			Defines the path where the extracted results are saved to a file.
		parameter : array_like
			Array of extracted FoM values to be saved into the file.

		"""		
		length = len(parameter)
		f = open(str(path), "w+")
		f.write("#Simulation id\t#Parameter\n")
		f.write("##############################\n###")
		for i in range(length):
			if (str(type(parameter)) is not '<class \'numpy.ndarray\'>'):
				parameter = np.array(parameter)
			arr = np.array2string(parameter[i], separator="\t")
			f.write("{0}\t{1}\n".format(i, arr))			
		f.close()
#----------------------------------------------------------------------------------------------------------------	
	def plot(self,fds1 , curve_high = None,curve_low = None, parameter_vt_high = None, parameter_vt_low = None,corriente_low = None,backend= None, save_plot = None):
		"""
		plot(fds1, parameter = None, method = None,cc_crit = None, curves = None, save = None, A=None, B=None)
			Class method that plots the extracted :math:`DIBL` values.

		Parameters
		----------
		fds1 : FoMpy Dataset
			Structure of data containing the most important parameters of a semiconductor's IV curve.
			Needed for generating the plot of any FoM.
		curve_high : array_like
			Array of data containing the IV curves at high drain bias.
		curve_low : array_like
			Array of data containing the IV curves at low drain bias.			
		parameter_vt_high : float
			Voltage value of high drain bias.
		parameter_vt_low : float
			Voltage value of low drain bias.
		corriente_low : float
			Current at the vth value extracted for the curve at low drain bias.
		backend : str
			String containing the name of the backend chosen to either plot or save the plots. The backends available are:
			'Agg', which only works whenever saving plots to files (non-GUI) and 'TkAgg' a GUI tools for visualizing the plots.
			'TkAgg' requires the package python3-tk installed in order to run.
		save_plot : bool
			If True the generated plot is save to the defined path.

		"""
		temp_plot = plotter()
		for i in range(len(fds1.dataset)):
			try:	
			# print(i)
				temp_plot.fomplot(i,fds1, fom = 'dibl', curve_high = curve_high[i], curve_low = curve_low[i], vth_high = parameter_vt_high[i], vth_low = parameter_vt_low[i], corriente_low = corriente_low[i],backend= backend, save_plot = save_plot)
			except (TypeError, ValueError):
				pass

#----------------------------------------------------------------------------------------------------------------	

def interpol(x = None ,y = None, n = None, strategy = None, d = None, s = None):
	"""
	Wrapper function for interpolating imported data from a semiconductor's IV curve.

	Parameters
	----------
	x : array_like, shape (n,)
		1-d array containing values of the independent variable.
	y : array_like
		Array containing values of the dependent variable.
		It can have arbitrary number of dimensions, but the length along axis
		must match the length of x. Values must be finite.
	strategy : str
		Keyword for defining the selected interpolation method: The list of available methods includes:
		'akima', 'pchip' and 'linear'.
	d : int
		Degree of the smoothing spline. Must be <= 5.
		Default is k=3, a cubic spline.
	s : float
		Positive smoothing factor used to choose the number of knots. 
	
	"""

	if(strategy is None) or (strategy is 'cubic_spline'):
		temp_interp = interpolator()
		x_interp, y_interp = temp_interp.spline_interpol(x, y, n, d, s)

	elif(strategy is 'akima'):
		temp_interp = interpolator()
		x_interp, y_interp = temp_interp.akima_interpol(x, y, n)

	elif(strategy is 'pchip'):
		temp_interp = interpolator()
		x_interp, y_interp = temp_interp.pchip_interpol(x, y, n)

	elif(strategy is 'linear'):
		temp_interp = interpolator()
		x_interp, y_interp = temp_interp.lin_interpol(x, y, n)

	return x_interp, y_interp
