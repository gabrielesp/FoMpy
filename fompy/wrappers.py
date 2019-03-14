
# -*- coding: utf-8 -*-
"""
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
wrappers.py 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This module includes all the wrapper functions used in the FoMpy library.
Code examples on how to use them can be seen either in the complete guide, the repository quickstart
or at the beginning of each file of the source code.

"""

from fompy.fds import daoFile
from fompy.fom import vth_ext, ioff_ext, ion_ext, ss_ext, dibl_ext
from fompy.conditioning import normalizer, filter_tool
from fompy.plots import plotter
from fompy.aux import checkPath
from glob import glob
import os

def dataset(path,filename_user = None, parser = None, save_to_file = None, interval = None, exclude = None, skiprows = None,  comments = None):
	"""
	Wrapper function that creates a FoMpy dataset.

	Parameters
	----------
	globstr : str
		Path to the file containing the IV curves
	parser : void
		Format to read the file
	save_to_file : str
		Path of the file to store all the FoMpy dataset
	interval : array_like
		List of two int values: start(index of the first simulation to load into the FompyDataset)
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
	FompyDataset
		Class containing the most important parameters of a semiconductor IV curve

	"""
	dao_dataset = daoFile()
	fds = dao_dataset.load(path, parser, interval, exclude, skiprows,  comments)
	if(save_to_file != None):
		dao_dataset.save(fds, save_to_file)
	return fds
#----------------------------------------------------------------------------------------------------------------
def extract(fds1, fds2 = None, fom = None, method = None, cc_criteria = None, vg_ext = None, print_fom = None,
 vg_start = None, vg_end = None):
	"""
	Wrapper function that extracts the most common figures of merit from a FoMpy dataset.

	Parameters
	----------
	fds1 : FoMpy Dataset
		Structure of data containing the most important parameters of a semiconductor's IV curve.
		Needed for the extraction of any FoM.
	fds2 : FoMpy Dataset
		additional structure of data containing the most important parameters of a semiconductor's IV curve.
		Needed for the extraction of the DIBL.
	fom : str
		Keyword indicating the desired FoMs to extract. The list of available FoMs includes:
		'vth', 'ioff', 'ion', 'ss', 'ratio', 'power' and 'dibl'.
	method : str
		Keyword indicating the desired method of extraction of the FoMs.The list of available methods includes:
		'SD', 'CC', 'TD' and 'LE'. If method is not defined the 'SD' is selected by default.
	cc_criteria : float, optional
		Float value used for the extraction of several FoMs using the constant current method.
	vg_ext : float
		Gate voltage value used to calculate a FoM like IOFF or ION.
	print_fom : bool
		If True all the FoMs are extracted from a FoMpy Dataset (except for the DIBL)
	vg_start : float
		Gate voltage defining the start of the interval in which the Subthreshold Swing is extracted.
	vg_end : float
		Gate voltage defining the end of the interval in which the Subthreshold Swing is extracted.

	Returns
	-------
	parameter : array_like
		1-d array containing the extracted FoMs.

	""" 
	if(str(fom) == 'vth'):
		temp_vth = vth_ext()
		if (method is not 'LE'):
			parameter_vth, curves= temp_vth.extraction(fds1, method, cc_criteria)
		else:
			parameter_vth, curves, A, B= temp_vth.extraction(fds1, method, cc_criteria)			
		return parameter_vth
	elif(str(fom) == 'ioff'):
		temp_ioff = ioff_ext()
		parameter_ioff, curves= temp_ioff.extraction(fds1,vg_ext = vg_ext)
		return parameter_ioff
	elif(str(fom) == 'ion'):
		temp_ion = ion_ext()
		if(vg_ext is not None) and (type(vg_ext) is float):
			parameter_ion, curves= temp_ion.extraction(fds1,vg_ext = vg_ext)
		else:
			temp_vth = vth_ext()
			if (method is not 'LE'):
				parameter_vth, curves= temp_vth.extraction(fds1, method, cc_criteria)
			else:
				parameter_vth, curves, A, B= temp_vth.extraction(fds1, method, cc_criteria)	
			parameter_ion, curves= temp_ion.extraction(fds1,vth = parameter_vth)
		return parameter_ion
	elif(str(fom) == 'ss'):
		temp_ss = ss_ext()
		parameter_ss, curves, _ = temp_ss.extraction(fds1, vg_start, vg_end)
		return parameter_ss					
	elif(fom == 'dibl'):
		temp = dibl_ext()
		parameter_dibl,curve_high, curve_low, vth_high, vth_low, corriente_low = temp.extraction(fds1,fds2, method)
		return parameter_dibl

	elif(print_fom == True):

		temp_vth = vth_ext()
		if (method is not 'LE'):
			parameter_vth_1, curves= temp_vth.extraction(fds1, method, cc_criteria)
			if(fds2 is not None):
				parameter_vth_2, curves= temp_vth.extraction(fds2, method, cc_criteria)

		else:
			parameter_vth_1, curves, A, B= temp_vth.extraction(fds1, method, cc_criteria)
			if(fds2 is not None):
				parameter_vth_2, curves, A, B= temp_vth.extraction(fds2, method, cc_criteria)
		
		temp_ioff = ioff_ext()
		parameter_ioff_1, curves= temp_ioff.extraction(fds1,vg_ext = vg_ext)
		if(fds2 is not None):
			parameter_ioff_2, curves= temp_ioff.extraction(fds2,vg_ext = vg_ext)

		temp_ion = ion_ext()
		if(vg_ext is not None) and (type(vg_ext) is float):
			parameter_ion_1, curves= temp_ion.extraction(fds1,vg_ext = vg_ext)
			if(fds2 is not None):
				parameter_ion_2, curves= temp_ion.extraction(fds2,vg_ext = vg_ext)
		else:
			temp_vth = vth_ext()
			if (method is not 'LE'):
				parameter_vth_1, curves= temp_vth.extraction(fds1, method, cc_criteria)
				if(fds2 is not None):
					parameter_vth_2, curves= temp_vth.extraction(fds2, method, cc_criteria)

			else:
				parameter_vth_1, curves, A, B= temp_vth.extraction(fds1, method, cc_criteria)
				if(fds2 is not None):
					parameter_vth_2, curves, A, B= temp_vth.extraction(fds2, method, cc_criteria)

			parameter_ion_1, curves= temp_ion.extraction(fds1,vth = parameter_vth_1)
			if(fds2 is not None):
				parameter_ion_2, curves= temp_ion.extraction(fds2,vth = parameter_vth_2)

		temp_ss = ss_ext()
		parameter_ss_1, curves, _ = temp_ss.extraction(fds1, vg_start, vg_end)
		if(fds2 is not None):
			parameter_ss_2, curves, _ = temp_ss.extraction(fds2, vg_start, vg_end)

		if((fds1 is not None) and (fds2 is not None)):
			temp = dibl_ext()
			parameter_dibl,curve_high, curve_low, vth_high, vth_low, corriente_low = temp.extraction(fds1,fds2, method)
			parameter_fds1 = (parameter_vth_1, parameter_ioff_1, parameter_ion_1, parameter_ss_1)
			parameter_fds2 = (parameter_vth_2, parameter_ioff_2, parameter_ion_2, parameter_ss_2)
			print('\nThe extracted values for fds1 are\nvth: {0}\nioff: {1}\nion: {2} \nss: {3}\n'.format(parameter_vth_1, parameter_ioff_1, parameter_ion_1, parameter_ss_1) )
			print('The extracted values for fds2 are\nvth: {0}\nioff: {1}\nion: {2} \nss: {3}\n'.format(parameter_vth_2, parameter_ioff_2, parameter_ion_2, parameter_ss_2) )
			print('The extracted value for the DIBL are\n{0}\n'.format(parameter_dibl) )
			return (parameter_fds1,parameter_fds2)
		else:
			print('\nThe extracted values for fds1 are\nvth: {0}\nioff: {1}\nion: {2} \nss: {3}\n'.format(parameter_vth_1, parameter_ioff_1, parameter_ion_1, parameter_ss_1) )
			return (parameter_vth_1, parameter_ioff_1, parameter_ion_1, parameter_ss_1)
	else:
		raise Exception('The figure of merit hasn\'t been defined!')
#----------------------------------------------------------------------------------------------------------------
def plot(fds1, fds2 = None, plot_type = None, fom = None, parameter=None,  method = None,bins = None, cc_criteria = None, vg_ext = None, vg_start = None, vg_end = None, save_plot = None):
	"""
	Wrapper function that plots the most common figures in semiconductor simulations.

	Parameters
	----------
	fds1 : FoMpy Dataset
		Structure of data containing the most important parameters of a semiconductor's IV curve.
		Needed for generating the plot of any FoM.
	fds2 : FoMpy Dataset
		additional structure of data containing the most important parameters of a semiconductor's IV curve.
		Needed for generating the plot of the calibration and the DIBL.
	plot_type : str
		Keyword indicating the type of plot to generate. The list of available plots includes:
		'iv', 'hist', 'qq', 'varplot', 'calib' and 'fomplot'.
	fom : str
		Keyword indicating the FoM to be plotted. The list of available methods includes:
		'vth', 'ioff', 'ion', 'ss', 'ratio', 'power' and 'dibl'.
	parameter : array_like
		Array of extracted FoM values to be plotted.
	method : str
		Keyword indicating the desired method of extraction of the FoMs. The list of available methods includes:
		'SD', 'CC', 'TD' and 'LE'. If method is not defined the 'SD' is selected by default.
	bins : int
		It defines the number of equal-width bins in the given range.
	cc_criteria : float
		Current criteria used to extract vth with the CC criteria for the fomplot.
	vg_ext : float
		Gate voltage value used to calculate a FoM like IOFF or ION.
	vg_start : float
		Gate voltage defining the start of the interval in which the Subthreshold Swing is extracted.
	vg_end : float
		Gate voltage defining the end of the interval in which the Subthreshold Swing is extracted.
	save_plot : bool
		If True the generated plot is save to the defined path.

	""" 
	plot = plotter()
	if(str(plot_type) is 'iv'):
		plot.iv(fds1, save_plot=save_plot)
	elif(str(plot_type) is 'hist'):
		plot.hist(parameter=parameter,bins=bins, save_plot=save_plot)			
	elif(str(plot_type) is 'qq'):
		plot.qq(parameter=parameter, save_plot=save_plot)	
	elif(str(plot_type) is 'varplot'):
		plot.varplot(fds1, save_plot=save_plot)	
	elif(str(plot_type) == 'calib'):
		plot = plotter()
		plot.calib(fds1,fds2, save_plot=save_plot)	
	elif(plot_type is None):
		if(str(fom) == 'vth'):
			temp_vth = vth_ext()
			if (method is not 'LE'):
				parameter_vth, curves= temp_vth.extraction(fds1, method, cc_criteria)
				temp_vth.plot(fds1 = fds1, parameter = parameter_vth,method=method,cc_criteria=cc_criteria, curves = curves, save_plot = save_plot)
			else:
				parameter_vth, curves, A, B= temp_vth.extraction(fds1, method, cc_criteria)
				temp_vth.plot(fds1 = fds1, parameter = parameter_vth,method=method,cc_criteria=cc_criteria, curves = curves, save_plot = save_plot, A =A, B=B)
		elif(str(fom) == 'ioff'):
			temp_ioff = ioff_ext()
			parameter_ioff, curves= temp_ioff.extraction(fds1,vg_ext)
			temp_ioff.plot(fds1 = fds1, parameter = parameter_ioff, curves = curves,vg_ext = vg_ext, save_plot = save_plot)
		elif(str(fom) == 'ion'):
			temp_ion = ion_ext()
			if(vg_ext is not None) and (type(vg_ext) is float):
				parameter_ion, curves= temp_ion.extraction(fds1,vg_ext = vg_ext)
				print(parameter_ion)
				temp_ion.plot(fds1 = fds1, parameter = parameter_ion, curves = curves,vg_ext = vg_ext, save_plot = save_plot)
			else:
				temp_vth = vth_ext()
				if (method is not 'LE'):
					parameter_vth, curves= temp_vth.extraction(fds1, method, cc_criteria)
				else:
					parameter_vth, curves, A, B= temp_vth.extraction(fds1, method, cc_criteria)	
				parameter_ion, curves= temp_ion.extraction(fds1,vth = parameter_vth)
				temp_ion.plot(fds1 = fds1, parameter = parameter_ion, curves = curves,parameter_vth=parameter_vth, save_plot = save_plot)
		elif(str(fom) == 'ss'):
			temp_ss = ss_ext()
			parameter_ss, curves, vt_sd_medio= temp_ss.extraction(fds1, vg_start, vg_end)
			temp_ss.plot(fds1 = fds1, parameter = parameter_ss, curves = curves, vg_start = vg_start, vg_end = vg_end, vt_sd_medio = vt_sd_medio, save_plot = save_plot)
		elif(str(fom) == 'dibl'):
			temp = dibl_ext()
			parameter_dibl,curve_high, curve_low, vth_high, vth_low, corriente_low  = temp.extraction(fds1,fds2, method)
			temp.plot(fds1, curve_high, curve_low, vth_high, vth_low, corriente_low, save_plot)
#----------------------------------------------------------------------------------------------------------------
def savetotxt(path,fom,parameter):
	"""
	Wrapper function that saves to a text file the extracted FoMs.

	Parameters
	----------
	path : str
		Defines the path where the extracted results are saved to a file.
	fom : str
		Keyword indicating the FoM to be plotted. The list of available methods includes:
		'vth', 'ioff', 'ion', 'ss', 'ratio', 'power' and 'dibl'.
	parameter : array_like
		Array of extracted FoM values to be plotted.

	""" 	
	if(str(fom) == 'vth'):
		temp_vth = vth_ext()
		temp_vth.save_results_to_file(path, parameter)
	elif(str(fom) == 'ioff'):
		temp_ioff = ioff_ext()
		temp_ioff.save_results_to_file(path, parameter)
	elif(str(fom) == 'ion'):
		temp_ion = ion_ext()
		temp_ion.save_results_to_file(path, parameter)
	elif(str(fom) == 'ss'):
		temp_ss = ss_ext()
		temp_ss.save_results_to_file(path, parameter)
	elif(str(fom) == 'dibl'):
		temp_dibl = dibl_ext()
		temp_dibl.save_results_to_file(path, parameter)
#----------------------------------------------------------------------------------------------------------------
def normalize(fds, norm):
	"""
	Wrapper function that normalizes the currents a semiconductor's IV curve.

	Parameters
	----------
	fds : FoMpy Dataset
		Structure of data containing the most important parameters of a semiconductor's IV curve.
	norm : float
		float value used to normalize the currents contained in the FoMpy dataset IV curves.

	""" 	
	fds.norm = norm
	temp_norm = normalizer()
	temp_norm.normalize(fds, norm)
#----------------------------------------------------------------------------------------------------------------
def filter(fds, theta_crit, show_theta = False):
	"""
	Wrapper class for filtering noisy data from a semiconductor's IV curve.

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

	temp_filter = filter_tool()
	fds.dataset = temp_filter.polar_filter(fds, theta_crit=theta_crit, show_theta=show_theta)
#----------------------------------------------------------------------------------------------------------------		
def version():
	"""
	(TODO)Function that prints the current installed version of FoMpy

	"""
	import pkg_resources  # part of setuptools
	version = pkg_resources.require("fompy")[0].version
