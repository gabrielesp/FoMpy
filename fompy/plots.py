# -*- coding: utf-8 -*-
"""
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
plots.py 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This module includes the routines used to generate several common figures
in semiconductor simulations.

Example
-------

After the FoMs have been extracted, FoMpy includes several useful visualizing routines, 
commonly used in semiconductor simulations. Code examples explaining how to use them can be seen below.

In order to generate a plot of the FoMs (fomplot)::

	import fompy
	path_file_high = './data/sim_FinFET_vd_high/'
	fds = fompy.dataset(path_file_JCJB, parser=fompy.JCJB)
	vth_array = fompy.extract(fds, fom = 'vth')
	fompy.plot(fds, fom = 'vth')

and FoMpy will generate a single graph for each simulation showing the IV curve with its correspondent FoM extraction criteria and the extracted FoM.

Most of the FoMs are plotted the same way, except for the DIBL. This is because, if we want to calculate the DIBL
two FoMpy Datasets are needed::

	import fompy
	path_file_JCJB = './data/sim_FinFET_vd_high/'
	path_file_low = './data/sim_FinFET_vd_low/'

	fds_hdb = fompy.dataset(path_file_JCJB, parser='JCJB')
	fds_ldb = fompy.dataset(path_file_low, parser='JCJB')
	fds_hdb.drain_bias_value = 0.7
	fds_ldb.drain_bias_value = 0.05

	fompy.plot(fds_hdb, fds_ldb, fom = 'dibl')

Additionally, other types of plots can be generated. The list of available plots includes: 'iv', 'hist', 'qq', 'varplot', 'calib' and 'fomplot'. These are some examples::

	path_file_var = './data/simulations/'
	fds_var = fompy.dataset(path_file_var, parser=fompy.JCJB)
	vth_array = fompy.extract(fds_var, fom = 'vth')

	fompy.plot(fds_var, type='hist', parameter=vth_array, bins=5)
	fompy.plot(fds_var, type='qq', parameter=vth_array)

	path_file_var = './data/simulations/'
	fds_var = fompy.dataset(path_file_var, parser=fompy.JCJB)
	fompy.plot(fds_var, type='varplot')


	path_file_JCJB = './data/sim_FinFET_vd_high/'
	path_file_low = './data/sim_FinFET_vd_low/'

	fds_hdb = fompy.dataset(path_file_JCJB, parser=fompy.JCJB)
	fds_hdb.print_parameters()
	fds_ldb = fompy.dataset(path_file_low, parser=fompy.JCJB)
	fds_ldb.print_parameters()

	norm_value = 35.8/10**9
	fompy.normalize(fds_hdb, norm_value)
	fompy.normalize(fds_ldb, norm_value)

	fompy.plot(fds_hdb,fds_ldb, type='calib')
"""

from abc import ABCMeta, abstractmethod
import numpy as np
from scipy import stats
from enum import Enum
import os
import probscale
from fompy.aux import find_closest, get_diff, checkPath
from fompy.conditioning import interpolator
import matplotlib.pyplot as plt


class plotStrategy(metaclass=ABCMeta):
	"""
	Abstract class containing several commonly used plots for semiconductor simulations.

	"""

	# def __init__(self, **kwargs):
	# 	self.type = None

	@abstractmethod
	def iv(self):
		pass

	@abstractmethod
	def hist(self):
		pass

	@abstractmethod
	def qq(self):
		pass

	@abstractmethod
	def varplot(self):
		pass

	@abstractmethod
	def calib(self):
		pass

	@abstractmethod
	def fomplot(self):
		pass	
class plotter(plotStrategy):
	"""
	Class used for generating common plots in semiconductor simulations.
	For extensive documentation on how to modify this code go to https://matplotlib.org/tutorials/index.html
	"""

	def iv(self,fds,backend = None, save_plot = None):
		"""
		Methods
		-------
		iv(fds, save_plot = None)
			Class method that filters data from a semiconductor's IV curve using Gaussian filtering.

		Parameters
		----------
		fds : FoMpy Dataset
			Structure of data containing the most important parameters of a semiconductor's IV curve.
		save_plot : path or None, optional
			Path indicating the folder where the user wishes to save the generated plots.

		"""

		import matplotlib
		if(backend is None) or (backend is 'Agg'):
			matplotlib.use('Agg')
		elif(backend is 'TkAgg'):
			matplotlib.use('TkAgg')
		else:
			pass

		for i in range(len(fds.dataset)):
			try:
				plt.close()
				fig, ax1 = plt.subplots()
				ax2 = ax1.twinx()
				ax1.plot(fds.dataset[i][:,0][:],fds.dataset[i][:,1][:],lw =0.4,color='red', label='Log')
				ax2.plot(fds.dataset[i][:,0][:],fds.dataset[i][:,1][:],lw =0.4,color='green', label='Lin')
				ax1.set_yscale('log')
				ax1.set_ylabel(r'$\mathrm{I_{D}} $',labelpad=20)
				ax1.set_xlabel(r'$\mathrm{V_{G}} $',labelpad=20)
				ax1.set_title('IV curve plot -> Curve %s' % str(i))
				handles2, labels2 = ax1.get_legend_handles_labels()
				labels2, ids2 = np.unique(labels2, return_index=True)
				ax1.legend(handles2, labels2, loc='best')
				ax2.set_yscale('linear')
				ax2.yaxis.set_visible(False)
				handles, labels = ax2.get_legend_handles_labels()
				labels, ids = np.unique(labels, return_index=True)
				ax2.legend(handles, labels, loc='lower right')
				plt.tight_layout()
				if(save_plot is not None):
					try:
						checkPath(save_plot)
						plt.savefig(save_plot+'iv_{0}.pdf'.format(str(i)), bbox_inches='tight', format='pdf', dpi=1200 )
					except NameError:
						print('No path to save the plots has been chosen')
				else:
					print('No path has been defined\nTrying to use GUI backend...')
					plt.show()
				plt.close()
			except(TypeError, ValueError):
				pass
		return


	def hist(self, bins = None, parameter = None,cont_parameter = None, backend = None, save_plot = None):
		"""
		Methods
		-------
		hist(bins = None, parameter = None, save_to_file = None):
			Plot a histogram.The return value is a tuple (n, bins, patches) or ([n0, n1, ...],
			bins, [patches0, patches1,...]) if the input contains multiple data.

		Parameters
		----------
		bins : int
			If an integer is given, bins + 1 bin edges are calculated and returned, consistent with 
			numpy.histogram
		parameter : array_like, shape (n,)
			Input values, this takes either a single array or a sequence of arrays
			which are not required to be of the same length.
		save_plot : path or None, optional
			Path indicating the folder where the user wishes to save the generated plots.
		"""

		import matplotlib
		if(backend is None) or (backend is 'Agg'):
			matplotlib.use('Agg')
		elif(backend is 'TkAgg'):
			matplotlib.use('TkAgg')
		else:
			pass

		if(len(parameter)<2):
			raise Exception('Not enough data to plot a histogram')

		parameter = [x for x in parameter if x == x] # Remove all NaN values from the parameter list
		plt.close()
		xmin = np.min(parameter)
		xmax = np.max(parameter)
		x_range = np.linspace(xmin*0.8, xmax*1.2, len(parameter))
		m, s = stats.norm.fit(parameter) # get mean and standard deviation
		pdf_g = stats.norm.pdf(x_range, m, s) # now get theoretical probability density function in our interval
		max_pdf = np.max(pdf_g)
		pdf_g = pdf_g/max_pdf
		weights = np.ones_like(parameter)/float(len(parameter))
		if(s==0):
			print('The normal fit has not been plotted as sigma is zero')
		if(bins==None):
			bins = 5 # default number of bines if parameter not defined
		if(cont_parameter is not None):
			plt.axvline(cont_parameter, color='blue', linestyle='dashed', linewidth=2, label='Ideal case') # plot the extracter parameter value for the ideal case
		w, bins, patches = plt.hist(parameter,bins,weights=weights, facecolor='green')
		plt.plot(x_range, np.max(w)*pdf_g, color = 'r', label='Normal Fit') # plot the normal distribution fit
		plt.axvline(np.mean(parameter), color='r', linestyle='dashed', linewidth=2)
		plt.title('Distribution parameters: '+r'$\mu='+str(("%0.3f"%m))+r', \sigma$='+str(("%0.3f"%s))+'\n'  )
		plt.xlabel('Parameter')
		plt.ylabel('Frecuency')
		plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		plt.tight_layout()
		if(save_plot is not None):
			try:
				checkPath(save_plot)
				plt.savefig(save_plot+'hist.pdf', bbox_inches='tight', format='pdf', dpi=1200 )
			except NameError:
				print('No path to save the plots has been chosen')
		else:
			print('No path has been defined\nTrying to use GUI backend...')
			plt.show()
		plt.close()

		return

	def qq(self, parameter,backend = None, save_plot = None):
		"""
		Methods
		-------
		qq( parameter, save_to_file = None):
			Plot a Quantile plot.Plotting positions are converted into quantiles
			or Z-scores based on a probability distribution

		Parameters
		----------

		parameter : array_like, shape (n,)
			Input values, this takes either a single array or a sequence of arrays
			which are not required to be of the same length.
		save_plot : path or None, optional
			Path indicating the folder where the user wishes to save the generated plots.		

		"""

		import matplotlib
		if(backend is None) or (backend is 'Agg'):
			matplotlib.use('Agg')
		elif(backend is 'TkAgg'):
			matplotlib.use('TkAgg')
		else:
			pass

		if(len(parameter)<2):
			raise Exception('Not enough data for a QQ-plot')
			
		parameter = [x for x in parameter if x == x] # Remove all NaN values from the parameter list
		plt.close()
		fig, ax = plt.subplots()
		slope, r = probscale.probplot(parameter, ax=ax, plottype='qq', bestfit=True, estimate_ci=False, return_best_fit_results=True,problabel='Standard Normal Quantiles',
		datalabel=r'Parameter',
		scatter_kws={'markersize': 10, 'linestyle': 'none', 'label': None})
		slope, intercept, rcoef, prob, sterrest = stats.linregress(r['x'],r['y'])
		# plt.text(1.5,max(r['y']), "$R^2_{"+str(method)+"}=%1.2f$" % rcoef**2)
		# plt.legend(loc='lower right', fontsize='medium')
		plt.title('QQ-plot')
		if(save_plot is not None):
			try:
				checkPath(save_plot)
				plt.savefig(save_plot+'qqplot.pdf', bbox_inches='tight', format='pdf', dpi=1200 )
			except NameError:
				print('No path to save the plots has been chosen')
		else:
			print('No path has been defined\nTrying to use GUI backend...')
			plt.show()
		plt.close()

		return

	def varplot(self,fds,backend = None, save_plot = None):
		"""
		Methods
		-------
		varplot(self,fds, save_plot = None):
			Plot all the IV curves for a common variability source.

		Parameters
		----------
		fds : FoMpy Dataset
			Structure of data containing the most important parameters of a semiconductor's IV curve.
		save_plot : path or None, optional
			Path indicating the folder where the user wishes to save the generated plots.

		"""
		import matplotlib
		if(backend is None) or (backend is 'Agg'):
			matplotlib.use('Agg')
		elif(backend is 'TkAgg'):
			matplotlib.use('TkAgg')
		else:
			pass

		plt.close()
		fig, ax1 = plt.subplots()	
		ax2 = ax1.twinx()		
		for i in range(len(fds.dataset)):
			try:		
				ax1.plot(fds.dataset[i][:,0][:],fds.dataset[i][:,1][:],lw =0.4,color='red', label='Log')
				ax2.plot(fds.dataset[i][:,0][:],fds.dataset[i][:,1][:],lw =0.4,color='green', label='Lin')
				ax1.set_yscale('log')
				ax1.set_ylabel(r'$\mathrm{I_{D}} $',labelpad=20)
				ax1.set_xlabel(r'$\mathrm{V_{G}} $',labelpad=20)
				ax1.set_title('Variability plot')
				ax2.yaxis.set_visible(False)	
				plt.tight_layout() 	
			except(TypeError, ValueError):
				pass					
		if(save_plot is not None):
			try:
				checkPath(save_plot)
				plt.savefig(save_plot+'varplot.pdf', bbox_inches='tight', format='pdf', dpi=1200 )
			except NameError:
				print('No path to save the plots has been chosen')
		else:
			print('No path has been defined\nTrying to use GUI backend...')
			plt.show()
		plt.close()
		return

	def calib(self,fds1,fds2,backend = None, save_plot = None):

		import matplotlib
		if(backend is None) or (backend is 'Agg'):
			matplotlib.use('Agg')
		elif(backend is 'TkAgg'):
			matplotlib.use('TkAgg')
		else:
			pass

		for i in range(len(fds1.dataset)):
			try:
				plt.close()
				fig, ax1 = plt.subplots()
				ax2 = ax1.twinx()
				ax1.plot(fds1.dataset[i][:,0][:],fds1.dataset[i][:,1][:],'o',linestyle= '-.',markerfacecolor="None",color='red', label=r'$\mathrm{V_{D}}=$'+str(fds1.drain_bias_value)+' V')
				ax2.plot(fds1.dataset[i][:,0][:],fds1.dataset[i][:,1][:],'o',linestyle= '-.',markerfacecolor="None",color='red')
				ax1.plot(fds2.dataset[i][:,0][:],fds2.dataset[i][:,1][:],'s',linestyle= '-.',markerfacecolor="None",color='blue', label=r'$\mathrm{V_{D}}=$'+str(fds2.drain_bias_value)+' V')
				ax2.plot(fds2.dataset[i][:,0][:],fds2.dataset[i][:,1][:],'s',linestyle= '-.',markerfacecolor="None",color='blue')
				ax1.set_yscale('log')
				ax1.set_ylabel('$I_{D} [\mu A/\mu m]$')
				ax1.set_xlabel('$V_{G} [V]$')
				ax1.set_title('Calibration plot')
				ax1.legend()
				plt.tight_layout()
			except(TypeError, ValueError):
				plt.close()	
		if(save_plot is not None):
			try:
				checkPath(save_plot)
				plt.savefig(save_plot+'calibration.pdf', bbox_inches='tight', format='pdf', dpi=1200 )
			except NameError:
				print('No path to save the plots has been chosen')
		else:
			print('No path has been defined\nTrying to use GUI backend...')
			plt.show()
		plt.close()
		return

	def fomplot(self,i,fds1, fom = None,  currents = None, voltages = None, parameter = None,method = None, cc_criteria = None,parameter_vth = None, vg_ext = None, curve_high = None, curve_low = None, vth_high = None, vth_low = None, corriente_low = None,backend = None, save_plot = None, A=None, B=None, vg_start = None, vg_end = None, vt_sd_medio = None):
		"""
		Methods
		-------
		fomplot(i, fom = None,  voltages = None, currents = None, parameter = None,method = None, cc_criteria = None,parameter_vth = None,
		vg_ext = None, curve_high = None, curve_low = None, vth_high = None, vth_low = None, corriente_low = None, save_plot = None, A=None, B=None, vg_start = None, vg_end = None, vt_sd_medio = None):

			Plot the most common figures of merit of a semiconductor's IV curve.

		Parameters
		----------
		fom : FoMpy Dataset
			Structure of data containing the most important parameters of a semiconductor's IV curve.
		currents : path or None, optional
			Path indicating the folder where the user wishes to save the generated plots.
		voltages : path or None, optional
			Path indicating the folder where the user wishes to save the generated plots.			
		backend : str
			String containing the name of the backend chosen to either plot or save the plots. The backends available are:
			'Agg', which only works whenever saving plots to files (non-GUI) and 'TkAgg' a GUI tools for visualizing the plots.'TkAgg' requires the package python3-tk installed in order to run.
		"""
		import matplotlib

		if(backend is None) or (backend is 'Agg'):
			matplotlib.use('Agg')
		elif(backend is 'TkAgg'):
			matplotlib.use('TkAgg')
		else:
			pass

		if(currents is not None):
			if(len(currents)<2):
				raise Exception('Not enough data for a variability curve')
		else:
			if(len(curve_high)<2) or (len(curve_low)<2):
				raise Exception('Not enough data for a variability curve')

		if(fom == None):
			raise Exception('The figure of merit to plot wasn\'t defined!')

		#Plot parameters
		curve_color = 'dodgerblue' #Color of the curve
		curve_label = 'Simulated Data' #Label of the curve
		x_label = r'$\mathrm{V_{G}} $' #Name of the x label
		y_label = r'$\mathrm{I_{D}} $' #Name of the y label
		label_pad = 20 #Padding applied to the x and y labels
		yscale = "log" #Type of scale["linear", "log", "symlog", "logit", ...]
		txt_fontsize = 15 #Fontsize of the plot
		txt_family = 'serif' #Font family of the plot
		txt_style = 'normal' #Font style of the plot

		if(fom == 'vth'):

			if(method == 'SD') or(method == None):

				title=r'$\mathrm{V_{TH}} $ SD extraction method ->Curve %s' % str(i) #Title of the plot
				SD_color = 'pink' #Color of the Second Derivative curve
				SD_linestyle = '--' #Line Style of the Second Derivative curve
				SD_label = 'SD ($d^{2}I_{D}/dV_{G}^{2}$)' #Label of the Second Derivative curve
				line_color = 'blue'
				line_style = '-.'
				line_label = r'$\mathrm{V_{TH}}$' #Label of the FoM vertical line	

				plt.close()
				fig, ax1 = plt.subplots()
				ax2 = ax1.twinx()
				curve = np.column_stack((voltages,currents))
				if(fds1.drain_bias_label=='High'):
					curve[:,1] = np.power( curve[:,1],0.5)
				ax1.plot(curve[:,0], curve[:,1],color = curve_color,label = curve_label)
				ax1.set_yscale(yscale)
				ax1.set_xlabel(x_label,labelpad=label_pad)
				ax1.set_ylabel(y_label,labelpad=label_pad)
				textstr = '$V_{T}=%.3f$'%(parameter)
				props = dict(boxstyle="square,pad=0.1", fc='w', alpha=1,ec='k', lw =2)
				ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=txt_fontsize,  va='top', bbox=props,family = txt_family, style = txt_style)
				d2_curve = get_diff(curve, order = 2)
				d2_interp_x, d2_interp_y = interpol(d2_curve[:,0], d2_curve[:,1],strategy = fds1.interpolation, n = 1000,s = 0)
				d2_interp = np.column_stack((d2_interp_x, d2_interp_y))
				ax2.plot(d2_interp[:, 0],d2_interp[:, 1], color=SD_color, linestyle= SD_linestyle ,label = SD_label)
				ax2.yaxis.set_visible(False) #The scale of the second derivative is removed
				ax1.axvline(x = parameter,color=line_color,label = line_label) #The value of the vth is represented with a vertical line
				ax1.legend(loc='best',fontsize = txt_fontsize)
				plt.title(title)
				plt.tight_layout()
				plt.show()
				if(save_plot is not None):
					try:
						checkPath(save_plot)
						plt.savefig(save_plot+'vth_SD_{0}.pdf'.format(str(i)), bbox_inches='tight', format='pdf', dpi=1200 )
					except NameError:
						print('No path to save the plots has been chosen')
				else:
					print('No path has been defined\nTrying to use GUI backend...')
					plt.show()
				plt.close()

			if(method == 'CC'):
				title=r'$\mathrm{V_{TH}} $ CC extraction method ->Curve %s' % str(i) #Title of the plot
				CC_color = 'pink' #Color of the CC curve
				CC_linestyle = '--' #Line Style of the CC curve
				CC_label = 'CC criteria' #Label of the CC curve
				line_color = 'blue'
				line_style = '-.'
				line_label = r'$\mathrm{V_{TH}}$' #Label of the FoM vertical line					

				plt.close()
				fig, ax1 = plt.subplots()			
				ax1.plot(voltages,currents,color = curve_color,label = curve_label)
				ax1.set_yscale(yscale)					
				ax1.set_xlabel(x_label,labelpad=label_pad)
				ax1.set_ylabel(y_label,labelpad=label_pad)
				textstr = '$V_{T}=%.3f$'%(parameter)
				props = dict(boxstyle="square,pad=0.1", fc='w', alpha=1,ec='k', lw =2)
				ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=txt_fontsize,  va='top', bbox=props,family = txt_family, style = txt_style)
				ax1.axhline(y=cc_criteria, color=CC_color, linestyle=CC_linestyle, label = CC_label)
				ax1.legend(loc='best',fontsize =txt_fontsize)
				plt.title(title)
				plt.tight_layout()				
				if(save_plot is not None):
					try:
						checkPath(save_plot)
						plt.savefig(save_plot+'vth_CC_{0}.pdf'.format(str(i)), bbox_inches='tight', format='pdf', dpi=1200 )
					except NameError:
						print('No path to save the plots has been chosen')
				else:
					print('No path has been defined\nTrying to use GUI backend...')
					plt.show()
				plt.close()

			if(method == 'TD'):
				title=r'$\mathrm{V_{TH}} $ TD extraction method ->Curve %s' % str(i) #Title of the plot
				TD_color = 'pink' #Color of the TD curve
				TD_linestyle = '--' #Line Style of the TD curve
				TD_label = 'SD ($d^{3}I_{D}/dV_{G}^{3}$)' #Label of the TD curve
				line_color = 'blue'
				line_style = '-.'
				line_label = r'$\mathrm{V_{TH}}$' #Label of the FoM vertical line	
				path_to_save = str(save_plot) + 'vth_TD'+'%03d'%i+'.pdf' #Path to save to plot to

				plt.close()
				fig, ax1 = plt.subplots()
				ax2 = ax1.twinx()
				curve = np.column_stack((voltages,currents))
				if(fds1.drain_bias_label=='High'):
					curve[:,1] = np.power( curve[:,1],0.5)
				ax1.plot(curve[:,0], curve[:,1],color = curve_color,label = curve_label)
				ax1.set_yscale(yscale)					
				ax1.set_xlabel(x_label,labelpad=label_pad)
				ax1.set_ylabel(y_label,labelpad=label_pad)
				textstr = '$V_{T}=%.3f$'%(parameter)
				props = dict(boxstyle="square,pad=0.1", fc='w', alpha=1,ec='k', lw =2)
				ax2.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=14,  va='top', bbox=props,family = txt_family, style = txt_style)

				d3_curve = get_diff(curve, order = 3)
				d3_interp_x, d3_interp_y = interpol(d3_curve[:,0], d3_curve[:,1],strategy = fds1.interpolation, n = 1000,s = 0)
				ax2.plot(d3_interp_x, d3_interp_y, color=TD_color, linestyle= TD_linestyle ,label = TD_label)
				ax2.yaxis.set_visible(False) #The scale of the second derivative is removed
				ax1.axvline(x = parameter,color=line_color, linestyle = line_style, label = line_label) #The value of the vth is represented with a vertical line
				ax1.legend(loc='best',fontsize = txt_fontsize)
				plt.title(title)
				plt.tight_layout()			
				if(save_plot is not None):
					try:
						checkPath(save_plot)
						plt.savefig(save_plot+'vth_TD_{0}.pdf'.format(str(i)), bbox_inches='tight', format='pdf', dpi=1200 )
					except NameError:
						print('No path to save the plots has been chosen')
				else:
					print('No path has been defined\nTrying to use GUI backend...')
					plt.show()
				plt.close()

			if(method == 'LE'):
				title=r'$\mathrm{V_{TH}} $ LE extraction method ->Curve %s' % str(i) #Title of the plot
				LE_color = 'pink' #Color of the LE curve
				LE_deriv_color = 'red' #Color of the LE curve
				LE_linestyle = '--' #Line Style of the LE curve
				LE_label = 'LE criteria' #Label of the LE curve
				LE_deriv_label ='1rst deriv' #Label of the first derivative
				line_color = 'blue'
				line_style = '-.'
				line_label = r'$\mathrm{V_{TH}}$' #Label of the FoM vertical line				

				plt.close()
				fig, ax1 = plt.subplots()
				ax2 = ax1.twinx()
				curve = np.column_stack((voltages, currents))
				if(fds1.drain_bias_label=='High'):
					curve[:,1] = np.power( curve[:,1],0.5)
				d1_curve = get_diff(curve, order = 1)
				d1_interp_x, d1_interp_y = interpol(d1_curve[:,0], d1_curve[:,1],strategy = fds1.interpolation, n = 1000,s = 0)
				x_interp, y_interp = interpol(voltages, currents,strategy = fds1.interpolation, n = 1000,s = 0)				
				ax1.plot(voltages,currents,color = curve_color,label = curve_label)
				fit = A*x_interp+B
				vt_index_le_fit = find_closest(fit,0)
				ax1.plot(x_interp, A*x_interp+B,color=LE_color,linestyle =LE_linestyle, label=LE_label)
				ax2.plot(d1_interp_x, d1_interp_y,color=LE_deriv_color,linestyle =LE_linestyle, label=LE_deriv_label)
				ax1.set_xlabel(x_label,labelpad=label_pad)
				ax1.set_ylabel(y_label,labelpad=label_pad)
				textstr = '$V_{T}=%.3f$'%(parameter)
				props = dict(boxstyle="square,pad=0.1", fc='w', alpha=1,ec='k', lw =2)
				ax2.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=14,  va='top', bbox=props,family = txt_family, style = txt_style)
				ax1.legend(loc='best')
				ax2.yaxis.set_visible(False)
				ax1.axvline(x = parameter,color=line_color, linestyle = line_style, label = line_label) #The value of the vth is represented with a vertical line
				ax1.axhline(y = 0,color=line_color, linestyle = line_style, label = line_label) #The value of the vth is represented with a vertical line				
				plt.title(title)
				plt.tight_layout()
				if(save_plot is not None):
					try:
						checkPath(save_plot)
						plt.savefig(save_plot+'vth_LE_{0}.pdf'.format(str(i)), bbox_inches='tight', format='pdf', dpi=1200 )
					except NameError:
						print('No path to save the plots has been chosen')
				else:
					print('No path has been defined\nTrying to use GUI backend...')
					plt.show()
				plt.close()

		if(fom == 'ioff'):
			title=r'$\mathrm{I_{OFF}[A]}$ extraction ->Curve %s' % str(i) #Title of the plot
			IOFF_color = 'pink' #Color of the IOFF curve
			IOFF_linestyle = '--' #Line Style of the IOFF curve
			IOFF_label = r'$\mathrm{I_{OFF}}[A]$' #Label of the IOFF curve
			line_color = 'blue'
			line_style = '-.'
			line_label = r'$\mathrm{I_{OFF}}$' #Label of the FoM vertical line
				
			plt.close()
			fig, ax1 = plt.subplots()
			ax1.set_yscale(yscale)	
			ax1.plot(voltages,currents,color = curve_color,label = curve_label)
			ax1.set_xlabel(x_label,labelpad=label_pad)
			ax1.set_ylabel(y_label,labelpad=label_pad)			
			ax1.axhline(y=float(parameter), color=IOFF_color, linestyle=IOFF_linestyle, label = IOFF_label)
			index = find_closest(currents, parameter)
			ax1.axvline(x = float(voltages[index]), color=line_color, linestyle=line_style)
			if(vg_ext is not None):
				ax1.axvline(x = vg_ext, color=line_color, linestyle=line_style, label = line_label)
			ax1.legend(loc='best',fontsize =25)
			plt.title(title)
			plt.tight_layout()

			if(save_plot is not None):
				try:
					checkPath(save_plot)
					plt.savefig(save_plot+'IOFF_{0}.pdf'.format(str(i)), bbox_inches='tight', format='pdf', dpi=1200 )
				except NameError:
					print('No path to save the plots has been chosen')
			else:
				print('No path has been defined\nTrying to use GUI backend...')
				plt.show()
			plt.close()

		if(fom == 'ion'):
			plt.close()
			fig, ax1 = plt.subplots()
			ax1.set_yscale("log")
			ax1.plot(voltages,currents,color='orangered',label ='Simulated Data')
			ax1.axhline(y=float(parameter), color='blue', linestyle='-.', label =r'$\mathrm{I_{ON}}[A]$')
			index = find_closest(currents, parameter)
			ax1.axvline(x = float(voltages[index]), color='blue', linestyle='-.')
			if(vg_ext is not None):
				ax1.axvline(x = vg_ext, color='blue', linestyle='-.')
			ax1.set_ylabel(r'$\mathrm{I_{D}} $',labelpad=20)
			ax1.set_xlabel(r'$\mathrm{V_{G}} $',labelpad=20)
			ax1.legend(loc='best',fontsize =25)
			plt.title(r'$\mathrm{I_{ON}[A]}$ extraction ->Curve %s' % str(i))
			plt.tight_layout()

			if(save_plot is not None):
				try:
					checkPath(save_plot)
					plt.savefig(save_plot+'ION_{0}.pdf'.format(str(i)), bbox_inches='tight', format='pdf', dpi=1200 )
				except NameError:
					print('No path to save the plots has been chosen')
			else:
				print('No path has been defined\nTrying to use GUI backend...')
				plt.show()
			plt.close()

		if(fom == 'ss'):
			plt.close()
			fig, ax1 = plt.subplots()
			#ax1.set_yscale("log")
			plt.plot(voltages, np.log10(currents), 'o')
			center = find_closest(voltages,(voltages[0]+voltages[-1])/2)
			if(vg_start is None) and (type(vg_end) is float):
				index = find_closest(voltages,vg_end)
				slope, intercept, _, _, _ = stats.linregress(voltages[0:index],np.log10(currents[0:index]))
			elif(vg_start is not None) and (vg_end is not None) and (type(vg_end) is float) and (type(vg_start) is float):
				index_start = find_closest(voltages,vg_start)
				index_end = find_closest(voltages,vg_end)
				slope, intercept, _, _, _ = stats.linregress(voltages[index_start:index_end],np.log10(currents[index_start:index_end]))
			elif(vg_start is not None) and (type(vg_start) is float):
				index_start = find_closest(voltages,vg_start)
				index = find_closest(voltages,vt_sd_medio)
				slope, intercept, _, _, _ = stats.linregress(voltages[index_start:index],np.log10(currents[index_start:index]))
			else:
				index = find_closest(voltages,vt_sd_medio)
				slope, intercept, _, _, _ = stats.linregress(voltages[0:index],np.log10(currents[0:index]))

			# print(1000/slope)
			dummy=np.linspace(0.0,voltages[center],1000)
			yss=dummy*slope+intercept
			plt.plot(dummy,yss,'-',color='r', label='SS')
			plt.ylabel(r'$\mathrm{I_{D}} $',labelpad=20)
			plt.xlabel(r'$\mathrm{V_{G}} $',labelpad=20)
			plt.title('Subthreshold Slope [mV/dec] ->Curve %s' % str(i))
			plt.tight_layout()

			if(save_plot is not None):
				try:
					checkPath(save_plot)
					plt.savefig(save_plot+'SS_{0}.pdf'.format(str(i)), bbox_inches='tight', format='pdf', dpi=1200 )
				except NameError:
					print('No path to save the plots has been chosen')
			else:
				print('No path has been defined\nTrying to use GUI backend...')
				plt.show()
			plt.close()

		if(fom == 'dibl'):
			vth_index_high = find_closest(curve_high[:,0],vth_high)
			corriente_high = curve_low[vth_index_high,1]
			vth_index_low = find_closest(curve_low[:,0],vth_low)
			corriente_low = curve_low[vth_index_low,1]

			plt.close()
			fig, ax1 = plt.subplots()
			ax1.set_yscale("log")
			ax1.plot(curve_low[:,0], curve_low[:,1], '-.', color='r')
			ax1.plot(curve_high[:,0], curve_high[:,1], '-.', color='b')
			ax1.axvline(x=vth_low, color='r', label ='Vt low')
			ax1.axhline(y=corriente_low, color='k')
			ax1.axvline(x=vth_high, color='b', label ='Vt high')
			ax1.set_xlabel(r'$\mathrm{V_{G}} $',labelpad=20)
			ax1.set_ylabel(r'$\mathrm{I_{D}} $',labelpad=20)
			plt.title('DIBL ->Curve %s' % str(i))
			ax1.legend(loc='lower right')
			plt.tight_layout()

			if(save_plot is not None):
				try:
					checkPath(save_plot)
					plt.savefig(save_plot+'DIBL_{0}.pdf'.format(str(i)), bbox_inches='tight', format='pdf', dpi=1200 )
				except NameError:
					print('No path to save the plots has been chosen')
			else:
				print('No path has been defined\nTrying to use GUI backend...')
				plt.show()
			plt.close()

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
