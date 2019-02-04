
from abc import ABCMeta, abstractmethod
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from enum import Enum
import os
import probscale

class plotStrategy(metaclass=ABCMeta):
	''
	def __init__(self, **kwargs):
		self.type = None

	@abstractmethod
	def hist(self):
		pass   


class plotter(plotStrategy):

	def hist(self, bines = None, parameter = None, save_to_file = None):

		if(len(parameter)<2):
			raise Exception('Not enough data to plot a histogram')

		plt.close()
		xmin = np.min(parameter)
		xmax = np.max(parameter)
		x_range = np.linspace(xmin*0.8, xmax*1.2, len(parameter))
		m, s = stats.norm.fit(parameter) # get mean and standard deviation
		pdf_g = stats.norm.pdf(x_range, m, s) # now get theoretical values in our interval
		weights = np.ones_like(parameter)/float(len(parameter))
		# plt.plot(x_range, pdf_g, color = 'r', label='Normal Fit') # plot it
		w, bins, patches = plt.hist(parameter,bines,weights=weights, facecolor='green')
		plt.axvline(np.mean(parameter), color='r', linestyle='dashed', linewidth=2)
		plt.title('Distribution parameters:'+'$\mu='+str(("%0.3f"%m))+', \sigma='+str(("%0.3f"%s))+'$\n'  )
		plt.xlabel('Parameter')
		plt.ylabel('Frecuency')
		# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		if(save_to_file == None):
			plt.show()
		else:
			checkPath(save_to_file)
			plt.savefig(save_to_file+'/hist.pdf', bbox_inches='tight', format='pdf', dpi=1200 )
		plt.close()

		return

	def qq(self, parameter = None, method = None, save_to_file = None):

		if(len(parameter)<2):
			raise Exception('Not enough data for a QQ-plot')

		plt.close()
		fig, ax = plt.subplots()
		slope, r = probscale.probplot(parameter, ax=ax, plottype='qq', bestfit=True, estimate_ci=False, return_best_fit_results=True,problabel='Standard Normal Quantiles',
		datalabel=r'Parameter',
		scatter_kws={'markersize': 10, 'linestyle': 'none', 'label': method})
		slope, intercept, rcoef, prob, sterrest = stats.linregress(r['x'],r['y'])
		plt.text(1.5,max(r['y']), "$R^2_{"+str(method)+"}=%1.2f$" % rcoef**2)
		plt.legend(loc='lower right', fontsize='medium')
		plt.title('QQ-plot')
		if(save_to_file == None):
			plt.show()
		else:
			checkPath(save_to_file)
			plt.savefig(save_to_file+'/qq.pdf', bbox_inches='tight', format='pdf', dpi=1200 )		
		plt.close()
		return

	def varplot(self,fds, save_to_file = None):

		plt.close()
		fig, ax1 = plt.subplots()
		ax2 = ax1.twinx()	

		for i in range(len(fds.dataset)):
			ax1.plot(fds.dataset[i][:,0][:],fds.dataset[i][:,1][:],lw =0.4,color='red', label='Log')
			ax2.plot(fds.dataset[i][:,0][:],fds.dataset[i][:,1][:],lw =0.4,color='green', label='Lin')
			ax1.set_yscale('log')
			ax1.set_ylabel('$I_{D} [\mu A/\mu m]$')
			ax1.set_xlabel('$V_{G} [V]$')
			ax1.set_title('Variability plot')
			handles2, labels2 = ax1.get_legend_handles_labels()
			labels2, ids2 = np.unique(labels2, return_index=True)
			ax1.legend(handles2, labels2, loc='best')
			ax2.set_yscale('linear')
			fig.tight_layout()
			handles, labels = ax2.get_legend_handles_labels()
			labels, ids = np.unique(labels, return_index=True)
			ax2.legend(handles, labels, loc='lower right')
		if(save_to_file == None):
			plt.show()
		else:
			checkPath(save_to_file)
			plt.savefig(save_to_file+'/varplot.pdf', bbox_inches='tight', format='pdf', dpi=1200 )	
		plt.close()
		return

	def calib(self,fds, save_to_file = None):

		plt.close()
		fig, ax1 = plt.subplots()
		ax2 = ax1.twinx()	

		for i in range(len(fds.dataset)):
			ax1.plot(fds.dataset[i][:,0][:],fds.dataset[i][:,1][:],lw =0.4,color='red', label='Log')
			ax2.plot(fds.dataset[i][:,0][:],fds.dataset[i][:,1][:],lw =0.4,color='green', label='Lin')
			ax1.set_yscale('log')
			ax1.set_ylabel('$I_{D} [\mu A/\mu m]$')
			ax1.set_xlabel('$V_{G} [V]$')
			ax1.set_title('Variability plot')
			handles2, labels2 = ax1.get_legend_handles_labels()
			labels2, ids2 = np.unique(labels2, return_index=True)
			ax1.legend(handles2, labels2, loc='best')
			ax2.set_yscale('linear')
			fig.tight_layout()
			handles, labels = ax2.get_legend_handles_labels()
			labels, ids = np.unique(labels, return_index=True)
			ax2.legend(handles, labels, loc='lower right')
		if(save_to_file == None):
			plt.show()
		else:
			checkPath(save_to_file)
			plt.savefig(save_to_file+'/varplot.pdf', bbox_inches='tight', format='pdf', dpi=1200 )	
		plt.close()
		return

	def fom(self, fom = None, method = None, parameter = None, save_to_file = None):

		if(len(parameter)<2):
			raise Exception('Not enough data for a variability curve')
		
		if(fom == None):
			raise Exception('The figure of merit to plot wasn\'t defined!')

		if(fom == 'vth'):
			fig, ax1 = plt.subplots()
			ax2 = ax1.twinx()

			ax1.plot(x,y,lw =0.4,color='red', label='Log HDB')
			ax2.plot(x,y,lw =0.4,color='red', label='Lin HDB')

			ax1.plot(x,y,lw =0.4,color='red', label='Log LDB')
			ax2.plot(x,y,lw =0.4,color='red', label='Lin LDB')
			
		if(fom == 'ioff'):
			pass
		if(fom == 'ion'):
			pass
		if(fom == 'ratio'):
			pass
		if(fom == 'power'):
			pass
		if(fom == 'dibl'):		
			pass					


		return


### checkPath :FUNCTION THAT CHECKS IF A FOLDER EXISTS AND IF NOT CREATES IT
### INPUT:
### cwd: directory to be checked/created
### OUTPUT:
### None

def checkPath(cwd):
	try:
		os.makedirs(cwd)
	except OSError:
		if not os.path.isdir(cwd):
			raise
