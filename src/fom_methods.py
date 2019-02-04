from src import dataset as dt
from abc import ABCMeta, abstractmethod
from enum import Enum
import importlib

	#        if (parameter.lower()== "ioff"):
	#           fom_output = iv_curve_data[0,1]
	#           if(plot==True):
	#                  plt.clf()
	#                  fig, ax1 = plt.subplots()
	#                  ax1.plot(voltages,np.log10(currents),color='orangered',label ='Simulated Data')
	#                  ax1.axhline(y=np.log10(float(fom_output)), color='blue', linestyle='-.', label ='CC criteria')
	#                  ax1.set_ylabel('$\mathrm{I_{D}} $',labelpad=20)
	#                  ax1.set_xlabel('$\mathrm{V_{G}} $',labelpad=20)
	#                  ax1.legend(loc='lower right',fontsize =25)

	#                  if(saveplot==True):
	#                         checkPath('./plots/ioff_output/')
	#                         plt.savefig(plotpath + 'ioff_'+"%03d"%i+'.pdf', bbox_inches='tight', format='pdf', dpi=1200 )
	#                  else:
	#                         plt.show()

	#                  plt.clf()

class _extractor(metaclass=ABCMeta):
	
	@abstractmethod
	def extract(self):
		pass
	@abstractmethod
	def save_results_to_file(self):
		pass
	@abstractmethod
	def plot(self):
		pass
	@abstractmethod
	def save_plot_to_file(self):
		pass
	@abstractmethod
	def test_extract(self):
		pass	

class ioff(_extractor):

	def extract(self,fds,method=None):

		if(method == dt.ioffMethod.default) or(method == None):
			length = len(fds.dataset)
			ioff_out = []
			for i in range(length):
				ioff_out.append(fds.dataset[i][0][1])
			
		return ioff_out

	def save_results_to_file(self):
		pass
	def plot(self):
		pass
	def save_plot_to_file(self):
		pass
	def test_extract(self):
		pass			


class parameter(Enum):
	'This class implements al figures of merit that can be extracted'
	vth, ion, ioff, ratio, power, dibl = range(6)

class vthMethod(Enum):
	'This class implements all the threshold voltage extraction methods'
	SD, CC, TD, LE = range(4)


class ioffMethod(Enum):
	'This class implements all the ioff extraction methods'
	default, crit= range(2)	




# def extraction_ioff( method=None, save_results_to_file = None, plot = None, save_plot_to_file = None, test_extract = None):
# 	if(method == dt.ioffMethod.default):
# 		length = len(self.dataset)
# 		for i in range(length):
# 			ioff_out = self.dataset[i][0]
		
# 	return ioff_out

# 
# def extraction_fom(FompyDataSet, parameter = None, ext_method = None, cc_crit = None, plot = None, plot_to_file = None):
# 	fom = []
# 	x = FompyDataset.dataset[0]
# 	y = FompyDataset.dataset[1]
# 	if(ext_method==vthMethod.SD):
		# x_interp,y_interp = interpol(voltages, currents, 1000,0)
		#               iv_curve_data = np.column_stack((x_interp, y_interp))
		#               if(bias.lower()=="high"):
		#                      iv_curve_data[:,1] = np.power( iv_curve_data[:,1],0.5)
		#               d1 = get_diff(iv_curve_data, order = 1,diff = 'central')
		#               d2 = get_diff(iv_curve_data, order = 2,diff = 'central')

		#               lower_bound=find_closest(d2[:,0],10*d2[2,0])
		#               higher_bound=find_closest(d2[2:,0],0.95*d2[-1,0])
		#               warning_interval_limit_low = int(np.round(lower_bound*1.5))
		#               warning_interval_limit_high = int(np.round(higher_bound*0.95))
		#               vt_temp = np.argmax(d2[lower_bound:higher_bound,1])

		#               if (d2[vt_temp+lower_bound, 0] <= d2[warning_interval_limit_low, 0]):
		#                      print('$V_{T}$ value outside of the confidence interval for simulation %i', i)
		#               if (d2[vt_temp+lower_bound, 0] >= d2[warning_interval_limit_high, 0]):
		#                      print('$V_{T}$ value outside of the confidence interval for simulation %i', i)
		#               vt_SD_float = d2[vt_temp+lower_bound,0]
		#               index_lower_bound=find_closest(d2[2:,0],0.05)
		#               index_higher_bound=find_closest(d2[index_lower_bound:,1],0)
		#               try:
		#                      vt_temp = np.argmax(d2[index_lower_bound:index_higher_bound+index_lower_bound,1])
		#                      vt_SD_float = d2[vt_temp+index_lower_bound,0]
		#               except ValueError:
		#                      print('Multiple indexes')
		#               parameter_output_vt = float(("%0.3f"%vt_SD_float))

		#               try:
		#                      v_shifted = parameter_output_vt + float(bias_voltage)
		#                      index_vt_on_SD=find_closest(iv_curve_data[:,0],v_shifted)
		#                      parameter_output_ion=iv_curve_data[index_vt_on_SD,1]
		#               except (NameError,TypeError) as e:
		#                      print('No bias voltage defined')
		#                      bias_voltage = 0
	# elif(ext_method==vthMethod.CC):
	# 	raise NotImplementedError
	# elif(ext_method==vthMethod.LE):
	# 	raise NotImplementedError
	# elif(ext_method==vthMethod.TD):
	# 	raise NotImplementedError
	# else:
	# 	raise Exception('No Vth extraction method has been defined!')
	# return fom

	# def extraction_vt(self,i,x,y,bias,method = None,bias_voltage=None,cc_crit=None,plotpath=None,plot=False,saveplot=False):

	#        i = i +1
	#        #print(i)
	#        voltages = x
	#        currents = y

	#        if(method.lower()=="sd"):

	#               x_interp,y_interp = interpol(voltages, currents, 1000,0)
	#               iv_curve_data = np.column_stack((x_interp, y_interp))
	#               if(bias.lower()=="high"):
	#                      iv_curve_data[:,1] = np.power( iv_curve_data[:,1],0.5)
	#               d1 = get_diff(iv_curve_data, order = 1,diff = 'central')
	#               d2 = get_diff(iv_curve_data, order = 2,diff = 'central')

	#               lower_bound=find_closest(d2[:,0],10*d2[2,0])
	#               higher_bound=find_closest(d2[2:,0],0.95*d2[-1,0])
	#               warning_interval_limit_low = int(np.round(lower_bound*1.5))
	#               warning_interval_limit_high = int(np.round(higher_bound*0.95))
	#               vt_temp = np.argmax(d2[lower_bound:higher_bound,1])

	#               if (d2[vt_temp+lower_bound, 0] <= d2[warning_interval_limit_low, 0]):
	#                      print('$V_{T}$ value outside of the confidence interval for simulation %i', i)
	#               if (d2[vt_temp+lower_bound, 0] >= d2[warning_interval_limit_high, 0]):
	#                      print('$V_{T}$ value outside of the confidence interval for simulation %i', i)
	#               vt_SD_float = d2[vt_temp+lower_bound,0]
	#               index_lower_bound=find_closest(d2[2:,0],0.05)
	#               index_higher_bound=find_closest(d2[index_lower_bound:,1],0)
	#               try:
	#                      vt_temp = np.argmax(d2[index_lower_bound:index_higher_bound+index_lower_bound,1])
	#                      vt_SD_float = d2[vt_temp+index_lower_bound,0]
	#               except ValueError:
	#                      print('Multiple indexes')
	#               parameter_output_vt = float(("%0.3f"%vt_SD_float))

	#               try:
	#                      v_shifted = parameter_output_vt + float(bias_voltage)
	#                      index_vt_on_SD=find_closest(iv_curve_data[:,0],v_shifted)
	#                      parameter_output_ion=iv_curve_data[index_vt_on_SD,1]
	#               except (NameError,TypeError) as e:
	#                      print('No bias voltage defined')
	#                      bias_voltage = 0
	#               if(plot==True):
	#                      plt.clf()
	#                      fig,ax1 = plt.subplots()
	#                      ax1.scatter(voltages,currents,color='orangered',label ='Simulated Data')
	#                      plot_1 = ax1.plot(voltages,currents,color='orangered',label ='Simulated Data')
	#                      ax1.set_ylim(np.min(currents), np.max(currents))
	#                      ax1.set_ylabel('$\mathrm{I_{D}} $',labelpad=20)
	#                      ax1.set_xlabel('$\mathrm{V_{G}} $',labelpad=20)
	#                      textstr = '$V_{T}=%.3f$'%(parameter_output_vt)
	#                      props = dict(boxstyle="square,pad=0.1", fc='w', alpha=1,ec='k', lw =2)
	#                      ax2 = ax1.twinx()
	#                      ax2.text(0.05, 0.55, textstr, transform=ax1.transAxes, fontsize=14,  va='top', bbox=props,family = 'serif', style = 'normal')
	#                      plot_2 = ax2.plot(d2[lower_bound:higher_bound,0],d2[lower_bound:higher_bound,1], label ='SD ($d^{2}I_{D}/dV_{G}^{2}$)', color = 'dodgerblue')
	#                      ax2.scatter(d2[lower_bound:higher_bound,0],d2[lower_bound:higher_bound,1], label ='SD ($d^{2}I_{D}/dV_{G}^{2}$)',s=0.1, color = 'dodgerblue')
	#                      ax2.set_ylim(np.min(d2[lower_bound:higher_bound,1]), np.max(d2[lower_bound:higher_bound,1]))
	#                      plots = plot_1+plot_2
	#                      labs = [l.get_label() for l in plots]
	#                      ax1.legend(plots, labs, loc=1,bbox_to_anchor=(1.017,0.225))
	#                      ax2.yaxis.set_visible(False)

	#                      if(saveplot==True):
	#                             checkPath(str(plotpath))
	#                             plt.savefig(plotpath + 'VT_SD_PLOT_'+"%03d"%i+'.pdf', bbox_inches='tight', format='pdf', dpi=1200 )
	#                      else:
	#                             plt.show()
	#                      plt.clf()
	#               if(bias.lower()=="high"):
	#                      iv_curve_data[:,1] = np.power( iv_curve_data[:,1],2)

	#        if(method.lower()=="cc"):
	#               x_interp,y_interp = interpol(voltages, currents, 1000,0)
	#               iv_curve_data = np.column_stack((x_interp, y_interp))
	#               try:
	#                      vt_index_cc=find_closest(iv_curve_data[:,1],float(cc_crit))
	#                      vt_CC_float =iv_curve_data[vt_index_cc,0]
	#                      parameter_output_vt = float(("%0.3f"%vt_CC_float))
	#               except (NameError,TypeError) as e:
	#                      print('No CC criteria defined')
	#                      parameter_output_vt = 0
	#               try:
	#                      v_shifted = parameter_output_vt + float(bias_voltage)
	#                      index_vt_on_cc=find_closest(iv_curve_data[:,0],v_shifted)
	#                      parameter_output_ion=iv_curve_data[index_vt_on_cc,1]
	#               except (NameError,TypeError) as e:
	#                      print('No bias voltage defined')
	#                      bias_voltage = 0

	#               if(plot==True):
	#                      plt.clf()
	#                      fig, ax1 = plt.subplots()

	#                      ax1.plot(iv_curve_data[:,0],np.log10(iv_curve_data[:,1]),color='k',label ='Simulated Data')

	#                      ax1.set_ylabel('$\mathrm{I_{D}} $',labelpad=20)
	#                      ax1.set_xlabel('$\mathrm{V_{G}} $',labelpad=20)
	#                      textstr = '$V_{T}=%.3f$'%(parameter_output_vt)
	#                      props = dict(boxstyle="square,pad=0.1", fc='w', alpha=1,ec='k', lw =2)
	#                      ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=14,  va='top', bbox=props,family = 'serif', style = 'normal')
	#                      # ax1.axhline(y=float(cc_crit)*(10**9/35.8), color='blue', linestyle='-.', label ='CC criteria')
	#                      ax1.axhline(y=np.log10(float(cc_crit)), color='blue', linestyle='-.', label ='CC criteria')
	#                      ax1.legend(loc='lower right',fontsize =25)

	#                      if(saveplot==True):
	#                             checkPath(str(plotpath))
	#                             plt.savefig(plotpath + 'VT_CC_PLOT_'+"%03d"%i+'.pdf', bbox_inches='tight', format='pdf', dpi=1200 )
	#                      else:
	#                             plt.show()
	#                      plt.clf()

	#        if(method.lower()=="td"):
	#               x_interp,y_interp = interpol(voltages, currents, 1000,0)
	#               iv_curve_data = np.column_stack((x_interp, y_interp))
	#               if(bias.lower()=="high"):
	#                      iv_curve_data[:,1] = np.power( iv_curve_data[:,1],0.5)

	#               d3 = get_diff(iv_curve_data, order = 3,diff = 'central')
 	# 				# checkPath('.././out/third_deriv/')
 	# 				# np.savetxt('.././out/third_deriv/third_deriv_'+str(i)+'.csv', d3, delimiter='\t')
	#               index_lower_bound=find_closest(d3[2:,0],0.1)
	#               index_higher_bound=find_closest(d3[index_lower_bound:,1],0)
	#               try:
	#                      vt_temp = np.argmax(d3[index_lower_bound:index_higher_bound+index_lower_bound,1])
	#                      vt_TD_float = d3[vt_temp+index_lower_bound,0]
	#               except ValueError:
	#                      print('Multiple indexes')
	#               parameter_output_vt = float(("%0.3f"%vt_TD_float))


	#               try:
	#                      v_shifted = parameter_output_vt + float(bias_voltage)
	#                      index_vt_on_td=find_closest(iv_curve_data[:,0],v_shifted)
	#                      parameter_output_ion=iv_curve_data[index_vt_on_td,1]
	#               except (NameError,TypeError) as e:
	#                      print('No bias voltage defined')
	#                      bias_voltage = 0

	#               if(plot==True):
	#                      plt.clf()
	#                      fig, ax1 = plt.subplots()
	#                      ax2 = ax1.twinx()
	#                      plot_1 = ax1.plot(voltages,currents,color='orangered',label ='Simulated Data')
	#                      ax1.set_ylabel('$\mathrm{I_{D}} $',labelpad=20)
	#                      ax1.set_xlabel('$\mathrm{V_{G}} $',labelpad=20)
	#                      textstr = '$V_{T}=%.3f$'%(parameter_output_vt)
	#                      props = dict(boxstyle="square,pad=0.1", fc='w', alpha=1,ec='k', lw =2)
	#                      ax2.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=14,  va='top', bbox=props,family = 'serif', style = 'normal')
	#                      plot_2 = ax2.plot(d3[index_lower_bound:index_higher_bound+index_lower_bound,0],d3[index_lower_bound:index_higher_bound+index_lower_bound,1], label ='SD ($d^{2}I_{D}/dV_{G}^{2}$)', color = 'dodgerblue')
	#                      plots = plot_1+plot_2
	#                      labs = [l.get_label() for l in plots]
	#                      ax1.legend(plots, labs, loc=1,bbox_to_anchor=(1.017,0.225))
	#                      ax2.yaxis.set_visible(False)

	#                      if(saveplot==True):
	#                             checkPath(str(plotpath))
	#                             plt.savefig(plotpath + 'VT_TD_PLOT_'+"%03d"%i+'.pdf', bbox_inches='tight', format='pdf', dpi=1200 )
	#                      else:
	#                             plt.show()

	#                      plt.clf()

	#               if(bias.lower()=="high"):
	#                      iv_curve_data[:,1] = np.power( iv_curve_data[:,1],2)

	#        if(method.lower()=="le"):
	#               def line(x, A, B):
	#                      return A*x + B
	#               x_interp,y_interp = interpol(voltages, currents, 1000,0)
	#               iv_curve_data = np.column_stack((x_interp, y_interp))
	#               if(bias.lower()=="high"):
	#                      iv_curve_data[:,1] = np.power( iv_curve_data[:,1],0.5)

	#               d1 = get_diff(iv_curve_data, order = 1,diff = 'central')
	#               d2 = get_diff(iv_curve_data, order = 2,diff = 'central')

	#               index_lower_bound=find_closest(d1[2:,0],0.05)
	#               index_higher_bound=find_closest(d1[index_lower_bound:,1],1)

	#               try:
	#                      vt_temp1 = np.argmax(d1[index_lower_bound:index_higher_bound,1])
	#                      index_d1max = vt_temp1+index_lower_bound
	#                      vt_temp_float = d1[index_d1max,0]
	#               except ValueError:
	#                      print('Multiple indexes')

	#               vt_temp = float(("%0.3f"%vt_temp_float))
	#               corriente_vt_temp =iv_curve_data[index_d1max,1]
	#               lim1 = int(index_d1max)
	#               lim2 = int(index_d1max*1.2)
	#               A,B = curve_fit(line, iv_curve_data[lim1:lim2,0],iv_curve_data[lim1:lim2,1])[0]
	#               fit = A*iv_curve_data[:,0]+B
	#               vt_index_le = find_closest(fit,0)

	#               if(bias.lower()=="high"):
	#                      vd_high = float(bias_voltage)
	#                      # vt_LE_float = iv_curve_data[vt_index_le,0]+vd_high/2
	#                      vt_LE_float = iv_curve_data[vt_index_le,0]

	#               if(bias.lower()=="low"):
	#                      vt_LE_float = iv_curve_data[vt_index_le,0]
 	# 					#corriente_spline_le =iv_curve_data[vt_index_le,1]
 	# 					#vt_vals = np.column_stack((vt_LE_float,corriente_spline_le))
 	# 					#checkPath('.././out/i_vt_le')
 	# 					#np.savetxt('.././out/i_vt_le/iv_'+str(i)+'.csv', vt_vals, delimiter='\t')

	#               parameter_output_vt = float(("%0.3f"%vt_LE_float))
	#               try:
	#                      v_shifted = parameter_output_vt + float(bias_voltage)
	#                      index_vt_on_le=find_closest(iv_curve_data[:,0],v_shifted)
	#                      parameter_output_ion=iv_curve_data[index_vt_on_le,1]
	#               except (NameError,TypeError) as e:
	#                      print('No bias voltage defined')
	#                      bias_voltage = 0


	#               if(plot==True):
	#                      plt.clf()
	#                      fig, ax1 = plt.subplots()
	#                      ax2 = ax1.twinx()
	#                      ax1.plot(iv_curve_data[:,0],iv_curve_data[:,1],color='red',label ='Interp Data')
	#                      ax1.plot(iv_curve_data[:,0], A*iv_curve_data[:,0]+B,color='blue',linestyle ='-.', label='Linear Extrapolation')
	#                      ax1.set_ylabel('$\mathrm{I_{D}} $',labelpad=20)
	#                      ax1.set_xlabel('$\mathrm{V_{G}} $',labelpad=20)
	#                      ax1.set_ylim([0, np.max(iv_curve_data[:,1])])
	#                      textstr = '$V_{T}=%.3f$'%(parameter_output_vt)
	#                      props = dict(boxstyle="square,pad=0.1", fc='w', alpha=1,ec='k', lw =2)
	#                      ax1.legend(loc='lower right')
	#                      ax2.yaxis.set_visible(False)

	#                      if(saveplot==True):
	#                             checkPath(str(plotpath))
	#                             plt.savefig(plotpath + 'VT_LE_PLOT_'+"%03d"%i+'.pdf', bbox_inches='tight', format='pdf', dpi=1200 )
	#                      else:
	#                             plt.show()
	#                      plt.clf()
	#               if(bias.lower()=="high"):
	#                      iv_curve_data[:,1] = np.power( iv_curve_data[:,1],2)

	#        return (parameter_output_vt, parameter_output_ion)

	# def extraction_fom(self,i,x,y,parameter,bias,vt_array = None,bias_voltage=None,plotpath=None,plot=False, saveplot=False):

	#        voltages = x
	#        currents = y


	#        if(parameter.lower()== "ss"):
	#           iv_curve_data = np.column_stack((voltages, currents))
	#           x=[iv_curve_data[0,0],iv_curve_data[1,0]]
	#           y=[iv_curve_data[0,1],iv_curve_data[1,1]]
	#           icero=iv_curve_data[0,1]
	#           icerouno=iv_curve_data[1,1]
	#           fom_output=50/(np.log10(icerouno)-np.log10(icero))
	#           if(plot==True):
	#                  plt.clf()
	#                  fig, ax1 = plt.subplots()
	#                  plt.plot(iv_curve_data[:,0], np.log10(iv_curve_data[:,1]), 'o')
	#                  slope, intercept, r_value, p_value, std_err = stats.linregress(iv_curve_data[:3,0],np.log10(iv_curve_data[:3,1]))
	#                  dummy=np.linspace(0.0,0.6,1000)
	#                  ss= np.asarray(fom_output)
	#                  yss=dummy*(1/ss)*1000+intercept
	#                  plt.plot(dummy,yss,'-',color='r', label='SS')
	#                  plt.ylabel('$\mathrm{I_{D}} $',labelpad=20)
	#                  plt.xlabel('$\mathrm{V_{G}} $',labelpad=20)

	#                  if(saveplot==True):
	#                         checkPath('./plots/ss_output/')
	#                         plt.savefig(plotpath + 'ss_'+"%03d"%i+'.pdf', bbox_inches='tight', format='pdf', dpi=1200 )
	#                  else:
	#                         plt.show()

	#                  plt.clf()

	#        #print(i)
	#        x_interp,y_interp = interpol(voltages, currents, 1000,0)
	#        iv_curve_data = np.column_stack((x_interp, y_interp))

	#        if (parameter.lower()== "ioff"):
	#           fom_output = iv_curve_data[0,1]
	#           if(plot==True):
	#                  plt.clf()
	#                  fig, ax1 = plt.subplots()
	#                  ax1.plot(voltages,np.log10(currents),color='orangered',label ='Simulated Data')
	#                  ax1.axhline(y=np.log10(float(fom_output)), color='blue', linestyle='-.', label ='CC criteria')
	#                  ax1.set_ylabel('$\mathrm{I_{D}} $',labelpad=20)
	#                  ax1.set_xlabel('$\mathrm{V_{G}} $',labelpad=20)
	#                  ax1.legend(loc='lower right',fontsize =25)

	#                  if(saveplot==True):
	#                         checkPath('./plots/ioff_output/')
	#                         plt.savefig(plotpath + 'ioff_'+"%03d"%i+'.pdf', bbox_inches='tight', format='pdf', dpi=1200 )
	#                  else:
	#                         plt.show()

	#                  plt.clf()


	#        if (parameter.lower()== "ratio"):

	#           ioff_ratio = iv_curve_data[0,1]
	#           vt = vt_array
	#           v_shifted=bias_voltage+vt
	#           index_vt_on=find_closest(iv_curve_data[:,0],v_shifted)
	#           ion_ratio = iv_curve_data[index_vt_on,1]
	#           fom_output = (ion_ratio/ioff_ratio)

	#           if(plot==True):
	#                  plt.clf()
	#                  fig, ax1 = plt.subplots()
	#                  plot_1 = ax1.plot(voltages,np.log10(currents),color='orangered',label ='Simulated Data')
	#                  ax1.axhline(y=np.log10(float(ioff_ratio)), color='blue', label ='Ioff')
	#                  ax1.axhline(y=np.log10(float(ion_ratio)), color='blue', label ='Ion')
	#                  ax1.set_ylabel('$\mathrm{I_{D}} $',labelpad=20)
	#                  ax1.set_xlabel('$\mathrm{V_{G}} $',labelpad=20)
	#                  ax1.legend(loc='lower right')

	#                  if(saveplot==True):
	#                         checkPath('./plots/ratio_output/')
	#                         plt.savefig(plotpath + 'ratio_'+"%03d"%i+'.pdf', bbox_inches='tight', format='pdf', dpi=1200 )
	#                  else:
	#                         plt.show()

	#                  plt.clf()

	#        if (parameter.lower()== "power"):

	#           ioff_power = iv_curve_data[0,1]
	#           fom_output=ioff_power*float(bias_voltage)

	#        return fom_output



#Function definitions

# def dibl_method(FompyDataset_HDB, FompyDataset_LDB, plot = None, plot_norm = None, norm_factor = None, plot_to_file = None, dibl_plot_path = None):

# 	dibl_results = []

# 	if(str(FompyDataset_HDB.drain_bias)!='drain_bias.High'):
# 		raise Exception('The first argument of the dibl method has to be a FompyDataSet simulated at high drain bias!')
# 	if(str(FompyDataset_LDB.drain_bias)!='drain_bias.Low'):
# 		raise Exception('The second argument of the dibl method has to be a FompyDataSet simulated at low drain bias!')	

	# for i in range(len(FompyDataset_HDB.dataset[0])):
		# x_high = FompyDataset_HDB.dataset
	# a = FompyDataset_HDB.dataset
	# print(a)
		# y_high = FompyDataset_HDB.dataset
		# print(y_high)

		# bias_voltage_high = FompyDataset_HDB.drain_bias
		# x_interp_high,y_interp_high = interpol(voltages_high, currents_high, 1000,0)
		# iv_curve_data_high = np.column_stack((x_interp_high, y_interp_high))

		# x_low = FompyDataset_LDB.dataset[0]
		# y_low = FompyDataset_LDB.dataset[1]
		# bias_voltage_low = FompyDataset_LDB.drain_bias
		# x_interp_low,y_interp_low = interpol(voltages_low, currents_low, 1000,0)
		# iv_curve_data_low = np.column_stack((x_interp_low, y_interp_low))

		# vt_low = FompyDataset_LDB.extraction_fom(parameter=fom.vth, ext_method=vthMethod.CC, cc_criteria=7.7e-7)
		# vt_index_low = find_closest(iv_curve_data_low[:,0],vt_low)
		# corriente_low = iv_curve_data_low[vt_index_low,1]
		# vt_low = iv_curve_data_low[vt_index_low,0]

		# vt_index_high = find_closest(iv_curve_data_high[:,1],corriente_low)
		# vt_high = iv_curve_data_high[vt_index_high,0]

		# print(vt_high)
		# print(vt_low)
		# vd_high=bias_voltage_high
		# vd_low=bias_voltage_low
		# print(vd_high)
		# print(vd_low)
		# dibl=(-(vt_high-vt_low)/(vd_high-vd_low))

		# #Pasamos a mV/V

		# dibl =dibl*1000
		# print(dibl)


		# if(plot!=None):
		# 	plt.clf()

		# 	if(plot_norm!=None):
		# 		if(norm_factor!=None)
		# 			iv_curve_data_low[:,1] = np.multply(iv_curve_data_low[:,1],norm_factor)
		# 		else:
		# 			raise Exception('Please specify the multiplying normalization factor')

		# 	fig, ax1 = plt.subplots()
		# 	ax1.plot(iv_curve_data_low[:,0], iv_curve_data_low[:,1], '-.', color='r')
		# 	ax1.plot(iv_curve_data_high[:,0], iv_curve_data_high[:,1], '-.', color='b')
		# 	ax1.axvline(x=vt_low, color='r', label ='Vth LDB')
		# 	ax1.axhline(y=corriente_low, color='k')
		# 	ax1.axvline(x=vt_high, color='b', label ='Vth HDB')
		# 	ax1.set_ylabel('$\mathrm{I_{DS}} $',labelpad=20)
		# 	ax1.set_xlabel('$\mathrm{V_{G}} $',labelpad=20)
		# 	ax1.legend(loc='best')

		# 	if(plot_to_file!=None):
		# 		checkPath(dibl_plot_path)
		# 		plt.savefig(dibl_plot_path + 'dibl_'+"%03d"%i+'.pdf', bbox_inches='tight', format='pdf', dpi=1200 )
		# 	else:
		# 		plt.show()

		# 	plt.clf()

	# return dibl_results
