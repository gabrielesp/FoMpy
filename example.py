import fompy

# ------------------------------------------------------------------------------------------#
# BASIC COMMANDS (Example: Vth extraction)

# path_file_JCJB = './data/sim_FinFET_vd_high/'
# fds = fompy.dataset(path_file_JCJB, parser=fompy.JCJB)
#fds.print_parameters()

# vth_array = fompy.extract(fds, fom = 'vth')
# print(vth_array)

# fompy.savetotxt('./results_vth.txt', 'vth', vth_array)

# fompy.plot(fds, fom = 'vth', save_plot='./vth_plots/sd/')


# ------------------------------------------------------------------------------------------#
# EXCLUDING CURVES FROM THE IMPORT

# path_file_JCJB = './data/sim_FinFET_vd_high/'
# fds = fompy.dataset(path_file_JCJB, parser=fompy.JCJB, exclude=[5,6])
# fds.print_parameters()
# fds = fompy.dataset(path_file_JCJB, parser=fompy.JCJB, interval=[0,4])
# fds.print_parameters()

# ------------------------------------------------------------------------------------------#
# DATA CONDITIONING

# path_file_JCJB = './data/sim_FinFET_vd_high/'
# fds = fompy.dataset(path_file_JCJB, parser=fompy.JCJB)
# norm_value = 35.8/10**9
# fompy.normalize(fds, norm_value)
# print(fds.dataset)

# fompy.filter(fds, theta_crit = 0.5, show_theta=True)

# ------------------------------------------------------------------------------------------#
# DIFFERENT EXTRACTION PARSERS 

# import os
# path = './data/default'
# fds = fompy.dataset(path, parser=fompy.file)
# print(fds.dataset)
# print('Number of simulations loaded into the dataset',fds.n_sims)

# import numpy as np
# #Here the arrays are defined
# arr1 =np.array([[0.00e+00, 1.00e-09],[1.00e-01, 2.20e-08],[2.00e-01, 3.20e-07],[3.00e-01, 2.74e-06],[4.00e-01, 9.90e-06],[5.00e-01, 2.20e-05],[6.00e-01, 3.22e-05],[7.00e-01, 4.16e-05],[8.00e-01, 5.23e-05],[9.00e-01, 6.04e-05],[1.00e+00, 6.60e-05]])
# arr2 =np.array([[0.00e+00, 1.00e-09],[1.00e-01, 2.15e-08],[2.00e-01, 3.18e-07],[3.00e-01, 2.72e-06],[4.00e-01, 9.85e-06],[5.00e-01, 2.12e-05],[6.00e-01, 3.16e-05],[7.00e-01, 4.10e-05],[8.00e-01, 5.46e-05],[9.00e-01, 6.15e-05],[1.00e+00, 6.57e-05]])
# arrays = np.stack((arr1, arr2)) #Here the arrays are put together
# fds = fompy.dataset(arr = arrays, parser=fompy.array)
# print(fds.dataset)
# print('Number of simulations loaded into the dataset',fds.n_sims)


#path_file_JCJB = './data/sim_FinFET_vd_high/'
#fds = fompy.dataset(path_file_JCJB, parser=fompy.JCJB)
#fds.print_parameters()
#print(fds.dataset)

# path_file_mc = './data/mc_data/'
# fds1 = fompy.dataset(path_file_mc, parser=fompy.MC)
# fds1.print_parameters()
# print(fds1.dataset)

# ------------------------------------------------------------------------------------------#
# SAVING EXTRACTION TO FILES

# path_to_save='./sim_data.dat'
# path_file_JCJB = './data/sim_FinFET_vd_high/'
# fds = fompy.dataset(path_file_JCJB, parser=fompy.JCJB, save_to_file = path_to_save)

# vth_array = fompy.extract(fds, fom = 'vth')
# fompy.savetotxt('./results_vth.txt', 'vth', vth_array)


# ------------------------------------------------------------------------------------------#
# VTH EXTRACTION & PLOT FULL LIST OF ARGUMENTS

# path_file_JCJB = './data/simulations/'
# fds = fompy.dataset(path_file_JCJB, parser=fompy.JCJB)
# vth_array = fompy.extract(fds, fom = 'vth')
# print(vth_array)
# fompy.plot(fds, fom = 'vth', save_plot='./vth_plots/sd/')
# fompy.plot(fds, fom = 'vth', backend='TkAgg')
# fompy.plot(fds, fom = 'vth')

# path_file_var = './data/simulations/'
# fds_var = fompy.dataset(path_file_var, parser=fompy.JCJB)
# vth_array = fompy.extract(fds_var, fom = 'vth')

# import numpy as np
# print('STDEV:',np.std(vth_array))
# print('MEAN:',np.mean(vth_array))

# path_file_JCJB = './data/sim_FinFET_vd_high/'
# fds = fompy.dataset(path_file_JCJB, parser=fompy.JCJB)
# vth_array_cc = fompy.extract(fds, fom = 'vth', method = 'CC', cc_criteria=1.5e-6)
# print(vth_array_cc)
# fompy.plot(fds, fom = 'vth',method = 'CC',cc_criteria = 1.5e-6, save_plot='./vth_plots/cc/')

# path_file_JCJB = './data/sim_FinFET_vd_high/'
# fds = fompy.dataset(path_file_JCJB, parser=fompy.JCJB)
# vth_array_td = fompy.extract(fds, fom = 'vth', method = 'TD')
# print(vth_array_td)
# fompy.plot(fds, fom = 'vth',method = 'TD', save_plot='./vth_plots/td/')

# path_file_JCJB = './data/simulations/'
# fds = fompy.dataset(path_file_JCJB, parser=fompy.JCJB)
# vth_array_le = fompy.extract(fds, fom = 'vth', method = 'LE')
# print(vth_array_le)
# fompy.plot(fds, fom = 'vth',method = 'LE', save_plot='./vth_plots/le/')

# ------------------------------------------------------------------------------------------#
# IOFF EXTRACTION & PLOT FULL LIST OF ARGUMENTS

# path_file_JCJB = './data/sim_FinFET_vd_high/'
# fds = fompy.dataset(path_file_JCJB, parser=fompy.JCJB)
# ioff_array = fompy.extract(fds, fom = 'ioff', vg_ext = 0.2)
# print(ioff_array)
# fompy.savetotxt('./results_ioff.txt', 'ioff', ioff_array)


# path_file_JCJB = './data/sim_FinFET_vd_high/'
# fds = fompy.dataset(path_file_JCJB, parser=fompy.JCJB)
# ioff_array = fompy.extract(fds, fom = 'ioff', vg_ext = 0.2)
# print(ioff_array)

# fompy.plot(fds, fom = 'ioff', vg_ext = 0.2, save_plot='./ioff_plots/')


# ------------------------------------------------------------------------------------------#
# ION EXTRACTION & PLOT FULL LIST OF ARGUMENTS

# path_file_JCJB = './data/sim_FinFET_vd_high/'
# fds = fompy.dataset(path_file_JCJB, parser=fompy.JCJB)
# ion_array_vg_defined = fompy.extract(fds, fom = 'ion',vg_ext = 0.7)
# print(ion_array_vg_defined)

# ion_array_default_SD = fompy.extract(fds, fom = 'ion')
# print(ion_array_default_SD)

# ion_array_LE = fompy.extract(fds, fom = 'ion',method = 'LE')
# print(ion_array_LE)

# ion_array_TD = fompy.extract(fds, fom = 'ion',method = 'TD')
# print(ion_array_TD)

# fompy.plot(fds, fom = 'ion', vg_ext = 0.5)
# fompy.plot(fds, fom = 'ion')
# fompy.plot(fds, fom = 'ion', method='TD')

# ------------------------------------------------------------------------------------------#
# SS EXTRACTION & PLOT FULL LIST OF ARGUMENTS

# path_file_JCJB = './data/sim_FinFET_vd_high/'
# fds = fompy.dataset(path_file_JCJB, parser=fompy.JCJB)

# ss_array = fompy.extract(fds, fom = 'ss', method = 'LE')
# print(ss_array)
# fompy.plot(fds, fom = 'ss', method='LE')

# ss_array = fompy.extract(fds, fom = 'ss', vg_start = 0.05)
# print(ss_array)
# fompy.plot(fds, fom = 'ss', vg_start = 0.05)

# ss_array = fompy.extract(fds, fom = 'ss', vg_start = 0.05, vg_end = 0.2)
# print(ss_array)
# fompy.plot(fds, fom = 'ss', vg_start = 0.05, vg_end = 0.2)

# ss_array = fompy.extract(fds, fom = 'ss', vg_end = 0.2)
# print(ss_array)
# fompy.plot(fds, fom = 'ss',vg_end = 0.2, save_plot='./ss/')

# ------------------------------------------------------------------------------------------#
# DIBL EXTRACTION & PLOT FULL LIST OF ARGUMENTS

#path_file_JCJB = './data/sim_FinFET_vd_high/'
#path_file_low = './data/sim_FinFET_vd_low/'


#fds_hdb = fompy.dataset(path_file_JCJB, parser=fompy.JCJB)
#fds_ldb = fompy.dataset(path_file_low, parser=fompy.JCJB)
#fds_hdb.drain_bias_value = 0.7
#fds_ldb.drain_bias_value = 0.05

#dibl_array = fompy.extract(fds_hdb, fds_ldb, fom = 'dibl', method = 'SD')
#print(dibl_array)

#fompy.plot(fds_hdb, fds_ldb, fom = 'dibl', save_plot='./dibl/')


# ------------------------------------------------------------------------------------------#
# FULL EXTRACTION & PLOT FULL LIST OF ARGUMENTS

# path_file_JCJB = './data/sim_FinFET_vd_high/'
# path_file_low = './data/sim_FinFET_vd_low/'

# fds_hdb = fompy.dataset(path_file_JCJB, parser=fompy.JCJB)
# fds_hdb.print_parameters()
# fds_ldb = fompy.dataset(path_file_low, parser=fompy.JCJB)
# fds_ldb.print_parameters()

# vth, ioff, ion, ss= fompy.extract(fds_hdb, print_fom=True)

# fds1_parameters, fds2_parameters= fompy.extract(fds_hdb, fds_ldb, print_fom=True)


# ------------------------------------------------------------------------------------------#
# GENERATING DIFFERENT PLOTS

# path_file_var = './data/simulations/'
# fds_var = fompy.dataset(path_file_var, parser=fompy.JCJB)
# vth_array = fompy.extract(fds_var, fom = 'vth')

# fompy.plot(fds, plot_type='hist', parameter=vth_array, bins=10,cont_parameter= 0.38, save_plot='./variability/')
# fompy.plot(fds, plot_type='qq', parameter=vth_array, save_plot='./variability/')

# path_file_var = './data/simulations/'
# fds_var = fompy.dataset(path_file_var, parser=fompy.JCJB)
# fompy.plot(fds, plot_type='varplot', save_plot='./variability/')


# path_file_JCJB = './data/sim_FinFET_vd_high/'
# path_file_low = './data/sim_FinFET_vd_low/'

# fds_hdb = fompy.dataset(path_file_JCJB, parser=fompy.JCJB)
# fds_hdb.print_parameters()
# fds_ldb = fompy.dataset(path_file_low, parser=fompy.JCJB)
# fds_ldb.print_parameters()

# norm_value = 35.8/10**9
# fompy.normalize(fds_hdb, norm_value)
# fompy.normalize(fds_ldb, norm_value)

# fompy.plot(fds_hdb,fds_ldb, plot_type='calib')
