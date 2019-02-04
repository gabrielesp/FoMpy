from src.aux_fun import *
from src.dataset import *
from src.fom_methods import *

# from src.conditioning import *
from src.plotter import *
import numpy as np

## IMPORTING THE DATA TO A FOMPY DATASET


####### ATTRIBUTE DEFINITION

## AS AN ARGUMENT IN THE CONSTRUCTOR
fompy_1 = FompyDataset(n_sims=300, norm = 1 , ext_method = vthMethod.SD, drain_bias = drain_bias.High, drain_bias_value=0.05)
# fompy_1.print_parameters()    #THE print_parameters() FUNCTION CAN BE CALLED TO PRINT ALL THE PARAMETERS TO THE TERMINAL

## ALSO THE FompyDataset CAN BE CONSTRUCTED WITHOUT THE PARAMETERS
## AND LATER ON UPDATE THEM
fompy_1.__dict__.update(n_sims=800)
fompy_1.print_parameters()

####### DATA IMPORTING

## THE DAO IS USED TO LOAD A DATASET FROM A FILE
## IF A LIST OF FILES WANTS TO BE LOADED THE USER CAN EITHER PASS
## A DIRECTORY PATH OR USE THE * IN UNIX SYSTEMS

path_file='./sim_FinFET_GER-CL40-RMS1/*/JCJB.dat.0.4'
dao_dataset = daoJCJB()
fompy_2 = dao_dataset.load(path_file)

# print(len(fompy_2.dataset))
# print(fompy_2.dataset[2][:])
# print(fompy_2.n_sims)

# The DAO also implements a save option

path_to_save=('./sim_data.dat')
dao_dataset = daoJCJB()
dao_dataset.save(fompy_2, path_to_save)

## This is a DAO wrapper for the FompyDataset

fompy_3 = JCJBtoDataset(path_file)
# print(fompy_3.dataset)
# print(fompy_3.n_sims)

datasettoFile(fompy_3, path_to_save)

####### PARAMETER EXTRACTION

# print(fompy_3.dataset)

temp = ioff()
ioff_fds3 = temp.extract(fompy_3,method=ioffMethod.default)
print(ioff_fds3)

## This is a wrapper for the extraction of the ioff values from the FompyDataset
ioff_fds1 = fompy_3.extraction(fom = parameter.ioff)
# print(ioff_fds3)

####### DATA VISUALIZATION
fig = plotter()
fig.hist(bines=15, parameter=ioff_fds3, save_to_file="./plots_CL_40_RMS1")
fig.qq(parameter=ioff_fds3,method='SD', save_to_file="./plots_CL_40_RMS1")
fig.varplot(fompy_3,save_to_file="./plots_CL_40_RMS1")


