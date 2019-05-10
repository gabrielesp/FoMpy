.. FoMpy documentation master file, created by
   sphinx-quickstart on Fri Feb  8 16:55:44 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to FoMpy's documentation!
=================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

	

=================================
Getting Started
=================================

FoMPy is an effective tool that extracts the main figures of merit (FoM) of a semiconductor's IV curve and provides useful statistical parameters for variability studies. It includes several methods to extract the threshold voltage. 

The source code of FoMpy can be downloaded from https://github.com/gabrielesp/FoMpy and is purely intended to run as a library for Linux systems. Note that the following instalation steps work on Debian-derived distributions: e.g. for Ubuntu 16.04 or later and CentOS (tested). Also, FoMpy has been proven to work in Python 3.5.4.

In the figure below the used can see the basic workflow behind the FoMpy library:


.. figure:: ./figs/simplified_diagram.jpg
    :width: 450px
    :align: center
    :height: 220px
    :alt: alternate text
    :figclass: align-center

    Diagram showing the basic methodology of the FoMpy library. After loading the data into a FoMpy Dataset, using the various tools implemented in the library, the user is able to process and extract important parameters from a given dataset.



Before using FoMpy you need to have installed **pip3** on your system. For Ubuntu, open up a terminal and type::

	sudo apt update

	sudo apt install python3-pip

The use of virtual environments is highly encouraged. The main purpose of Python virtual environments is to create an isolated environment for Python projects so that no modular dependency issues with other projects can appear. In order to use them run the following commands in a terminal::

	#Install virtual environments
	sudo apt install python3-venv 

	#Create and name the environment "venv"
	python3 -m venv .venv

	#Activate the venv
	source .venv/bin/activate

Note that as of this moment you're inside a virtual environment (Notice (.venv) $ in the terminal) with a limited/isolated version of python and therefore you will have to install all the packages you need for that particular project (including the ones you may have installed in the system as they may not be installed in the virtual environment).

**Instalation of FoMpy via pip3**

Run the following command in a terminal::

	pip3 install --extra-index-url https://test.pypi.org/simple/ fompy

and check the library is installed by importing it from a **python3 terminal**::

 	import fompy

Unless an error comes up, FoMpy is now installed on your virtual environment.


**Note: Most of the packages will be installed automatically during the FoMpy instalation. If you experience some issue, you can try to install the needed modules them yourself by typing in a terminal**::
	
	pip3 install setuptools
	pip3 install pytest
	pip3 install numpy
	pip3 install scipy
	pip3 install probscale
	pip3 install matplotlib
	sudo apt-get install python3-tk #optional


A simple example is included with the code so the user can test some basic commands and check the library
works as intended in their systems. After grabbing this repostiroty::

	git clone https://gitlab.citius.usc.es/gabriel.espineira/FoMPy/
	cd FoMPy-master

in the directory FoMpy-master, a file called ``example.py`` with command examples and a folder containing ensembles of simulated IV curves are included inside the path './data'. 

In order to test it comment and uncomment the lines that you want to run inside example.py and in a **python3 terminal** type::

	python3 example.py
	
Next, a full list of all the components of the FoMpy library will be presented.

=================================
Modules
=================================

.. automodule:: fompy.aux
    :members:
.. automodule:: fompy.conditioning
    :members:
.. automodule:: fompy.fds
    :members:
.. automodule:: fompy.fom
    :members:
.. automodule:: fompy.plots
    :members:
.. automodule:: fompy.wrappers
    :members:


