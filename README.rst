FoMpy: A figure of merit extraction tool for semiconductor device simulations
===============

.. What is this?
.. +++++++++++++

* A `FoMpy <https://gitlab.citius.usc.es/gabriel.espineira/FoMPy/>`__ tutorial introduction.
.. * By `Gabriel Espiñeira <https://github.com/>`__.
.. * February 20, 2019.

You will learn how to use the basic capabilities of the FoMpy library.

1. Introduction
---------------

FoMPy is an effective tool that extracts the main figures of merit (FoM) of a semiconductor's IV curve and provides useful statistical parameters for variability studies. It includes several methods to extract the threshold voltage. 

In the figure below the user can see the basic workflow behind the FoMpy library:


.. figure:: ./docs/figs/simplified_diagram.jpg
    :width: 600px
    :align: center
    :height: 400px
    :alt: alternate text
    :figclass: align-center

After loading the data into a FoMpy Dataset, using the various tools implemented in the library, the user is able to process and extract important parameters of the given curves.

2. Installation
---------------

The source code of FoMpy can be downloaded from https://gitlab.citius.usc.es/gabriel.espineira/FoMPy/ and is purely intended to run as a library for Linux systems (Debian, Ubuntu, FreeBSD, ...). FoMpy uses several external libraries that need to be installed in order to be able to use the full functionality of the tool. A list of these libraries can be seen below:

setuptools>=21.0.0 

numpy>=1.10.0

scipy>=0.8.0

matplotlib>=3.0.2

First you need to have installed pip on your system. Open up a terminal and type::

	sudo apt update

	sudo apt install python3-pip

Once the installation is complete, verify the installation by checking the pip version::

	pip3 --version

The use of virtual enviroments is highly encouraged. In order to use them run the following commands in a terminal::


	sudo apt install python3-venv
	python3 -m venv .venv
	source .venv/bin/activate

Note that as of this moment you're inside a virtual enviroment (Notice (.venv) $ in the terminal) with a limited version of python and therefore you will have to install all the packages you need (including the ones mentioned above because they are installed in the system, not the virtual enviroment).

.. and then::
	
.. 	python -m pip3 install -r requirements.txt

To check the version available of a package on your system type in a terminal::

	python3 -c "import numpy; print(numpy.version.version)"

and the version that is currently installed of numpy on your system will be printed.


.. Run in a terminal again::

.. 	pip install <library>
..	sudo apt install python3-tk


.. Via pip (recommended)
Via pip (recommended)

Run the following command in a terminal::

	pip3 install --extra-index-url https://test.pypi.org/simple/ fompy

and check the library is installed by importing it from a **python3 terminal**::

 	import fompy

Unless an error comes up, FoMpy is now installed on your virtual enviroment.

.. Via conda (not working)
 
.. Run the following command in a terminal::
 
.. 	conda search fompy

.. 	conda install fompy
 
.. Via source code (not working)
 
.. Go to https://github.com/ and download the project. Go to the parent folder and run::
 
..	make
	
.. make install


3. Quickstart 
-------------

In this section the user can learn the most basic yet powerful commands implemented in the FoMpy library. In order to do so either start by reading the basic commands or 
download and try the exampled provided in the repository explained at the end of this page.

Basic commands
+++++++++++++++++

A bunch of useful FoMpy commands are now provided. Supported tools include fompy.extract, fompy.plot or fompy.savetotxt. Here are some quick examples of the core capabilities of FoMpy:

In order to load a FoMpy Dataset run inside a **python3 terminal**::

	import fompy

FoMpy implements an importing tools that allows the user to extract the data from various sources
(from a file, an array stored in memory, a Sentaurus output file, etc). Inside the folder './data/' the user has to store all simulations in individual folders (i.e. './data/sim_1/current_file_1.txt', './data/sim_2/current_file_2.txt', etc)::

	path_data = './data'
	fds = fompy.dataset(path_data, parser=fompy.JCJB)

Note that the defined path has to point to the parent directory of the folders containing our single curve files.

After running this, a Fompy Dataset is created and the IV curves are stored inside it.
They can be accessed by calling the dataset attribute::

	print(fds.dataset)

Now that the Fompy Dataset has been implemented several other parameters can be defined like the
number of simulations (fds.n_sims) or a value for normalizing the curves (fds.norm)., the default extraction
method (fds.ext_method), the drain bias for the ensemble of curves (fds.drain_bias), the drain bias value
(fds.drain_bias_value) and the default interpolation method (fds.interpolation). All these parameters can be defined/updated
like the following example (Note that some of them will be defined automatically, like the number of simulations,
once the IV curves are loaded)::

	fds.drain_bias_value = 0.66

Also a predefined function can be called in order to print the current value of the attributes of the selected Fompy Dataset::

	fds.print_parameters()

The most important capability of Fompy is that it allows the user to extract the most common figures of merit (LATEX FOM)
of a semiconductor's IV curve using different methodologies. In order to extract these FoM the user has to call the 
function extract. The following example extracts the threshold voltage values (LATEX VTH) of the curves in the Fompy Dataset::

	vth_array = fompy.extract(fds, fom = 'vth')

and write the results to a file::

	fompy.savetotxt('./results_vth.txt', 'vth', vth_array)

Note that since no extraction method has been defined the library uses the second derivative method ('SD') as a default. 
This can be changed to oterh commonly used methods like the constant current method, the third derivative or the linear extrapolation (See further instructions
on how to choose this in the full documentation).

FoMpy also has built-in several plotting capabilities to be able to check the extraction results. A simple plot
of the threshold voltage with the 'SD' method and the second derivative of the curve goes as follows::

	fompy.plot(fds, fom = 'vth', save_plot='./vth_plots/sd/')

Note that the plots have been saved to the path './vth_plots/sd/', keeping the indexing of the curves as stored in the Fompy Dataset.



Repository Example
+++++++++++++++++++


A simple example is included with the code so the user can test some basic commands and check the library
works as intended in their systems. After grabbing this repostiroty::

	git clone https://gitlab.citius.usc.es/gabriel.espineira/FoMPy/
	cd FoMPy-master

in the directory FoMpy-master, a file called ``example.py`` with command examples and a folder containing ensembles of simulated IV curves are included inside the path './data'. 

In order to test it comment and uncomment the lines that you want to run inside example.py and in a **python3 terminal** type::

	python3 example.py

**Further documentation on the FoMpy library can be found inside ./docs/_build/latex/FoMpy.pdf**


