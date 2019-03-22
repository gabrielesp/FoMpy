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

In the figure below the used can see the basic workflow behind the FoMpy library:


.. figure:: ./figs/simplified_diagram.jpg
    :width: 600px
    :align: center
    :height: 300px
    :alt: alternate text
    :figclass: align-center

    Diagram showing the basic methodology of the FoMpy library. After loading the data into a FoMpy Dataset, using the various tools implemented in the library, the user is able to process and extract important parameters from a given dataset.

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

