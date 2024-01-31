#######################
Lagrangian Trajectories
#######################
info:
Lagrangian trajectories are a simplified tool that uses the three components of the wind to interpolate designed for performing Lagrangian analysis in atmospheric science.
It uses wind data to trace the air mass motion. This tool requires Input data in a netCDF file in gridded format. The seeding points are the points in x,y, and z directions (longitude, latitude, pressure) of which the lagrangian trajectories are computed in forward or backward.

###########
References
###########
- `Rivi`ere et al. (2021) <https://doi.org/10.5194/wcd-2-1011-2021>`_:
- `Wimmer. M et al. (2022) <https://doi.org/10.5194/wcd-3-863-2022>`_:

###########
Authors
###########
- Vinita Deshmukh 
- Meryl Wimmer 
- Philippe Arbogast
- Léo Ducongé
- Gwendal Rivière
- Sébastien Fromang


###########
License
###########
This tool is open-source under the [your chosen license]. Feel free to use and modify it according to the license terms.
| Please cite **LagraTraj** in your publication: https://github.com/Vinita-D/LagranTraj/*.

============
Installation
============

Using pip
---------

Ideally install it in a virtual environment (development version, master).

.. code:: bash

    pip install git+https://github.com/Vinita-D/LagranTraj/*

Run the tests:

.. code-block:: bash

    python -m pytest

==========
Tutorial
==========

Example: Computing trajectories
---------------------------------------








