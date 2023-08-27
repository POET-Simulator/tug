.. Tug documentation master file, created by
   sphinx-quickstart on Mon Aug 14 11:30:23 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Tug's documentation!
===============================

Welcome to the documentation of the TUG project, a simulation program
for solving transport equations in one- and two-dimensional uniform
grids using cell centered finite differences.

Diffusion
-----------

TUG can solve diffusion problems with heterogeneous and anisotropic
diffusion coefficients. The partial differential equation expressing
diffusion reads:

.. math::
   \frac{\partial C}{\partial t} =  \nabla \cdot \left[ \mathbf{\alpha} \nabla C \right]

In 2D, the equation reads:
   
.. math::
   \frac{\partial C}{\partial t} =  \frac{\partial}{\partial x}\left[ \alpha_x \frac{\partial C}{\partial x}\right] + \frac{\partial}{\partial y}\left[ \alpha_y \frac{\partial C}{\partial y}\right]

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

Table of Contents
^^^^^^^^^^^^^^^^^
.. toctree::    

    :maxdepth: 2
    self
    installation
    theory
    user
    developer
    examples
    visualization
