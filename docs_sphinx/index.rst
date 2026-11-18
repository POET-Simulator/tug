.. Tug documentation master file, created by
   sphinx-quickstart on Mon Aug 14 11:30:23 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Tug's documentation!
===============================

Welcome to the documentation of the TUG project, a simulation program
for solving transport equations in one- and two-dimensional uniform
grids using cell centered finite differences.

---------
Diffusion
---------

TUG can solve diffusion problems with heterogeneous and anisotropic
diffusion coefficients. The partial differential equation expressing
diffusion reads:

.. math::
   \frac{\partial C}{\partial t} =  \nabla \cdot \left[ \mathbf{\alpha} \nabla C \right]

In 2D, the equation reads:
   
.. math::
   \frac{\partial C}{\partial t} =  \frac{\partial}{\partial x}\left[ \alpha_x \frac{\partial C}{\partial x}\right] + \frac{\partial}{\partial y}\left[ \alpha_y \frac{\partial C}{\partial y}\right]

---------
Advection
---------

TUG simulates advective solute transport in porous media driven by Darcy flow velocities.
The fluid flow velocities :math:`\mathbf{q} = (q_x, q_y)^T` are governed by Darcy's law:

.. math::
   \mathbf{q} = -\mathbf{K} \nabla h

where :math:`\mathbf{K}` denotes the hydraulic conductivity tensor and :math:`h` is the hydraulic head.
The conservative transport equation expressing solute advection in a porous medium with porosity :math:`\phi` reads:

.. math::
   \phi \frac{\partial C}{\partial t} + \nabla \cdot \left( \mathbf{q} C \right) = 0

In two dimensions, this takes the form:

.. math::
   \phi \frac{\partial C}{\partial t} + \frac{\partial (q_x C)}{\partial x} + \frac{\partial (q_y C)}{\partial y} = 0

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
    contributors
