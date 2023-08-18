Developer Guide
===============

=========================
Class Diagram of user API
=========================

The following graphic shows the class diagram of the user API. The FTCS and 
BTCS functionalities are externally outsourced and not visible to the user.

.. image:: images/class_diagram.svg
    :width: 2000
    :alt: Class diagram for the user API

====================================================
Activity Diagram for run routine in simulation class
====================================================

The following activity diagram represents the actions when the run method is called within the simulation class. 
For better distinction, the activities of the calculation methods FTCS and BTCS are shown in two separate activity diagrams.

.. image:: images/activity_diagram_run.svg
    :width: 2000
    :alt: Activity diagram for the run method in the simulation class


**Activity Diagram for FTCS method**

.. image:: images/activity_diagram_FTCS.svg
    :width: 400
    :alt: Activity diagram for the FTCS method


**Activity Diagram for BTCS method**

.. image:: images/activity_diagram_BTCS.svg
    :width: 400
    :alt: Activity diagram for the BTCS method
