/**
 * @file Boundary.hpp
 * @brief 
 * 
 * 
 */
#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include <Eigen/Core>
#include <cstddef>
#include "Grid.hpp"

using namespace std;
using namespace Eigen;

enum BC_TYPE {
    BC_TYPE_CLOSED,
    BC_TYPE_CONSTANT
};

enum BC_SIDE {
    BC_SIDE_LEFT,
    BC_SIDE_RIGHT,
    BC_SIDE_TOP,
    BC_SIDE_BOTTOM
};

/**
 * This class defines the boundary conditions of individual boundary elements.
 * These can be flexibly used and combined later in other classes.
 * The class serves as an auxiliary class for structuring the Boundary class.
 */
class BoundaryElement {
    public:
        // bc type closed
      /**
       * @brief Construct a new Boundary Element object for the closed case.
       *        The boundary type is here automatically set to the type
       *        BC_TYPE_CLOSED, where the value takes NaN.
       */
      BoundaryElement();

      /**
       * @brief Construct a new Boundary Element object for the constant case.
       *        The boundary type is automatically set to the type
       * BC_TYPE_CONSTANT.
       *
       * @param value Value of the constant concentration to be assumed at the
       *              corresponding boundary element. 
       */
      BoundaryElement(double value);

      /**
       * @brief Allows changing the boundary type of a corresponding
       *        BoundaryElement object.
       *
       * @param type Type of boundary condition. Either BC_TYPE_CONSTANT or
                     BC_TYPE_CLOSED. 
       */
      void setType(BC_TYPE type);
      
      /**
       * @brief Sets the value of a boundary condition for the constant case.
       * 
       * @param value Concentration to be considered constant for the 
       *              corresponding boundary element.
       */
      void setValue(double value);

      /**
       * @brief Return the type of the boundary condition, i.e. whether the
       *        boundary is considered closed or constant.
       *
       * @return BC_TYPE Type of boundary condition, either BC_TYPE_CLOSED or
                 BC_TYPE_CONSTANT.
       */
      BC_TYPE getType();

      /**
       * @brief Return the concentration value for the constant boundary condition.
       * 
       * @return double Value of the concentration.
       */
      double getValue();

    private:
        BC_TYPE type;
        double value;
};


/**
 * This class implements the functionality and management of the boundary
 * conditions in the grid to be simulated.
 * This class implements the functionality and management of the boundary
 * conditions in the grid to be simulated.
 */
class Boundary {
    public:
      /**
       * @brief Creates a boundary object based on the passed grid object and
       *        initializes the boundaries as closed.
       *
       * @param grid Grid object on the basis of which the simulation is to take place.
       */
      Boundary(Grid grid);

      /**
       * @brief Sets all elements of the specified boundary side to the boundary
       *        condition closed.
       *
       * @param side Side to be set to closed, e.g. BC_SIDE_LEFT.
       */
      void setBoundarySideClosed(BC_SIDE side);

      /**
       * @brief Sets all elements of the specified boundary side to the boundary
       *        condition constant. Thereby the concentration values of the
       *        boundaries are set to the passed value.
       *
       * @param side Side to be set to constant, e.g. BC_SIDE_LEFT.
       * @param value Concentration to be set for all elements of the specified page.
       */
      void setBoundarySideConstant(BC_SIDE side, double value);

      /**
       * @brief Specifically sets the boundary element of the specified side
       *        defined by the index to the boundary condition closed.
       *
       * @param side Side in which an element is to be defined as closed.
       * @param index Index of the boundary element on the corresponding
       *              boundary side. Must index an element of the corresponding side.
       */
      void setBoundaryElementClosed(BC_SIDE side, int index);

      /**
       * @brief Specifically sets the boundary element of the specified side
       *        defined by the index to the boundary condition constant with the
                given concentration value.
       * 
       * @param side Side in which an element is to be defined as constant.
       * @param index Index of the boundary element on the corresponding
       *              boundary side. Must index an element of the corresponding side.
       * @param value Concentration value to which the boundary element should be set.
       */
      void setBoundaryElementConstant(BC_SIDE side, int index, double value);

      /**
       * @brief Returns the boundary condition of a specified side as a vector
       *        of BoundarsElement objects.
       *
       * @param side Boundary side from which the boundaryconditions are to be returned.
       * @return vector<BoundaryElement> Contains the boundary conditions as BoundaryElement objects.
       */
      vector<BoundaryElement> getBoundarySide(BC_SIDE side);

      /**
       * @brief Returns the boundary condition of a specified element on a given side.
       * 
       * @param side Boundary side in which the boundary condition is located.
       * @param index Index of the boundary element on the corresponding
       *              boundary side. Must index an element of the corresponding side.
       * @return BoundaryElement Boundary condition as a BoundaryElement object.
       */
      BoundaryElement getBoundaryElement(BC_SIDE side, int index);

      /**
       * @brief Returns the type of a boundary condition, i.e. either BC_TYPE_CLOSED or
                BC_TYPE_CONSTANT.
       * 
       * @param side Boundary side in which the boundary condition type is located.
       * @param index Index of the boundary element on the corresponding
       *              boundary side. Must index an element of the corresponding side.
       * @return BC_TYPE Boundary Type of the corresponding boundary condition.
       */  
      BC_TYPE getBoundaryElementType(BC_SIDE side, int index);

      /**
       * @brief Returns the concentration value of a corresponding
       *        BoundaryElement object if it is a constant boundary condition.
       *
       * @param side Boundary side in which the boundary condition value is
       *             located.
       * @param index Index of the boundary element on the corresponding
       *              boundary side. Must index an element of the corresponding
       *              side.
       * @return double Concentration of the corresponding BoundaryElement object.
       */
      double getBoundaryElementValue(BC_SIDE side, int index);

    private:
        Grid grid;
        
        vector<vector<BoundaryElement>> boundaries;
};

#endif

