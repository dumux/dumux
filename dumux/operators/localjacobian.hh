// $Id$

#ifndef DUNE_LOCALJACOBIAN_HH
#define DUNE_LOCALJACOBIAN_HH

#include<iostream>
#include<vector>
#include<set>
#include<map>
#include<stdio.h>
#include<stdlib.h>

#include<dune/common/timer.hh>
#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/fixedarray.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/grid.hh>
#include<dune/istl/operators.hh>
#include"dumux/operators/localstiffnessextended.hh"
#include<dune/disc/operators/boundaryconditions.hh>
#include"dumux/fvgeometry/fvelementgeometry.hh"

/**
* @file
* @brief  defines a class for piecewise linear finite element functions
* @author Peter Bastian
*/

/*! @defgroup DISC_Operators Operators
@ingroup DISC
@brief

@section D1 Introduction
<!--=================-->

To be written
*/

namespace Dune
{
/** @addtogroup DISC_Operators
*
* @{
*/
/**
* @brief base class for assembling local jacobian matrices
*
*/


/*! @brief Base class for local assemblers

This class serves as a base class for local assemblers. It provides
space and access to the local jacobian matrix. The actual assembling is done
in a derived class via the virtual assemble method.

The template parameters are:
- Imp  The implementation of the Interface with Barton-Nackman
- Grid    A grid type
- Scalar   The field type used in the elements of the jacobian matrix
- numEq    number of degrees of freedom per node (system size)
*/
template<class Imp, class Grid, class Scalar, int numEq>
class LocalJacobian : public LocalStiffness<Imp, Grid, Scalar, numEq>
{
    // grid types
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    enum {dim=Grid::dimension};

public:
    // types for matrics, vectors and boundary conditions
    typedef LocalStiffness<Imp, Grid, Scalar, numEq> ThisLocalStiffness;
    typedef FieldVector<Scalar,numEq> SolutionVector;                        // one entry in the global vectors
    typedef array<BoundaryConditions::Flags,numEq> BCBlockType; // componentwise boundary conditions
    typedef FVElementGeometry<Grid> ElementGeometry;
    enum {SIZE=8};

    void computeElementData(const Element& element) {
        return this->getImp().computeElementData(element);
    }

    virtual void updateVariableData(const Element& element, const SolutionVector* sol, int i, bool old = false)
    {
        return this->getImp().updateVariableData(element, sol, i, old);
    }

    void updateVariableData(const Element& element, const SolutionVector* sol, bool old = false)
    {
        return this->getImp().updateVariableData(element, sol, old);
    }

    template<class TypeTag>
    void assemble (const Element& element, int k = 1)
    {
        fvGeom.update(element);
        computeElementData(element);

        int size = element.template count<dim>();

        // set to Zero
        for (int i=0; i < size; i++) {
            this->bctype[i].assign(BoundaryConditions::neumann);
            this->b[i] = 0;
            this->def[i] = 0;
        }

        setLocalSolution(element);

        updateStaticData(element, u);
        bool old = true;
        updateVariableData(element, uold, old);
        updateVariableData(element, u);

        localDefect<TypeTag>(element, u);

        SolutionVector bTemp[size];
        for (int i=0; i<size; i++)
        {
            for (int equationnumber = 0; equationnumber < numEq; equationnumber++)
                if (this->bctype[i][equationnumber]==BoundaryConditions::neumann)
                    bTemp[i][equationnumber] = this->def[i][equationnumber];
                else
                    bTemp[i][equationnumber] = 0;
        }

        if (analytic) {
            analyticJacobian<TypeTag>(element, u);
        }
        else {
            SolutionVector defu[size];
            for (int i = 0; i < size; i++)
                defu[i] = def[i];
            SolutionVector uPlusEps[size];
            SolutionVector uMinusEps[size];

            for (int j = 0; j < size; j++)
                for (int comp = 0; comp < numEq; comp++)
                {
                    Scalar eps = std::max(fabs(1e-5*u[j][comp]), 1e-5);
                    for (int i = 0; i < size; i++) {
                        uPlusEps[i] = u[i];
                        uMinusEps[i] = u[i];
                    }
                    uPlusEps[j][comp] += eps;
                    uMinusEps[j][comp] -= eps;

                    updateVariableData(element, uPlusEps, j);

                    // calculate the defect without taking into account BCs
                    // ASSUMES that BCs do not depend on the solution
                    bool withoutBC = false;
                    localDefect<TypeTag>(element, uPlusEps, withoutBC);
                    SolutionVector defuPlusEps[size];
                    for (int i = 0; i < size; i++)
                        defuPlusEps[i] = def[i];

                    updateVariableData(element, uMinusEps, j);
                    localDefect<TypeTag>(element, uMinusEps, withoutBC);

                    updateVariableData(element, u, j);

                    Scalar oneByEps = 0.5/eps;
                    for (int i = 0; i < size; i++)
                        for (int compi = 0; compi < numEq; compi++)
                            this->A[i][j][compi][comp] = oneByEps*(defuPlusEps[i][compi] - def[i][compi]);
                }
        }

        //        for (int i = 0; i < size; i++)
        //            for (int compi = 0; compi < numEq; compi++) {
        //                for (int j = 0; j < size; j++) {
        //                    for (int compj = 0; compj < numEq; compj++)
        //                        std::cout << std::setw(9) << this->A[i][j][compi][compj] << ", ";
        //                    std::cout << "\t";
        //                }
        //                std::cout << std::endl;
        //            }


        for (int i=0; i<size; i++)
        {
            for (int equationnumber = 0; equationnumber < numEq; equationnumber++)
                if (this->bctype[i][equationnumber]==BoundaryConditions::neumann) {
                    this->b[i][equationnumber] = bTemp[i][equationnumber];
                }
        }

        return;
    }

    template<class TypeTag>
    void localDefect (const Element& element, const SolutionVector* sol, bool withBC = true)
    {
        this->getImp().template localDefect<TypeTag>(element, sol, withBC);
    }

    void setLocalSolution (const Element& element)
    {
        this->getImp().setLocalSolution(element);
    }

    template<class TypeTag>
    void assembleBC (const Element& element)
    {
        this->getImp().template assembleBC<TypeTag>(element);
    }

    template<class TypeTag>
    void analyticJacobian (const Element& element, const SolutionVector* sol)
    {
        this->getImp().template analyticJacobian<TypeTag>(element, sol);
    }

    virtual void updateStaticData (const Element& element, SolutionVector* sol)
    {
        return;
    }

    virtual void clearVisited ()
    {
        return;
    }

    //! access defect for each dof
    /*! Access defect for each degree of freedom. Elements are
    undefined without prior call to the assemble method.
    */
    const SolutionVector& defect (int i) const
    {
        return def[i];
    }

    void setDt(double d)
    {
        this->getImp().setDt(d);
    }

    virtual double getDt()
    {
        return this->getImp().getDt();
    }

    ElementGeometry fvGeom;
    SolutionVector def[SIZE];
    SolutionVector u[SIZE];
    SolutionVector uold[SIZE];
    bool analytic;
};


/** @} */

}
#endif
