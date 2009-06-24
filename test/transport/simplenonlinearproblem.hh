// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   Copyright (C) 2009 by Yufei Cao                                         *
 *   Institute of Applied Analysis and Numerical Simulation                  *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUNE_SIMPLENONLINEARPROBLEM_HH
#define DUNE_SIMPLENONLINEARPROBLEM_HH

#include <dumux/transport/transportproblem.hh>

#include <dumux/material/matrixproperties.hh>
#include <dumux/material/property_baseclasses.hh>

namespace Dune
{

/** \todo Please doc me! */

template<class Grid, class Scalar>
class HomogeneousNonlinearSoil: public HomogeneousSoil<Grid,Scalar>
{
public:
    enum {dim=Grid::dimension, numEq=1};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef FieldVector<Scalar,dim> LocalPosition;
    typedef FieldVector<Scalar,dim> GlobalPosition;

    virtual double porosity(const GlobalPosition& globalPos, const Element& e, const LocalPosition& localPos) const
    {
        return 1.0;
    }

    virtual std::vector<double> paramRelPerm(const GlobalPosition& globalPos, const Element& e, const LocalPosition& localPos, const double T) const
    {
        std::vector<double> param(2);

        // Brooks-Corey law parameters
        param[0] = 2.; // lambda
        param[1] = 0.; // entry-pressure


        return param;
    }

    virtual typename Matrix2p<Grid,Scalar>::modelFlag relPermFlag(const GlobalPosition& globalPos, const Element& e, const LocalPosition& localPos) const
    {
        // Brooks-Corey law
        return Matrix2p<Grid,Scalar>::brooks_corey;
    }

    HomogeneousNonlinearSoil()
    : HomogeneousSoil<Grid,Scalar>()
    {}
};

//! \ingroup transportProblems
//! @brief example class for a transport problem
template<class GridView, class Scalar, class VariableClass>
class SimpleNonlinearProblem : public TransportProblem<GridView, Scalar, VariableClass> {
    enum {dim=GridView::dimension, numEq=1};
    typedef typename GridView::Grid Grid;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef FieldVector<Scalar,dim> LocalPosition;
    typedef FieldVector<Scalar,dim> GlobalPosition;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

private:
    Scalar left_;
    Scalar right_;

public:
    BoundaryConditions::Flags bctypeSat (const GlobalPosition& globalPos, const Element& element,
                                         const LocalPosition& localPos) const
    {
        if (globalPos[0] > right_-1E-8 || globalPos[0] < left_+1e-8)
            return Dune::BoundaryConditions::dirichlet;
        else
            return Dune::BoundaryConditions::neumann;
    }

    Scalar dirichletSat (const GlobalPosition& globalPos, const Element& element,
                         const LocalPosition& localPos) const
    {
        if (globalPos[0] < left_+1e-8)
            return 1;
        else
            return 0;
    }

    virtual Scalar neumannSat (const GlobalPosition& globalPos, const Element& element,
                               const LocalPosition& localPos, Scalar helpFactor) const
    {
        return 0;
    }

    Scalar initSat (const GlobalPosition& globalPos, const Element& element,
                    const LocalPosition& localPos) const
    {
        return 0;
    }

    virtual Scalar neumann(const GlobalPosition& globalPos, const Element& element,
            const IntersectionIterator& intersectionIt,
                               const LocalPosition& localPos) const
    {
        return neumannSat(globalPos, element, localPos, 0);
    }

    BoundaryConditions::Flags bctype (const GlobalPosition& globalPos, const Element& element,
                const IntersectionIterator& intersectionIt,
                                   const LocalPosition& localPos) const
    {
        return bctypeSat(globalPos, element, localPos);
    }

    SimpleNonlinearProblem(VariableClass& variables, Fluid& wettingPhase, Fluid& nonwettingPhase,
    		Matrix2p<Grid,Scalar>& soil, TwoPhaseRelations<Grid,Scalar>& materialLaw,
    		GlobalPosition& left = 0, GlobalPosition& right = 1)
    : TransportProblem<GridView, Scalar, VariableClass>(variables, wettingPhase, nonwettingPhase, soil, materialLaw),
	left_(left[0]), right_(right[0])
    {}

    SimpleNonlinearProblem(VariableClass& variables, TwoPhaseRelations<Grid,Scalar>& materialLaw,
    		GlobalPosition& left = 0, GlobalPosition& right = 1)
    : TransportProblem<GridView, Scalar, VariableClass>(variables, materialLaw),
	left_(left[0]), right_(right[0])
    {}
};

}
#endif
