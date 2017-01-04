// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/**
 * \file
 * \brief Definition of a problem, for the linear elasticity problem:
 * Problem definition for the deformation of an elastic solid.
 */
#ifndef DUMUX_ELASTICMATRIXPROBLEM_HH
#define DUMUX_ELASTICMATRIXPROBLEM_HH

#include <dumux/geomechanics/elastic/model.hh>
#include <dumux/implicit/fem/properties.hh>
#include <dumux/implicit/fem/problem.hh>

#include "elasticspatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class ElasticMatrixProblem;

namespace Properties
{
NEW_TYPE_TAG(ElasticMatrixProblem, INHERITS_FROM(FemModel, Elastic, ElSpatialParams));

// Set the grid type
SET_TYPE_PROP(ElasticMatrixProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(ElasticMatrixProblem, Problem, Dumux::ElasticMatrixProblem<TypeTag>);
}

/*!
 * \ingroup ElasticBoxProblems
 * \ingroup ImplicitTestProblems
 *
 * \brief Problem definition for the deformation of an elastic matrix.
 *
 * The problem defined here leads to the following linear distribution of the
 * solid displacement: u(x,y,z) = 1/E * (x,0,-nu*z) which for the given grid
 * The numerical results can be verified analytically.
 *
 * The 3D domain given in linearelastic.dgf spans from (0,0,0) to (10,1,10).
 *
 * Dirichlet boundary conditions (u=0.0) are applied for the displacement in y-direction
 * in the whole domain, for the displacement in x-direction on the left boundary (x < eps)
 * and for the displacement in z-direction for the lower left edge (x<eps && z<eps).
 * On the remaining boundaries Neumann boundary conditions are applied.
 * The normal stress applied on each boundary results from solving the momentum balance
 * analytically for the solid displacement function u(x,y,z) = 1/E * (x,0,-nu*z).
 * This leads to the normal stresses: \f$ \boldsymbol{\sigma_{xx}} = 2\,\frac{\mu}{E} + \frac{\lambda}{E} ( 1-\nu)\f$,
 *                                       \f$ \boldsymbol{\sigma_{yy}} = \frac{\lambda}{E} ( 1-\nu)\f$,
 *                                       \f$ \boldsymbol{\sigma_{zz}} = -2\,\frac{\mu \nu}{E} + \frac{\lambda}{E}\f$.
 * The shear stresses are set to zero.
 *
 * This problem uses the \ref ElasticModel model.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_elastic -parameterFile ./test_elastic.input</tt>
 */

template <class TypeTag>
class ElasticMatrixProblem : public ImplicitFemProblem<TypeTag>
{
    using ParentType = ImplicitFemProblem<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using IpData = typename GET_PROP_TYPE(TypeTag, FemIntegrationPointData);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    // copy some indices for convenience
    enum
    {
        // indices of the primary variables
        uxIdx = Indices::uxIdx,
        uyIdx = Indices::uyIdx,
        uzIdx = Indices::uzIdx,

        // indices of the equations
        momentumXEqIdx = Indices::momentumXEqIdx,
        momentumYEqIdx = Indices::momentumYEqIdx,
        momentumZEqIdx = Indices::momentumZEqIdx
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    ElasticMatrixProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView) {}

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const char *name() const
    { return "elasticmatrix";}
    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        values.setAllNeumann();

        if(globalPos[0] < eps_)
        {
            values.setDirichlet(uyIdx);
            values.setDirichlet(uxIdx);
            if(dim == 3 && globalPos[2] < eps_)
                values.setDirichlet(uzIdx);
        }

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Dirichlet
     *        boundary segment.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann
     *        boundary segment..
     */
    PrimaryVariables neumann(const Element& element,
                             const Intersection& intersection,
                             const ElementSolutionVector& elemSol,
                             const IpData& ipData) const
    {
        if (ipData.ipGlobal()[0] < this->bBoxMax()[0] - eps_)
            return PrimaryVariables(0.0);

        // interpolate solution
        PrimaryVariables ipSol(0.0);
        for (int i = 0; i < elemSol.size(); ++i)
            ipSol += elemSol[i]*ipData.shapeValues(i);

        // get Lame parameters
        const auto& lameParams = this->spatialParams().lameParams(element, ipSol);
        Scalar lambda = lameParams[0];
        Scalar mu = lameParams[1];
        Scalar E = this->spatialParams().E(element, ipSol);
        Scalar nu = this->spatialParams().nu(element, ipSol);

        // calculate values of sigma in normal direction
        Dune::FieldMatrix<Scalar, dim, dim> sigma(0);
        sigma[0][0] = 2.0*mu + lambda*(1 - nu);
        sigma[1][1] = lambda*(1 - nu);
        if (dim == 3)
            sigma[2][2] = -2.0*mu*nu + lambda*(1 - nu);

        sigma *= -1.0/E;

        // determine normal vector of current face
        auto normal = intersection.unitOuterNormal(intersection.geometryInInside().local(ipData.ipLocal()));

        // use stress in normal direction as boundary condition
        PrimaryVariables values(0.0);
        sigma.mv(normal, values);

        return PrimaryVariables({-10.0, 0});
    }
    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     */
    PrimaryVariables sourceAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return initial_(globalPos); }

    // \}

private:
    // the internal method for the initial condition
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    { return PrimaryVariables(0.0); }

    static constexpr Scalar eps_ = 3e-6;
};
} //end namespace

#endif
