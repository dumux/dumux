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

// Quadrature order
SET_INT_PROP(ElasticMatrixProblem, FemBasisOrder, 1);

SET_TYPE_PROP(ElasticMatrixProblem, LinearSolver, UMFPackBackend<TypeTag>);
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
    using SecondaryVariables = typename GET_PROP_TYPE(TypeTag, SecondaryVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using StressTensor = Dune::FieldMatrix<Scalar, dim, dimWorld>;

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
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const char *name() const
    { return "elasticmatrix";}

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
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
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param ipData Data on the shape values and gradients at the integration point
     * \param secVars The primary/secondary variables evaluated at the integration point
     *
     */
    PrimaryVariables source(const Element& element,
                            const IpData& ipData,
                            const SecondaryVariables& secVars) const
    {
        using std::sin;
        using std::cos;

        static const Scalar pi = 3.14159265358979323846;
        const auto ipGlobal = ipData.ipGlobal();
        const auto x = ipGlobal[0];
        const auto y = ipGlobal[1];

        // the lame parameters
        const auto lambda = secVars.lambda();
        const auto mu = secVars.mu();

        // precalculated products
        const Scalar pi_2 = 2.0*pi;
        const Scalar pi_2_square = pi_2*pi_2;
        const Scalar cos_2pix = cos(pi_2*x);
        const Scalar sin_2pix = sin(pi_2*x);
        const Scalar cos_2piy = cos(pi_2*y);
        const Scalar sin_2piy = sin(pi_2*y);

        const Scalar dE11_dx = -2.0*sin_2piy;
        const Scalar dE22_dx = pi_2_square*cos_2pix*cos_2piy;
        const Scalar dE11_dy = pi_2*(1.0-2.0*x)*cos_2piy;
        const Scalar dE22_dy = -1.0*pi_2_square*sin_2pix*sin_2piy;
        const Scalar dE12_dy = 0.5*pi_2_square*(cos_2pix*cos_2piy - (x-x*x)*sin_2piy);
        const Scalar dE21_dx = 0.5*((1.0-2*x)*pi_2*cos_2piy - pi_2_square*sin_2pix*sin_2piy);

        // compute exact divergence of sigma
        PrimaryVariables divSigma(0.0);
        divSigma[Indices::momentum(uxIdx)] = lambda*(dE11_dx + dE22_dx) + 2*mu*(dE11_dx + dE12_dy);
        divSigma[Indices::momentum(uyIdx)] = lambda*(dE11_dy + dE22_dy) + 2*mu*(dE21_dx + dE22_dy);

        return divSigma;
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Evaluate the exact displacement to this problem at a given position.
     */
    PrimaryVariables exactSolution(const GlobalPosition& globalPos) const
    {
        using std::sin;

        static const Scalar pi = 3.14159265358979323846;
        const auto x = globalPos[0];
        const auto y = globalPos[1];

        PrimaryVariables exact(0.0);
        exact[uxIdx] = (x-x*x)*sin(2*pi*y);
        exact[uyIdx] = sin(2*pi*x)*sin(2*pi*y);
        return exact;
    }

    /*!
     * \brief Evaluate the exact displacement gradient to this problem at a given position.
     */
    StressTensor exactGradient(const GlobalPosition& globalPos) const
    {
        using std::sin;
        using std::cos;

        static const Scalar pi = 3.14159265358979323846;
        const auto x = globalPos[0];
        const auto y = globalPos[1];

        StressTensor exactGrad(0.0);
        exactGrad[uxIdx][0] = (1-2*x)*sin(2*pi*y);
        exactGrad[uxIdx][1] = (x - x*x)*2*pi*cos(2*pi*y);
        exactGrad[uyIdx][0] = 2*pi*cos(2*pi*x)*sin(2*pi*y);
        exactGrad[uyIdx][1] = 2*pi*sin(2*pi*x)*cos(2*pi*y);
        return exactGrad;
    }

private:
    static constexpr Scalar eps_ = 3e-6;
};
} //end namespace

#endif
