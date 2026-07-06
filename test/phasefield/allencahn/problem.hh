// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup AllenCahnTests
 * \brief Test problem for the Allen-Cahn model discretized with PQ2 (hybrid CVFE) elements.
 */
#ifndef DUMUX_TEST_PHASEFIELD_ALLENCAHN_PROBLEM_HH
#define DUMUX_TEST_PHASEFIELD_ALLENCAHN_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/problem.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/concepts/ipdata_.hh>

#include <dumux/discretization/cvfe/hybrid/elementvariables.hh>
#include <dumux/discretization/cvfe/elementsolution.hh>
#include <dumux/phasefield/allencahn/flux.hh>
#include <dumux/phasefield/freeenergy/doublewell.hh>

namespace Dumux {

/*!
 * \ingroup AllenCahnTests
 * \brief Test problem for the Allen-Cahn model on PQ2 (hybrid CVFE) elements.
 *
 * The problem sets homogeneous Neumann (no-flux) boundary conditions on the
 * entire domain boundary, and adds the derivative of a double-well free
 * energy (scaled by the mobility) as a source term.
 */
template<class TypeTag>
class AllenCahnTestProblem : public Dumux::Experimental::Problem<TypeTag>
{
    using ParentType = Dumux::Experimental::Problem<TypeTag>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<ModelTraits::numEq()>;
    using OldElementSolution = Dumux::CVFEElementSolution<FVElementGeometry, PrimaryVariables>;

public:
    AllenCahnTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , freeEnergy_(getParam<Scalar>("Problem.EnergyScale"))
    {
        mobility_ = getParam<Scalar>("Problem.Mobility");
        surfaceTension_ = getParam<Scalar>("Problem.SurfaceTension");
    }

    /*!
     * \brief Update the stored previous-time-level solution, used by
     *        `source()` to evaluate the explicit (concave) part of the
     *        convex-concave split double-well reaction term. Must be called
     *        once with the initial solution before the first time step, and
     *        again after each converged time step (mirroring the `oldSol`
     *        bookkeeping already done in `main.cc`).
     */
    void updateOldSolution(const SolutionVector& sol)
    { oldSol_ = sol; }

    /*!
     * \brief Specifies which kind of boundary condition should be used for
     *        which equation on a given boundary position. We use a
     *        homogeneous Neumann (no-flux) condition everywhere.
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllFluxBoundary();
        return values;
    }

    /*!
     * \brief Evaluate the source term at the given interpolation point.
     *
     * Adds (minus) the mobility-scaled derivative of the double-well free
     * energy, using a convex-concave (Eyre) splitting: the convex part is
     * evaluated at the current (unknown) phase field — implicit, contributing
     * to the Jacobian — and the concave part at the previous time step's
     * phase field — explicit, a constant contribution as far as the current
     * Newton solve is concerned. This makes the discrete free energy
     * unconditionally non-increasing, regardless of time step size. The
     * interpolation point may either coincide with a control-volume local dof
     * (in which case we can look up the nodal phase-field value directly) or
     * with a finite-element quadrature point within the element interior (in
     * which case we interpolate the phase field from the shape functions).
     */
    template<class ElementVariables, class IpData>
    NumEqVector source(const FVElementGeometry& fvGeometry,
                       const ElementVariables& elemVars,
                       const IpData& ipData) const
    {
        NumEqVector values(0.0);

        Scalar phi = 0.0;
        Scalar phiOld = 0.0;

        OldElementSolution oldElemSol;
        oldElemSol.update(fvGeometry.element(), oldSol_, fvGeometry);

        if constexpr (Dumux::Concept::LocalDofIpData<IpData>)
        {
            phi = elemVars[ipData].phaseField();
            phiOld = oldElemSol[ipData.localDofIndex()][Indices::phaseFieldIdx];
        }
        else
        {
            const auto& ipCache = cache(elemVars, ipData);
            phi = AllenCahnFluxFunctionContext(fvGeometry, elemVars, ipCache).phaseField();

            phiOld = 0.0;
            const auto& shapeValues = ipCache.shapeValues();
            for (const auto& localDof : localDofs(fvGeometry))
                phiOld += shapeValues[localDof.index()][0]*oldElemSol[localDof.index()][Indices::phaseFieldIdx];
        }

        values[Indices::phaseFieldEqIdx] = mobility_*(-freeEnergy_.convexDerivative(phi) + freeEnergy_.concaveDerivative(phiOld));
        return values;
    }

    Scalar mobility() const
    { return mobility_; }

    Scalar surfaceTension() const
    { return surfaceTension_; }

    //! The free-energy density f(phi), for energy-dissipation checks.
    Scalar freeEnergyDensity(Scalar phi) const
    { return freeEnergy_.value(phi); }

private:
    Scalar mobility_;
    Scalar surfaceTension_;
    Dumux::PhaseField::FreeEnergy::DoubleWell<Scalar> freeEnergy_;
    SolutionVector oldSol_;
};

} // end namespace Dumux

#endif
