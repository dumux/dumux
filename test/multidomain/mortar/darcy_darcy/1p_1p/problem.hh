// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \ingroup MultiDomainTests
 * \brief The problem for the single-phase darcy-darcy mortar-coupling test
 */
#ifndef DUMUX_MORTAR_DARCY_ONEP_TEST_PROBLEM_HH
#define DUMUX_MORTAR_DARCY_ONEP_TEST_PROBLEM_HH

#include <bitset>
#include <memory>

#include <dune/common/exceptions.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/discretization/method.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/multidomain/mortar/couplingmanager.hh>
#include <dumux/multidomain/mortar/porousmediumflow/traceassembly.hh>

namespace Dumux {

#ifndef PERMEABILITY_JUMP
#define PERMEABILITY_JUMP 0
#endif

/*!
 * \brief The exact solution of the coupled problem.
 *
 * Without a jump this is the harmonic x*y. With one it is piecewise linear across an
 * interface at x = 1/2, chosen so that both the pressure and the normal flux are continuous
 * there for permeabilities 1 and 2: p = 4x/3 on the left and 2x/3 + 1/3 on the right, giving
 * an interface pressure of 2/3 and a normal flux of -4/3 on either side.
 */
template<class GlobalPosition>
double mortarExactSolution(const GlobalPosition& p)
{
#if PERMEABILITY_JUMP
    return p[0] < 0.5 ? 4.0*p[0]/3.0 : 2.0*p[0]/3.0 + 1.0/3.0;
#else
    return p[0]*p[1];
#endif
}

/*!
 * \brief The exact normal flux across an axis-aligned interface whose normal points along the
 *        positive axis given by \a normalAxis: \f$-K\,\nabla p\cdot n\f$, evaluated with the
 *        permeability of the side the normal points away from.
 */
template<class GlobalPosition>
double mortarExactFlux(const GlobalPosition& p, int normalAxis)
{
#if PERMEABILITY_JUMP
    // the solution varies along x only, so nothing crosses an interface normal to y
    return normalAxis == 0 ? -4.0/3.0 : 0.0;
#else
    // the gradient of x*y is (y, x)
    return -p[1 - normalAxis];
#endif
}


/*!
 * \ingroup OnePTests
 * \ingroup MultiDomainTests
 * \brief The problem for the single-phase darcy-darcy mortar-coupling test
 */
template<class TypeTag>
class DarcyProblem
: public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;

    static constexpr bool isCVFE = DiscretizationMethods::isCVFE<typename GridGeometry::DiscretizationMethod>;

public:
    using MortarCouplingManager = Mortar::SubDomainCouplingManager<
        GridGeometry,
        GetPropType<TypeTag, Properties::MortarGrid>,
        GetPropType<TypeTag, Properties::MortarSolutionVector>
    >;

    DarcyProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                 std::shared_ptr<const MortarCouplingManager> mortarCouplingManager,
                 const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    , mortar_(std::move(mortarCouplingManager))
    {}

    template<typename GridVariables>
    auto assembleTraceVariables(std::size_t mortarId,
                                const GridVariables& gridVariables,
                                const SolutionVector& x) const
    {
        // in natural mode the mortar datum was imposed as a flux, so the conjugate value
        // trace is read back; in essential mode it is the other way around. The face
        // pressure is reconstructed from the cell value of a cell-centred scheme.
        if (mortar_->couplingMode() == Mortar::CouplingMode::natural)
        {
            if constexpr (isCVFE)
                DUNE_THROW(Dune::NotImplemented, "The value trace of a control-volume finite element subdomain in natural mode");
            else
                return Mortar::assembleTrace(*mortar_, mortarId, *this, gridVariables, x,
                                             Mortar::reconstructedFacePressure(*mortar_));
        }
        return Mortar::assembleTrace(*mortar_, mortarId, *this, gridVariables, x,
                                     Mortar::darcyMassFlux<FluxVariables>());
    }

    template<typename SubControlEntity>
    BoundaryTypes boundaryTypes(const Element& element, const SubControlEntity& sce) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        if (mortar_->isCoupled(element, sce) && mortar_->couplingMode() == Mortar::CouplingMode::natural)
            values.setAllNeumann();
        return values;
    }

    template<typename ElementVolumeVariables, typename ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        if (mortar_->isCoupled(element, scvf))
            return NumEqVector(mortar_->traceAt(element, scvf));
        return NumEqVector(0.0);
    }

    static constexpr bool enableInternalDirichletConstraints()
    { return !isCVFE; }

    std::bitset<PrimaryVariables::dimension>
    hasInternalDirichletConstraint(const Element& element, const SubControlVolume& scv) const
    {
        std::bitset<PrimaryVariables::dimension> constraint;
        if (mortar_->isPinnedDof(element, scv))
            constraint.set(0);
        return constraint;
    }

    PrimaryVariables internalDirichlet(const Element& element, const SubControlVolume& scv) const
    { return PrimaryVariables(0.0); }

    template<typename SubControlEntity>
    PrimaryVariables dirichlet(const Element& element, const SubControlEntity& sce) const
    {
        if (mortar_->isCoupled(element, sce))
            return mortar_->traceAt(element, sce);
        if (mortar_->isHomogeneous())
            return PrimaryVariables(0);
        if constexpr (isCVFE)
            return PrimaryVariables(mortarExactSolution(sce.dofPosition()));
        else
            return PrimaryVariables(mortarExactSolution(sce.ipGlobal()));
    }

private:
    std::shared_ptr<const MortarCouplingManager> mortar_;
};

} // end namespace Dumux

#endif
