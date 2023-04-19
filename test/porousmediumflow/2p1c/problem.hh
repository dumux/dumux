// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPOneCTests
 * \brief Non-isothermal steam injection test problem for the 2p1cni model.
 */

#ifndef DUMUX_STEAM_INJECTIONPROBLEM_HH
#define DUMUX_STEAM_INJECTIONPROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup TwoPOneCTests
 * \brief Non-isothermal 2D problem where steam is injected on the lower left side of the domain.
 *
 */
template <class TypeTag>
class InjectionProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryTypes = Dumux::BoundaryTypes<ModelTraits::numEq()>;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    // copy some indices for convenience
    enum {
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        conti0EqIdx = Indices::conti0EqIdx,
        energyEqIdx = Indices::energyEqIdx,

        // phase state
        liquidPhaseOnly = Indices::liquidPhaseOnly
    };

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    InjectionProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    { FluidSystem::init(); }


    //! \copydoc Dumux::FVProblem::source()
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    { return NumEqVector(0.0); }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;

        if(globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_ || globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
           bcTypes.setAllDirichlet();
        else
           bcTypes.setAllNeumann();

         return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
       return initialAtPos(globalPos);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The sub-control volume face
     *
     * For this method, the \a values parameter stores the flux
     * in normal direction of each phase. Negative values mean influx.
     * E.g. for the mass balance that would the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    NumEqVector neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        const auto& ipGlobal = scvf.ipGlobal();

        if (ipGlobal[0] < eps_)
        {
            if (ipGlobal[1] > 2.0 - eps_ && ipGlobal[1] < 3.0 + eps_)
            {
                const Scalar massRate = 1e-1;
                values[conti0EqIdx] = -massRate;
                values[energyEqIdx] = -massRate * 2690e3;
            }
        }
        return values;
    }

    /*!
     * \brief Evaluates the initial values for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        const Scalar densityW = 1000.0;
        values[pressureIdx] = 101300.0 + (this->gridGeometry().bBoxMax()[1] - globalPos[1])*densityW*9.81; // hydrostatic pressure
        values[switchIdx] = 283.13;

        values.setState(liquidPhaseOnly);

        return values;
    }

private:

    static constexpr Scalar eps_ = 1e-6;
};
} // end namespace Dumux

#endif
