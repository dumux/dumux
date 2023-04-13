// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief A test problem for the  one-phase two-component pore network model.
 */
#ifndef DUMUX_PNM1P2C_PROBLEM_HH
#define DUMUX_PNM1P2C_PROBLEM_HH

// base problem
#include <dumux/porousmediumflow/problem.hh>
// Pore network model
#include <dumux/porenetwork/1pnc/model.hh>

#include <dumux/common/boundarytypes.hh>

namespace Dumux
{
template <class TypeTag>
class PNMOnePTwoCProblem;

template <class TypeTag>
class PNMOnePTwoCProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Labels = GetPropType<TypeTag, Properties::Labels>;
    enum {
        // primary variable indices
        conti0EqIdx = Indices::conti0EqIdx,
        pressureIdx = Indices::pressureIdx,
        transportEqIdx = 1,

#if !ISOTHERMAL
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,
#endif

    };

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<GridView::dimension>::Entity;

public:
    template<class SpatialParams>
    PNMOnePTwoCProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<SpatialParams> spatialParams)
    : ParentType(gridGeometry, spatialParams)
    {
        // get some parameters from the input file
        inletPressure_ = getParam<Scalar>("Problem.InletPressure");
        outletPressure_ = getParam<Scalar>("Problem.OutletPressure");
        inletMoleFraction_ = getParam<Scalar>("Problem.InletMoleFraction");
        outletMoleFraction_ = getParam<Scalar>("Problem.OutletMoleFraction");

        FluidSystem::init();
    }

     /*!
     * \name Boundary conditions
     */
    // \{
    //! Specifies which kind of boundary condition should be used for
    //! which equation for a finite volume on the boundary.
    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const
    {
        BoundaryTypes bcTypes;
        if (isInletPore_(scv) || isOutletPore_(scv))
            bcTypes.setAllDirichlet();
#if !ISOTHERMAL
        bcTypes.setDirichlet(temperatureIdx);
#endif
        return bcTypes;
    }

        /*!
     * \brief Evaluate the boundary conditions for a Dirichlet
     *        control volume.
     *
     * \param values The Dirichlet values for the primary variables
     * \param vertex The vertex (pore body) for which the condition is evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichlet(const Element& element,
                               const SubControlVolume& scv) const
    {
        PrimaryVariables values(0.0);
        if(isInletPore_(scv))
        {
            values[pressureIdx] = inletPressure_;
            values[transportEqIdx] = inletMoleFraction_;
        }
        else
        {
            values[pressureIdx] = outletPressure_;
            values[transportEqIdx] = outletMoleFraction_;
        }

#if !ISOTHERMAL
         if(isInletPore_(scv))
            values[temperatureIdx] = 273.15 +25;
         else
            values[temperatureIdx] = 273.15 +20;
#endif
         return values;
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initial(const Vertex& vertex) const
    {
        PrimaryVariables values(0.0);
            values[pressureIdx] = 1e5;
#if !ISOTHERMAL
            values[temperatureIdx] = 273.15 +20;
#endif
            return values;
    }

private:

    bool isInletPore_(const SubControlVolume& scv) const
    {
        return this->gridGeometry().poreLabel(scv.dofIndex()) == Labels::inlet;
    }

    bool isOutletPore_(const SubControlVolume& scv) const
    {
        return this->gridGeometry().poreLabel(scv.dofIndex()) == Labels::outlet;
    }

    Scalar inletPressure_;
    Scalar outletPressure_;
    Scalar inletMoleFraction_;
    Scalar outletMoleFraction_;
};
} //end namespace Dumux

#endif
