// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief A test problem for the one-phase pore network model.
 */
#ifndef DUMUX_PNM1P_PROBLEM_HH
#define DUMUX_PNM1P_PROBLEM_HH

// base problem
#include <dumux/porousmediumflow/problem.hh>
// Pore network model
#include <dumux/porenetwork/1p/model.hh>

#include <dumux/common/boundarytypes.hh>

namespace Dumux {
template <class TypeTag>
class PNMOnePProblem;

template <class TypeTag>
class PNMOnePProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Labels = GetPropType<TypeTag, Properties::Labels>;

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<GridView::dimension>::Entity;

public:
    template<class SpatialParams>
    PNMOnePProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<SpatialParams> spatialParams)
    : ParentType(gridGeometry, spatialParams)
    {
        testHydrostaticPressure_ = getParamFromGroup<bool>(this->paramGroup(), "Problem.EnableGravity");
        inletPressure_ = getParam<Scalar>("Problem.InletPressure");
        outletPressure_ = getParam<Scalar>("Problem.OutletPressure");
    }

    /*!
     * \name Simulation steering
     */
    // \{

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
        bcTypes.setDirichlet(Indices::temperatureIdx);
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
            values[Indices::pressureIdx] = inletPressure_;
        else
            values[Indices::pressureIdx] = outletPressure_;
#if !ISOTHERMAL
         if(isInletPore_(scv))
            values[Indices::temperatureIdx] = 273.15 +25;
         else
            values[Indices::temperatureIdx] = 273.15 +20;
#endif
         return values;
    }


    // \}

    /*!
     * \name Volume terms
     */
    // \{
    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub control volume
     *
     * For this method, the return parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
    template<class ElementVolumeVariables>
    PrimaryVariables source(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume &scv) const
    {
        return PrimaryVariables(0);
    }

    // \}

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initial(const Vertex& vertex) const
    {
        PrimaryVariables values(0.0);
            values[Indices::pressureIdx] = 1e5;
#if !ISOTHERMAL
            values[Indices::temperatureIdx] = 273.15 +20;
#endif
            return values;
    }

    // \}

private:

    bool isInletPore_(const SubControlVolume& scv) const
    {
        return this->gridGeometry().poreLabel(scv.dofIndex()) == Labels::inlet;
    }

    bool isOutletPore_(const SubControlVolume& scv) const
    {
        if (testHydrostaticPressure_)
            return false;

        return this->gridGeometry().poreLabel(scv.dofIndex()) == Labels::outlet;
    }

    bool testHydrostaticPressure_;
    Scalar inletPressure_;
    Scalar outletPressure_;
};
} //end namespace Dumux

#endif
