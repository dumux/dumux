// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
/*!
 * \file
 *
 * \brief A test problem for the one-phase pore network model.
 */
#ifndef DUMUX_PNM1P_PROBLEM_HH
#define DUMUX_PNM1P_PROBLEM_HH

#include <dune/foamgrid/foamgrid.hh>

// base problem
#include <dumux/porousmediumflow/problem.hh>
// Pore network model
#include <dumux/porenetworkflow/1p/model.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

namespace Dumux
{
template <class TypeTag>
class PNMOnePProblem;

namespace Properties
{

namespace TTag {
#if ISOTHERMAL
struct PNMOnePProblem { using InheritsFrom = std::tuple<PNMOneP>; };
#else
struct PNMOnePProblem { using InheritsFrom = std::tuple<PNMOnePNI>; };
#endif
}

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PNMOnePProblem>
{ using type = Dumux::PNMOnePProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PNMOnePProblem>
{
    using Scalar = GetPropType<TypeTag, Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Dumux::Components::SimpleH2O<Scalar> > ;
};

// the grid
template<class TypeTag>
struct Grid<TypeTag, TTag::PNMOnePProblem>
{ using type = Dune::FoamGrid<1, 3>; };

template<class TypeTag>
struct SinglePhaseTransmissibilityLaw<TypeTag, TTag::PNMOnePProblem>
{
    using Scalar = GetPropType<TypeTag, Scalar>;
    using type = TransmissibilityAzzamDullien<Scalar>;
};

// use the incompressible local residual (provides analytic jacobian)
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::PNMOnePProblem> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

}

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
    }

    /*!
     * \name Simulation steering
     */
    // \{

#if ISOTHERMAL
    /*!
     * \brief Return the temperature within the domain in [K].
     *
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10°C
    // \}
#endif
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
        else // neuman for the remaining boundaries
            bcTypes.setAllNeumann();
#if !ISOTHERMAL
        bcTypes.setDirichlet(Indices::temperatureIdx);
#endif
        return bcTypes;
    }

        /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param vertex The vertex (pore body) for which the condition is evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichlet(const Element& element,
                               const SubControlVolume& scv) const
    {
        PrimaryVariables values(0.0);
        if(isInletPore_(scv))
            values[Indices::pressureIdx] = 1.0e5;
        else
            values[Indices::pressureIdx] = 0.9e5;
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
};
} //end namespace

#endif
