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
#ifndef DUMUX_PNM2P_PROBLEM_HH
#define DUMUX_PNM2P_PROBLEM_HH

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porenetworkflow/2p/model.hh>

// spatial params
#include <dumux/material/spatialparams/porenetwork/porenetwork2p.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/porenetworkflow/common/utilities.hh>


namespace Dumux
{
template <class TypeTag>
class DrainageProblem;

namespace Properties
{
// Create new type tags
namespace TTag {
#if ISOTHERMAL
struct DrainageProblem { using InheritsFrom = std::tuple<PNMTwoP>; };
#else
struct DrainageProblem { using InheritsFrom = std::tuple<PNMTwoPNI>; };
#endif
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DrainageProblem> { using type = Dumux::DrainageProblem<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DrainageProblem>
  {
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = Dumux::FluidSystems::H2OAir<Scalar, Dumux::Components::SimpleH2O<Scalar>>;
  };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DrainageProblem> { using type = Dune::FoamGrid<1, 3>; };

} // end namespace Properties

template <class TypeTag>
class DrainageProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Labels = GetPropType<TypeTag, Properties::Labels>;
    enum {
        pwIdx = Indices::pressureIdx,
        snIdx = Indices::saturationIdx,
        nPhaseIdx = FluidSystem::phase1Idx,

#if !ISOTHERMAL
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,
#endif
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<GridView::dimension>::Entity;

public:
    template<class SpatialParams>
    DrainageProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<SpatialParams> spatialParams)
    : ParentType(gridGeometry, spatialParams)
    {
        vtpOutputFrequency_ = getParam<int>("Problem.VtpOutputFrequency");
        useFixedPressureAndSaturationBoundary_ = getParam<bool>("Problem.UseFixedPressureAndSaturationBoundary", false);
        pc_ = getParam<Scalar>("Problem.CapillaryPressure");
    }

    /*!
     * \name Simulation steering
     */
    // \{

    /*!
     * \name Problem parameters
     */
    // \{

    bool shouldWriteOutput(const int timeStepIndex, const GridVariables& gridVariables) const
    {
        if (vtpOutputFrequency_ < 0)
            return true;

        if (vtpOutputFrequency_ == 0)
            return (timeStepIndex == 0 || gridVariables.gridFluxVarsCache().invasionState().hasChanged());
        else
            return (timeStepIndex % vtpOutputFrequency_ == 0 || gridVariables.gridFluxVarsCache().invasionState().hasChanged());
    }

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
    BoundaryTypes boundaryTypes(const Element &element, const SubControlVolume &scv) const
    {
        BoundaryTypes bcTypes;

        // If a global phase pressure difference (pn,inlet - pw,outlet) with fixed saturations is specified, use a Dirichlet BC here
        if (useFixedPressureAndSaturationBoundary_ && isInletPore_(scv))
            bcTypes.setAllDirichlet();
        else if (isOutletPore_(scv))
           bcTypes.setAllDirichlet();

        else // neuman for the remaining boundaries
           bcTypes.setAllNeumann();
#if !ISOTHERMAL
        bcTypes.setDirichlet(temperatureIdx);
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
    PrimaryVariables dirichlet(const Element &element,
                               const SubControlVolume &scv) const
    {
        PrimaryVariables values(0.0);
        values[pwIdx] = 1e5;
        values[snIdx] = 0.0;

        // If a global phase pressure difference (pn,inlet - pw,outlet) is specified and the saturation shall also be fixed, apply:
        // pw,inlet = pw,outlet = 1e5; pn,outlet = pw,outlet + pc(S=0) = pw,outlet; pn,inlet = pw,inlet + pc_
        if (useFixedPressureAndSaturationBoundary_ && isInletPore_(scv))
            values[snIdx] = 1.0 - GetPropType<TypeTag, Properties::SpatialParams>::MaterialLaw::sw(this->spatialParams().materialLawParams(element, scv, 0/*dummy value*/), pc_);

#if !ISOTHERMAL
        if(isInletPore_(scv))
            values[temperatureIdx] = 273.15 + 15;
        else
            values[temperatureIdx] = 273.15 + 10;
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
    PrimaryVariables source(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume &scv) const
    {
        PrimaryVariables values(0.0);

        // If we do not want to use global phase pressure difference with fixed saturations and pressures,
        // we can instead only fix the non-wetting phase pressure and allow the wetting phase saturation to changle freely
        // by applying a Nitsche-type boundary condition which tries to minimize the difference between the present pn and the given value
        if (!useFixedPressureAndSaturationBoundary_ && isInletPore_(scv))
            values[snIdx] = (elemVolVars[scv].pressure(nPhaseIdx) - (1e5 + pc_)) * 1e8;

        return values;
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
        values[pwIdx] = 1e5;

        // get global index of pore
        const auto dofIdxGlobal = this->gridGeometry().vertexMapper().index(vertex);
        if(isInletPore_(dofIdxGlobal))
            values[snIdx] = 0.5;
        else
            values[snIdx] = 0.0;

#if !ISOTHERMAL
        values[temperatureIdx] = 273.15 + 10;
#endif
        return values;
    }

    /*!
     * \brief Evaluate the initial invasion state of a pore throat
     *
     * Returns true for a invaded throat and false elsewise.
     */
    bool initialInvasionState(const Element &element) const
    { return false; }

    // \}

private:

    bool isInletPore_(const SubControlVolume& scv) const
    {
        return isInletPore_(scv.dofIndex());
    }

    bool isInletPore_(const std::size_t dofIdxGlobal) const
    {
        return this->gridGeometry().poreLabel(dofIdxGlobal) == Labels::inlet;
    }

    bool isOutletPore_(const SubControlVolume& scv) const
    {
        return this->gridGeometry().poreLabel(scv.dofIndex()) == Labels::outlet;
    }

    int vtpOutputFrequency_;
    bool useFixedPressureAndSaturationBoundary_;
    Scalar pc_;
};
} //end namespace

#endif
