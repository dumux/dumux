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
 * \brief A test problem for the two-phase pore network model.
 */
#ifndef DUMUX_PNM_2P_NC_PROBLEM_HH
#define DUMUX_PNM_2P_NC_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/parameters.hh>
#include <dumux/porenetwork/2pnc/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/components/air.hh>

namespace Dumux {

template <class TypeTag>
class DrainageProblem;

template <class TypeTag>
class DrainageProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Labels = GetPropType<TypeTag, Properties::Labels>;
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
        source_ = getParam<Scalar>("Problem.Source");
        inletPressure_ = getParam<Scalar>("Problem.InletPressure", 1e5);
        outletPressure_ = getParam<Scalar>("Problem.OutletPressure", 1e5);
#if !ISOTHERMAL
        inletTemperature_ = getParam<Scalar>("Problem.InletTemperature", 283.15);
        outletTemperature_ = getParam<Scalar>("Problem.OutletTemperature", 283.15);
#endif
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

     /*!
     * \name Boundary conditions
     */
    // \{
    //! Specifies which kind of boundary condition should be used for
    //! which equation for a finite volume on the boundary.
    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const
    {
        BoundaryTypes bcTypes;

#if ISOTHERMAL
        if (isInletPore_(scv))
#if EVAPORATION
            bcTypes.setAllDirichlet();
#else
            bcTypes.setAllNeumann();
#endif
#else
        if (isInletPore_(scv))
            bcTypes.setDirichlet(Indices::temperatureIdx);
#endif
        else if (isOutletPore_(scv))
            bcTypes.setAllDirichlet();

        return bcTypes;
    }


    //! Evaluate the boundary conditions for a Dirichlet control volume.
    PrimaryVariables dirichlet(const Element& element,
                               const SubControlVolume& scv) const
    {
        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = 1e5;
        values[Indices::switchIdx] = 0.0;

        // If a global phase pressure difference (pn,inlet - pw,outlet) is specified and the saturation shall also be fixed, apply:
        // pw,inlet = pw,outlet = 1e5; pn,outlet = pw,outlet + pc(S=0) = pw,outlet; pn,inlet = pw,inlet + pc_
        if (isInletPore_(scv))
        {
            values.setState(Indices::secondPhaseOnly);
            values[Indices::pressureIdx] = inletPressure_;
            values[Indices::switchIdx] = 0.0;
#if !ISOTHERMAL
            values[Indices::temperatureIdx] = inletTemperature_;
#endif
        }
        else if (isOutletPore_(scv))
        {
            values.setState(Indices::firstPhaseOnly);
            values[Indices::pressureIdx] = outletPressure_;
            values[Indices::switchIdx] = 0.0;
#if !ISOTHERMAL
            values[Indices::temperatureIdx] = outletTemperature_;
#endif
        }

        return values;
    }


    // \}

    /*!
     * \name Volume terms
     */
    // \{

    //! Evaluate the source term for all phases within a given sub-control-volume.
    PrimaryVariables source(const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume& scv) const
    {
        PrimaryVariables values(0.0);
        if (isInletPore_(scv))
            values[Indices::conti0EqIdx+1] = source_/sumInletPoresVolume_;
        return values;
    }
    // \}

    //! Evaluate the initial value for a control volume.
    PrimaryVariables initial(const Vertex& vertex) const
    {
        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = outletPressure_;

        // get global index of pore
        const auto dofIdxGlobal = this->gridGeometry().vertexMapper().index(vertex);
        if (isInletPore_(dofIdxGlobal))
        {
#if EVAPORATION
            values.setState(Indices::secondPhaseOnly);
#else
            values.setState(Indices::firstPhaseOnly);
#endif
            values[Indices::switchIdx] = 0.0;
        }
        else
        {
            values.setState(Indices::firstPhaseOnly);
            values[Indices::switchIdx] = 0.0;
        }

#if !ISOTHERMAL
        values[Indices::temperatureIdx] = outletTemperature_;
#endif

        return values;
    }

    //!  Evaluate the initial invasion state of a pore throat
    bool initialInvasionState(const Element& element) const
    { return false; }

    // \}

    //! Loop over the scv in the domain to calculate the sum volume of inner inlet pores
    void calculateSumInletVolume()
    {
        sumInletPoresVolume_ = 0.0;
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bind(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                if (isInletPore_(scv))
                        sumInletPoresVolume_ += this->gridGeometry().poreVolume(scv.dofIndex())/this->gridGeometry().coordinationNumber(scv.dofIndex());
            }
        }
    }

private:

    bool isInletPore_(const SubControlVolume& scv) const
    {
        return isInletPore_(scv.dofIndex());
    }

    bool isInletPore_(const std::size_t dofIdxGlobal) const
    {
        return this->gridGeometry().poreLabel(dofIdxGlobal) == 4;
    }

    bool isOutletPore_(const SubControlVolume& scv) const
    {
        return this->gridGeometry().poreLabel(scv.dofIndex()) == 3;
    }

    int vtpOutputFrequency_;
    bool useFixedPressureAndSaturationBoundary_;
    Scalar pc_;
    Scalar source_;
    Scalar inletPressure_;
    Scalar outletPressure_;
    Scalar sumInletPoresVolume_;
#if !ISOTHERMAL
    Scalar inletTemperature_;
    Scalar outletTemperature_;
#endif
};
} //end namespace Dumux

#endif
