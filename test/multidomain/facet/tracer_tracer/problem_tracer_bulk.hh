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
 * \ingroup FacetTests
 * \brief The problem for the bulk domain of the tracer facet coupling test.
 */

#ifndef DUMUX_TEST_TPFAFACETCOUPLING_TRACER_BULK_PROBLEM_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_TRACER_BULK_PROBLEM_HH

#include <dune/alugrid/grid.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/multidomain/facet/box/properties.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/properties.hh>
#include <dumux/multidomain/facet/cellcentered/mpfa/properties.hh>

#include <dumux/porousmediumflow/tracer/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/fluidsystems/base.hh>

#include "spatialparams_tracer.hh"
#include "tracerfluidsystem.hh"
#include "tracermodeltraits.hh"

namespace Dumux {

template <class TypeTag>
class TracerBulkProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct TracerTestBulk { using InheritsFrom = std::tuple<Tracer>; };

// define the type tags
struct TracerBulkTpfa { using InheritsFrom = std::tuple<CCTpfaFacetCouplingModel, TracerTestBulk>; };
struct TracerBulkMpfa { using InheritsFrom = std::tuple<CCMpfaFacetCouplingModel, TracerTestBulk>; };
struct TracerBulkBox { using InheritsFrom = std::tuple<BoxFacetCouplingModel, TracerTestBulk>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TracerTestBulk> { using type = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>; };

//! Overwrite the advection type property
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::TracerBulkTpfa> { using type = StationaryVelocityField<GetPropType<TypeTag, Properties::Scalar>>; };
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::TracerBulkBox> { using type = StationaryVelocityField<GetPropType<TypeTag, Properties::Scalar>>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TracerTestBulk> { using type = TracerBulkProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TracerTestBulk>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TracerSpatialParams<GridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::TracerTestBulk> { static constexpr bool value = false; };

//! set the model traits (with disabled diffusion)
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::TracerTestBulk>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = TracerTestModelTraits<FluidSystem::numComponents, getPropValue<TypeTag, Properties::UseMoles>()>;
};

// use the test-specific fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TracerTestBulk> { using type = TracerFluidSystem<TypeTag>; };
} // end namespace Properties


/*!
 * \ingroup FacetTests
 * \brief The problem for the bulk domain of the tracer facet coupling test.
 */
template <class TypeTag>
class TracerBulkProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    //! property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using typename ParentType::SpatialParams;

    TracerBulkProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                      std::shared_ptr<SpatialParams> spatialParams,
                      std::shared_ptr<CouplingManager> couplingManager,
                      const std::string& paramGroup = "")
    : ParentType(gridGeometry, spatialParams, paramGroup)
    , couplingManagerPtr_(couplingManager)
    , initialMassFraction_(getParamFromGroup<Scalar>(paramGroup, "Problem.ContaminationMassFraction"))
    {
        // stating in the console whether mole or mass fractions are used
        const auto problemName = getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        std::cout<< "problem " << problemName << " uses " << (useMoles ? "mole" : "mass") << " fractions" << '\n';
        problemName_  =  getParamFromGroup<std::string>(this->paramGroup(), "Vtk.OutputName") + "_" + problemName;
    }

    //! The problem name.
    const std::string& name() const
    { return problemName_; }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    //! Specifies the type of interior boundary condition at a given position.
    BoundaryTypes interiorBoundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        if (isInContaminatedRegion_(globalPos))
            return PrimaryVariables(initialMassFraction_);
        return PrimaryVariables(0.0);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The sub control volume face
     *
     * Negative values mean influx.
     * E.g. for the mass balance that would the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    template<class ElementVolumeVariables, class ElementFluxVarsCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVarsCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        // get the volume flux on this segment
        const auto flux = this->spatialParams().volumeFlux(element, fvGeometry, elemVolVars, scvf);
        if (flux > 0.0)
        {
            const auto& insideVolVars = elemVolVars[fvGeometry.scv(scvf.insideScvIdx())];
            const auto tracerFlux = insideVolVars.massFraction(/*phaseIdx*/0, /*compIdx*/0)*flux;
            return NumEqVector(tracerFlux);
        }

        return NumEqVector(0.0);
    }

    //! Returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

private:
    bool isInContaminatedRegion_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > 0.75 && globalPos[0] < 1.3
               && globalPos[1] > 5.3 && globalPos[1] < 5.8;
    }

    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    Scalar initialMassFraction_;
    std::string problemName_;
};

} // end namespace Dumux

#endif
