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
 * \ingroup TwoPTests
 * \brief The properties for the incompressible 2p-boxdfm test.
 */

#ifndef DUMUX_INCOMPRESSIBLE_TWOPBOXDFM_TEST_PROBLEM_HH
#define DUMUX_INCOMPRESSIBLE_TWOPBOXDFM_TEST_PROBLEM_HH

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/material/components/trichloroethene.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/porousmediumflow/boxdfm/model.hh>
#include <dumux/porousmediumflow/2p/incompressiblelocalresidual.hh>

#include "spatialparams.hh"

#ifndef GRIDTYPE
#define GRIDTYPE Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>;
#endif

namespace Dumux {

// forward declarations
template<class TypeTag> class TwoPTestProblem;

namespace Properties {

// we need to derive first from twop and then from the box-dfm Model
// because the flux variables cache type of TwoP is overwritten in BoxDfmModel
// Create new type tags
namespace TTag {
struct TwoPIncompressibleBoxDfm { using InheritsFrom = std::tuple<BoxDfmModel, TwoP>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPIncompressibleBoxDfm> { using type = GRIDTYPE; };

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPIncompressibleBoxDfm> { using type = TwoPTestProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TwoPIncompressibleBoxDfm>
{
private:
    using FVG = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = TwoPTestSpatialParams<FVG, Scalar>;
};

// the local residual containing the analytic derivative methods
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::TwoPIncompressibleBoxDfm> { using type = TwoPIncompressibleLocalResidual<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TwoPIncompressibleBoxDfm>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Trichloroethene<Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

// Enable caching
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TwoPIncompressibleBoxDfm> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::TwoPIncompressibleBoxDfm> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::TwoPIncompressibleBoxDfm> { static constexpr bool value = false; };

// Enable the box-interface solver
template<class TypeTag>
struct EnableBoxInterfaceSolver<TypeTag, TTag::TwoPIncompressibleBoxDfm> { static constexpr bool value = true; };
} // end namespace Properties

/*!
 * \ingroup TwoPTests
 * \brief The incompressible 2p-boxdfm test problem.
 */
template<class TypeTag>
class TwoPTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

    // some indices for convenience
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum
    {
        pressureH2OIdx = Indices::pressureIdx,
        saturationDNAPLIdx = Indices::saturationIdx,
        contiDNAPLEqIdx = Indices::conti0EqIdx + FluidSystem::comp1Idx,
        waterPhaseIdx = FluidSystem::phase0Idx,
        dnaplPhaseIdx = FluidSystem::phase1Idx
    };

    static constexpr Scalar eps_ = 1e-6;
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    TwoPTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    /*!
     * \brief Returns how much the domain is extruded at a given sub-control volume.
     *
     * Here, we extrude the fracture scvs by half the aperture
     */
    template<class ElementSolution>
    Scalar extrusionFactor(const Element& element,
                           const SubControlVolume& scv,
                           const ElementSolution& elemSol) const
    {
        // In the box-scheme, we compute fluxes etc element-wise,
        // thus per element we compute only half a fracture !!!
        static const Scalar aHalf = getParam<Scalar>("SpatialParams.FractureAperture")/2.0;

        if (scv.isOnFracture())
            return aHalf;

        return 1.0;
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        if (globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_ || globalPos[0] < 1e-6)
            values.setAllDirichlet();
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        auto values = initialAtPos(globalPos);
        if (globalPos[0] < 1e-6)
            values[saturationDNAPLIdx] = 0.5;
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

    /*!
     * \brief Evaluates the initial values for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;

        // pressure gradient from left to right
        values[pressureH2OIdx] = 2e5 - 1e5*globalPos[0]/this->gridGeometry().bBoxMax()[0];
        values[saturationDNAPLIdx] = 0;
        return values;
    }

    //! Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
    Scalar temperature() const { return 293.15; /* 10°C */ }
};

} // end namespace Dumux

#endif
