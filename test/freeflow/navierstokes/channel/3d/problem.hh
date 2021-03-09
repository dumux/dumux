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
 * \ingroup NavierStokesTests
 * \brief Channel flow test for the staggered grid (Navier-)Stokes model.
 *
 * The channel is either modeled in 3D or in 2D, using an additional wall friction term
 * to mimic the 3D behavior of the flow.
 */

#ifndef DUMUX_3D_CHANNEL_PROBLEM_HH
#define DUMUX_3D_CHANNEL_PROBLEM_HH

#include <dune/common/float_cmp.hh>
#include <dune/grid/yaspgrid.hh>

#ifndef GRID_DIM
#define GRID_DIM 3
#endif

#if HAVE_DUNE_SUBGRID && GRID_DIM == 3
#include <dune/subgrid/subgrid.hh>
#endif

#include <dumux/discretization/staggered/freeflow/properties.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/model.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

namespace Dumux {

template <class TypeTag>
class ThreeDChannelTestProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct ThreeDChannelTest { using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ThreeDChannelTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::ThreeDChannelTest>
{
    static constexpr int dim = GRID_DIM;

    using HostGrid = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, dim> >;

#if HAVE_DUNE_SUBGRID && GRID_DIM == 3
    using type = Dune::SubGrid<HostGrid::dimension, HostGrid>;
#else
    using type = HostGrid;
#endif
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::ThreeDChannelTest> { using type = ThreeDChannelTestProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::ThreeDChannelTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::ThreeDChannelTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::ThreeDChannelTest> { static constexpr bool value = true; };
} // end namespace Properties

/*!
 * \brief Test problem for the one-phase (Navier-) Stokes model in a 3D or pseudo 3D channel.
 *
 * Flow from left to right in a three-dimensional channel is considered. At the inlet (left)
 * and outlet (right) fixed values for pressure are set.
 * The channel is confined by solid walls at all other sides of the domain which corresponds
 * to no-slip/no-flow conditions.
 * The value of an analytical solution for the given flow configuration is furthermore provided.
 * For sake of efficiency, the 3D problem can be reduced to a two-dimensional one by including
 * an additional wall friction term to the momentum balance (Flekkoy et al., 1995 \cite flekkoy1995a).
 */
template <class TypeTag>
class ThreeDChannelTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryTypes = Dumux::NavierStokesBoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Element = typename GridView::template Codim<0>::Entity;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CellCenterPrimaryVariables = GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>;
    using FacePrimaryVariables = GetPropType<TypeTag, Properties::FacePrimaryVariables>;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

    static constexpr bool enablePseudoThreeDWallFriction = !(GRID_DIM == 3);

public:
    ThreeDChannelTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        deltaP_ = getParam<Scalar>("Problem.DeltaP");
        height_ = getParam<Scalar>("Problem.Height");
        rho_ = getParam<Scalar>("Component.LiquidDensity");
        nu_ = getParam<Scalar>("Component.LiquidKinematicViscosity");

        if(dim == 3 && !Dune::FloatCmp::eq(height_, this->gridGeometry().bBoxMax()[2]))
            DUNE_THROW(Dune::InvalidStateException, "z-dimension must equal height");

        if(enablePseudoThreeDWallFriction)
            extrusionFactor_ = 2.0/3.0 * height_;
        else
            extrusionFactor_ = 1.0;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C


    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub-control volume face.
     */
    using ParentType::source;
    template<class ElementVolumeVariables, class ElementFaceVariables>
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const ElementFaceVariables& elemFaceVars,
                       const SubControlVolumeFace &scvf) const
    {
        auto source = NumEqVector(0.0);

#if GRID_DIM != 3
            static const Scalar height = getParam<Scalar>("Problem.Height");
            static const Scalar factor = getParam<Scalar>("Problem.PseudoWallFractionFactor", 8.0);
            source[scvf.directionIndex()] = this->pseudo3DWallFriction(scvf, elemVolVars, elemFaceVars, height, factor);
#endif

        return source;
    }

    Scalar extrusionFactorAtPos(const GlobalPosition& pos) const
    { return extrusionFactor_; }


    // \}
    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        // set a fixed pressure at the inlet and outlet
        if (isOutlet_(globalPos) || isInlet_(globalPos))
            values.setDirichlet(Indices::pressureIdx);
        else
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
            if(dim == 3)
                values.setDirichlet(Indices::velocityZIdx);
        }

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        if(isInlet_(globalPos))
            values[Indices::pressureIdx] = 1e5 + deltaP_;

        if(isOutlet_(globalPos))
            values[Indices::pressureIdx] = 1e5;

        return values;
    }

    // \}

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        return values;
    }

    //! Returns the analytical solution for the flux through the rectangular channel
    Scalar analyticalFlux() const
    {
        const Scalar h = height_;
        const Scalar w = this->gridGeometry().bBoxMax()[1];
        const Scalar L = this->gridGeometry().bBoxMax()[0];

        const Scalar mu = nu_*rho_;

        return h*h*h * w * deltaP_ / (12*mu*L) * (1.0 - 0.630 * h/w);
    }

private:

    bool isInlet_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < eps_;
    }

    bool isOutlet_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
    }

    static constexpr Scalar eps_=1e-6;
    Scalar deltaP_;
    Scalar extrusionFactor_;
    Scalar height_;
    Scalar rho_;
    Scalar nu_;
};
} // end namespace Dumux

#endif
