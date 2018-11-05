// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief Channel flow test for the staggered grid (Navier-)Stokes model.
 *
 * The channel is either modeled in 3D or in 2D, using an additional wall friction term
 * to mimic the 3D behavior of the flow.
 *
 */
#ifndef DUMUX_3D_CHANNEL_PROBLEM_HH
#define DUMUX_3D_CHANNEL_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/model.hh>
#include <dune/common/float_cmp.hh>

#ifndef DIM_3D
#define DIM_3D 0
#endif

namespace Dumux
{


template <class TypeTag>
class ThreeDChannelTestProblem;

namespace Properties
{
NEW_TYPE_TAG(ThreeDChannelTest, INHERITS_FROM(StaggeredFreeFlowModel, NavierStokes));

// the fluid system
SET_PROP(ThreeDChannelTest, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
#if DIM_3D
SET_TYPE_PROP(ThreeDChannelTest, Grid, Dune::YaspGrid<3>);
#else
SET_TYPE_PROP(ThreeDChannelTest, Grid, Dune::YaspGrid<2>);
#endif

// Set the problem property
SET_TYPE_PROP(ThreeDChannelTest, Problem, ThreeDChannelTestProblem<TypeTag> );

SET_BOOL_PROP(ThreeDChannelTest, EnableFVGridGeometryCache, true);
SET_BOOL_PROP(ThreeDChannelTest, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(ThreeDChannelTest, EnableGridVolumeVariablesCache, true);
}

/*!
 * \brief  Test problem for the one-phase model:
   \todo doc me!
 */
template <class TypeTag>
class ThreeDChannelTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using Element = typename GridView::template Codim<0>::Entity;

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);

    static constexpr bool enablePseudoThreeDWallFriction = !DIM_3D;

public:
    ThreeDChannelTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry), eps_(1e-6)
    {
        deltaP_ = getParam<Scalar>("Problem.DeltaP");
        height_ = getParam<Scalar>("Problem.Height");
        rho_ = getParam<Scalar>("Component.LiquidDensity");
        nu_ = getParam<Scalar>("Component.LiquidKinematicViscosity");

        if(dim == 3 && !Dune::FloatCmp::eq(height_, this->fvGridGeometry().bBoxMax()[2]))
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
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C


    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume face.
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

#if !DIM_3D
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
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
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
     * \brief Evaluate the initial value for a control volume.
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
        const Scalar w = this->fvGridGeometry().bBoxMax()[1];
        const Scalar L = this->fvGridGeometry().bBoxMax()[0];

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
        return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_;
    }

    Scalar eps_;
    Scalar deltaP_;
    Scalar extrusionFactor_;
    Scalar height_;
    Scalar rho_;
    Scalar nu_;
};
} //end namespace

#endif
