// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Channel flow test for the staggered grid (Navier-)Stokes model.
 *
 * The channel is either modeled in 3D or in 2D, using an additional wall friction term
 * to mimic the 3D behavior of the flow.
 */

#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_3D_CHANNEL_PROBLEM_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_3D_CHANNEL_PROBLEM_HH

#include <dune/common/float_cmp.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
#include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>
#include <dumux/freeflow/navierstokes/mass/1p/advectiveflux.hh>

namespace Dumux {

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
template <class TypeTag, class BaseProblem>
class ThreeDChannelTestProblem : public BaseProblem
{
    using ParentType = BaseProblem;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    static constexpr int dim = GridGeometry::GridView::dimension;
    static constexpr bool enablePseudoThreeDWallFriction = dim != 3;

public:
    ThreeDChannelTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        deltaP_ = getParam<Scalar>("Problem.DeltaP");
        rho_ = getParam<Scalar>("Component.LiquidDensity");
        nu_ = getParam<Scalar>("Component.LiquidKinematicViscosity");

        height_ = getParam<Scalar>("Problem.Height");
        if(dim == 3 && !Dune::FloatCmp::eq(height_, this->gridGeometry().bBoxMax()[2]))
            DUNE_THROW(Dune::InvalidStateException, "z-dimension must equal height");
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub-control volume
     */
    template<class ElementVolumeVariables>
    Sources source(const Element& element,
                   const FVElementGeometry& fvGeometry,
                   const ElementVolumeVariables& elemVolVars,
                   const SubControlVolume& scv) const
    {
        auto source = Sources(0.0);

        if constexpr (ParentType::isMomentumProblem() && enablePseudoThreeDWallFriction)
        {
            static const Scalar height = getParam<Scalar>("Problem.Height");
            static const Scalar factor = getParam<Scalar>("Problem.PseudoWallFractionFactor", 8.0);
            source[scv.dofAxis()] = this->pseudo3DWallFriction(element, fvGeometry, elemVolVars, scv, height, factor);
        }

        return source;
    }

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

        if constexpr (ParentType::isMomentumProblem())
        {
            values.setAllDirichlet();
            if (isOutlet_(globalPos) || isInlet_(globalPos))
                values.setAllNeumann();
        }
        else
            values.setNeumann(Indices::conti0EqIdx);

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    DirichletValues dirichletAtPos(const GlobalPosition &globalPos) const
    {
        // no-flow/no-slip
        return DirichletValues(0.0);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    BoundaryFluxes neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const SubControlVolumeFace& scvf) const
    {
        BoundaryFluxes values(0.0);
        const auto& globalPos = scvf.ipGlobal();

        if constexpr (ParentType::isMomentumProblem())
        {
            const auto p = isInlet_(globalPos) ? 1e5 + deltaP_ : 1e5;
            values = NavierStokesMomentumBoundaryFluxHelper::fixedPressureMomentumFlux(
                *this, fvGeometry, scvf, elemVolVars, elemFluxVarsCache, p
            );
        }
        else
        {
            values = NavierStokesScalarBoundaryFluxHelper<AdvectiveFlux<ModelTraits>>::scalarOutflowFlux(
                *this, element, fvGeometry, scvf, elemVolVars
            );
        }

        return values;
    }

    // \}

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
    { return globalPos[0] < eps_; }

    bool isOutlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    static constexpr Scalar eps_=1e-6;
    Scalar deltaP_;
    Scalar height_;
    Scalar rho_;
    Scalar nu_;
};

} // end namespace Dumux

#endif
