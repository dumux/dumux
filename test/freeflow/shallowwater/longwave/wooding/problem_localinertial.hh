// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief Rainfall runoff on Wooding's V-catchment with the local-inertial shallow water equations
 */
#ifndef DUMUX_TEST_FREEFLOW_SHALLOWWATER_LONGWAVE_WOODING_PROBLEM_LOCALINERTIAL_HH
#define DUMUX_TEST_FREEFLOW_SHALLOWWATER_LONGWAVE_WOODING_PROBLEM_LOCALINERTIAL_HH

#include <array>
#include <cmath>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/flux/shallowwater/riemannproblem.hh>
#include <dumux/freeflow/shallowwater/problem.hh>

namespace Dumux {

/*!
 * \ingroup ShallowWaterTests
 * \brief Rainfall runoff on Wooding's V-catchment with the local-inertial shallow water equations
 *
 * The same catchment as in Dumux::WoodingProblem, solved with local acceleration retained, as
 * a reference for the long-wave models. Water leaves the catchment only where the channel meets
 * the lower boundary, as a free overfall at critical depth; the rest of the boundary is a wall.
 */
template<class TypeTag>
class WoodingLocalInertialProblem : public ShallowWaterProblem<TypeTag>
{
    using ParentType = ShallowWaterProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using BoundaryTypes = Dumux::BoundaryTypes<ModelTraits::numEq()>;
    using Indices = typename ModelTraits::Indices;

public:
    WoodingLocalInertialProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , rainFallRate_(getParam<Scalar>("Problem.RainFallRate"))
    , outletMaxY_(getParam<Scalar>("Problem.OutletMaxY"))
    , initialWaterDepth_(getParam<Scalar>("Problem.InitialWaterDepth"))
    {}

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition&) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Boundary fluxes from a ghost state: the mirrored cell state at walls, and the
     *        critical-depth state of a free overfall at the outlet
     *
     * The local-inertial momentum term is a free-surface difference, so the ghost state that
     * determines the mass flux also drives it: it vanishes along the walls, where the ghost
     * mirrors the cell, and it is the drawdown of the overfall at the outlet.
     */
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        const auto& insideVolVars = elemVolVars[fvGeometry.scv(scvf.insideScvIdx())];
        const auto& nxy = scvf.unitOuterNormal();
        const auto gravity = this->spatialParams().gravity(scvf.center());

        std::array<Scalar, 3> outer{{insideVolVars.waterDepth(),
                                     -insideVolVars.velocity(0),
                                     -insideVolVars.velocity(1)}};
        if (isOutlet(scvf.ipGlobal()))
        {
            using std::hypot, std::sqrt, std::cbrt, std::max;
            const auto depth = insideVolVars.waterDepth();
            const auto speed = hypot(insideVolVars.velocity(0), insideVolVars.velocity(1));
            if (speed >= sqrt(gravity*max(depth, 1e-12)))
                outer = {{depth, insideVolVars.velocity(0), insideVolVars.velocity(1)}};
            else
            {
                const auto criticalDepth = cbrt(depth*speed*depth*speed/gravity);
                const auto rescale = sqrt(gravity*criticalDepth)/max(speed, 1e-12);
                outer = {{criticalDepth,
                          insideVolVars.velocity(0)*rescale,
                          insideVolVars.velocity(1)*rescale}};
            }
        }

        const auto riemannFlux = ShallowWater::riemannProblem(
            insideVolVars.waterDepth(), outer[0],
            insideVolVars.velocity(0), outer[1],
            insideVolVars.velocity(1), outer[2],
            insideVolVars.bedSurface(), insideVolVars.bedSurface(),
            gravity, nxy
        );

        const auto force = gravity*insideVolVars.waterDepth()*0.5*(outer[0] - insideVolVars.waterDepth());

        NumEqVector flux(0.0);
        flux[Indices::massBalanceIdx] = riemannFlux[0];
        flux[Indices::momentumXBalanceIdx] = force*nxy[0];
        flux[Indices::momentumYBalanceIdx] = force*nxy[1];
        return flux;
    }

    template<class ElementVolumeVariables>
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        NumEqVector source(0.0);
        source[Indices::massBalanceIdx] = rainFallRate();

        const auto& volVars = elemVolVars[scv];
        const auto stress = this->spatialParams().frictionLaw(element, scv).bottomShearStress(volVars);
        source[Indices::momentumXBalanceIdx] = -stress[0]/volVars.density();
        source[Indices::momentumYBalanceIdx] = -stress[1]/volVars.density();
        return source;
    }

    //! A thin film, since the velocity is undefined in a dry cell
    PrimaryVariables initialAtPos(const GlobalPosition&) const
    {
        PrimaryVariables values(0.0);
        values[Indices::waterdepthIdx] = initialWaterDepth_;
        return values;
    }

    void setRainActive(bool active)
    { rainActive_ = active; }

    Scalar rainFallRate() const
    { return rainActive_ ? rainFallRate_ : 0.0; }

    bool isOutlet(const GlobalPosition& globalPos) const
    { return globalPos[1] <= outletMaxY_ && this->spatialParams().catchment().inChannel(globalPos); }

    bool inChannel(const GlobalPosition& globalPos) const
    { return this->spatialParams().catchment().inChannel(globalPos); }

private:
    Scalar rainFallRate_;
    Scalar outletMaxY_;
    Scalar initialWaterDepth_;
    bool rainActive_ = false;
};

} // end namespace Dumux

#endif
