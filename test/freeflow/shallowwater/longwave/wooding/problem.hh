// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief Rainfall runoff on Wooding's V-catchment with the long-wave models
 */
#ifndef DUMUX_TEST_FREEFLOW_SHALLOWWATER_LONGWAVE_WOODING_PROBLEM_HH
#define DUMUX_TEST_FREEFLOW_SHALLOWWATER_LONGWAVE_WOODING_PROBLEM_HH

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>

#include <dumux/freeflow/shallowwater/longwave/discharge.hh>

namespace Dumux {

/*!
 * \ingroup ShallowWaterTests
 * \brief Rainfall runoff on Wooding's V-catchment with the long-wave models
 *
 * The V-shaped catchment of Wooding (1965): two planes of 800 m by 1000 m drain sideways into a
 * channel of 20 m width, which drains along its length towards the outlet. Rain falls at a
 * constant rate on the whole catchment for a given time.
 *
 * The channel is resolved in two dimensions. Water leaves the catchment only where the channel
 * meets the lower boundary, under normal flow; the rest of the boundary is closed. Once the
 * catchment has reached equilibrium, the outflow equals the rainfall on the catchment.
 */
template<class TypeTag>
class WoodingProblem : public FVProblemWithSpatialParams<TypeTag>
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using BoundaryTypes = Dumux::BoundaryTypes<ModelTraits::numEq()>;
    using Indices = typename ModelTraits::Indices;

    static constexpr int dimWorld = GridGeometry::GridView::dimensionworld;

public:
    WoodingProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , rainFallRate_(getParam<Scalar>("Problem.RainFallRate"))
    , outletMaxY_(getParam<Scalar>("Problem.OutletMaxY"))
    {}

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition&) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Normal flow at the outlet, where the free surface parallels the bed; closed elsewhere
     */
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector flux(0.0);
        if (!isOutlet(scvf.ipGlobal()))
            return flux;

        const auto& fluxVarCache = elemFluxVarsCache[scvf];
        Dune::FieldVector<Scalar, dimWorld> gradZ(0.0);
        for (const auto& scv : scvs(fvGeometry))
            gradZ.axpy(elemVolVars[scv].bedSurface(), fluxVarCache.gradN(scv.indexInElement()));

        flux[Indices::massBalanceIdx] = LongWave::normalFlowDischarge(
            gradZ, scvf.unitOuterNormal(),
            elemVolVars[scvf.insideScvIdx()].waterDepth(),
            this->spatialParams().manningN(element)
        );
        return flux;
    }

    NumEqVector sourceAtPos(const GlobalPosition&) const
    {
        NumEqVector source(0.0);
        source[Indices::massBalanceIdx] = rainFallRate();
        return source;
    }

    PrimaryVariables initialAtPos(const GlobalPosition&) const
    { return PrimaryVariables(0.0); }

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
    bool rainActive_ = false;
};

} // end namespace Dumux

#endif
