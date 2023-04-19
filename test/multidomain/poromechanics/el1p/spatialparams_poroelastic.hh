// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Definition of the spatial parameters for the poro-elastic problem.
 */
#ifndef DUMUX_POROELASTIC_SPATIAL_PARAMS_HH
#define DUMUX_POROELASTIC_SPATIAL_PARAMS_HH

#include <dumux/geomechanics/lameparams.hh>
#include <dumux/geomechanics/poroelastic/fvspatialparams.hh>
#include <dumux/material/fluidmatrixinteractions/porositydeformation.hh>

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief Definition of the spatial parameters for the poro-elastic
 *        sub-problem in the coupled poro-mechanical el1p problem.
 */
template<class Scalar, class GridGeometry, class CouplingManager>
class PoroElasticSpatialParams : public FVPoroElasticSpatialParams< GridGeometry,
                                                                    Scalar,
                                                                    PoroElasticSpatialParams<Scalar, GridGeometry, CouplingManager> >
{
    using ThisType = PoroElasticSpatialParams<Scalar, GridGeometry, CouplingManager>;
    using ParentType = FVPoroElasticSpatialParams<GridGeometry, Scalar, ThisType>;

    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
public:
    //! Export the type of the lame parameters
    using LameParams = Dumux::LameParams<Scalar>;

    PoroElasticSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                             std::shared_ptr<CouplingManager> couplingManagerPtr)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManagerPtr)
    , initPorosity_(getParam<Scalar>("SpatialParams.InitialPorosity"))
    {
        // Young's modulus [Pa]
       Scalar E = 6.e9;
       // Poisson's ratio [-]
       Scalar nu = 0.2;
       // Lame parameters [Pa]
       lameParams_.setLambda( (E * nu) / ((1 + nu)*(1 - 2 * nu)) );
       lameParams_.setMu( E / (2 * (1 + nu)) );
    }

    //! Defines the Lame parameters.
    const LameParams& lameParamsAtPos(const GlobalPosition& globalPos) const
    { return lameParams_; }

    //! Returns the porosity of the porous medium.
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const

    {
        return PorosityDeformation<Scalar>::evaluatePorosity(this->gridGeometry(), element, scv, elemSol, initPorosity_);
    }

    /*!
     * \brief Returns the effective fluid density.
     */
    Scalar effectiveFluidDensity(const Element& element,
                                 const SubControlVolume& scv) const
    {
        // get porous medium flow volume variables from coupling manager
        const auto& pmFlowVolVars = couplingManager().getPMFlowVolVars(element);
        return pmFlowVolVars.density();
    }

    /*!
     * \brief Returns the effective pore pressure.
     */
    template<class ElementVolumeVariables, class FluxVarsCache>
    Scalar effectivePorePressure(const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& elemVolVars,
                                 const FluxVarsCache& fluxVarsCache) const
    {
        // get porous medium flow volume variables from coupling manager
        const auto& pmFlowVolVars = couplingManager().getPMFlowVolVars(element);
        return pmFlowVolVars.pressure();
    }

    //! Returns the Biot coefficient of the porous medium.
    Scalar biotCoefficientAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    //! Returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    std::shared_ptr<const CouplingManager> couplingManager_;
    Scalar initPorosity_;
    LameParams lameParams_;
};

} // end namespace Dumux

#endif
