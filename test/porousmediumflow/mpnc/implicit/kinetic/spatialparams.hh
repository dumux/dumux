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
 * \ingroup MPNCTests
 * \brief spatialparameters for the kinetic test-case of the mpnc model. "Poor-mans" coupling of free-flow and porous medium.
 *
 */
#ifndef DUMUX_EVAPORATION_ATMOSPHERE_SPATIALPARAMS_HH
#define DUMUX_EVAPORATION_ATMOSPHERE_SPATIALPARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fvnonequilibrium.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchtenoftemperature.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/common/parameters.hh>

// material laws for interfacial area
#include <dumux/material/fluidmatrixinteractions/2pia/efftoabslawia.hh>
#include <dumux/material/fluidmatrixinteractions/2pia/awnsurfacepolynomial2ndorder.hh>
#include <dumux/material/fluidmatrixinteractions/2pia/awnsurfacepcmaxfct.hh>
#include <dumux/material/fluidmatrixinteractions/2pia/awnsurfaceexpswpcto3.hh>

namespace Dumux {

/**
 * \brief Definition of the spatial parameters for the evaporation atmosphere Problem (using a "poor man's coupling")
 */
template<class FVGridGeometry, class Scalar>
class EvaporationAtmosphereSpatialParams
: public FVNonEquilibriumSpatialParams<FVGridGeometry, Scalar,
                                       EvaporationAtmosphereSpatialParams<FVGridGeometry, Scalar>>
{
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ThisType = EvaporationAtmosphereSpatialParams<FVGridGeometry, Scalar>;
    using ParentType = FVNonEquilibriumSpatialParams<FVGridGeometry, Scalar, ThisType>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr auto dimWorld = GridView::dimensionworld;
public:
    //! export the type used for the permeability
    using PermeabilityType = Scalar;
    //! export the material law type used
    using MaterialLaw = EffToAbsLaw<RegularizedBrooksCorey<Scalar>>;
    //! convenience aliases of the law parameters
    using MaterialLawParams = typename MaterialLaw::Params;

    //! export the types used for interfacial area calculations
    using EffectiveIALawAws = AwnSurfacePolynomial2ndOrder<Scalar>;
    using EffectiveIALawAwn = AwnSurfacePcMaxFct<Scalar>;
    using EffectiveIALawAns = AwnSurfaceExpSwPcTo3<Scalar>;
    using AwnSurface = EffToAbsLawIA<EffectiveIALawAwn, MaterialLawParams>;
    using AwsSurface = EffToAbsLawIA<EffectiveIALawAws, MaterialLawParams>;
    using AnsSurface = EffToAbsLawIA<EffectiveIALawAns, MaterialLawParams>;

    using AwnSurfaceParams = typename AwnSurface::Params;
    using AwsSurfaceParams = typename AwsSurface::Params;
    using AnsSurfaceParams = typename AnsSurface::Params;

    EvaporationAtmosphereSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        heightPM_               = getParam<std::vector<Scalar>>("Grid.Positions1")[1];
        heightDomain_           = getParam<std::vector<Scalar>>("Grid.Positions1")[2];

        porosityPM_                 = getParam<Scalar>("SpatialParams.PorousMedium.porosity");
        intrinsicPermeabilityPM_    = getParam<Scalar>("SpatialParams.PorousMedium.permeability");

        porosityFF_                 = getParam<Scalar>("SpatialParams.FreeFlow.porosity");
        intrinsicPermeabilityFF_    = getParam<Scalar>("SpatialParams.FreeFlow.permeability");

        aWettingNonWettingA1_ = getParam<Scalar>("SpatialParams.soil.aWettingNonWettingA1");
        aWettingNonWettingA2_ = getParam<Scalar>("SpatialParams.soil.aWettingNonWettingA2");
        aWettingNonWettingA3_ = getParam<Scalar>("SpatialParams.soil.aWettingNonWettingA3");

        aNonWettingSolidA1_ = getParam<Scalar>("SpatialParams.soil.aNonWettingSolidA1");
        aNonWettingSolidA2_ = getParam<Scalar>("SpatialParams.soil.aNonWettingSolidA2");
        aNonWettingSolidA3_ = getParam<Scalar>("SpatialParams.soil.aNonWettingSolidA3");

        BCPd_           = getParam<Scalar>("SpatialParams.soil.BCPd");
        BClambda_       = getParam<Scalar>("SpatialParams.soil.BClambda");
        Swr_            = getParam<Scalar>("SpatialParams.soil.Swr");
        Snr_            = getParam<Scalar>("SpatialParams.soil.Snr");

        characteristicLengthFF_   = getParam<Scalar>("SpatialParams.FreeFlow.meanPoreSize");
        characteristicLengthPM_   = getParam<Scalar>("SpatialParams.PorousMedium.meanPoreSize");

        factorEnergyTransfer_ = getParam<Scalar>("SpatialParams.PorousMedium.factorEnergyTransfer");
        factorMassTransfer_ = getParam<Scalar>("SpatialParams.PorousMedium.factorMassTransfer");

        // residual saturations
        materialParamsFF_.setSwr(0.0);
        materialParamsFF_.setSnr(0.00);

        materialParamsPM_.setSwr(Swr_);
        materialParamsPM_.setSnr(Snr_);

        // pc / kr parameters
        materialParamsPM_.setLambda(BClambda_);
        materialParamsPM_.setPe(BCPd_);

        // for making pc == 0 in the FF
        materialParamsFF_.setLambda(42.);
        materialParamsFF_.setPe(0.);

        // determine maximum capillary pressure for wetting-nonwetting surface
        using TwoPLaw = EffToAbsLaw<RegularizedBrooksCorey<Scalar>>;
        pcMax_ = TwoPLaw::pc(materialParamsPM_, /*sw = */0.0);
        aWettingNonWettingSurfaceParams_.setPcMax(pcMax_);

        // wetting-non wetting: surface which goes to zero on the edges, but is a polynomial
        aWettingNonWettingSurfaceParams_.setA1(aWettingNonWettingA1_);
        aWettingNonWettingSurfaceParams_.setA2(aWettingNonWettingA2_);
        aWettingNonWettingSurfaceParams_.setA3(aWettingNonWettingA3_);

        // non-wetting-solid
        aNonWettingSolidSurfaceParams_.setA1(aNonWettingSolidA1_);
        aNonWettingSolidSurfaceParams_.setA2(aNonWettingSolidA2_);
        aNonWettingSolidSurfaceParams_.setA3(aNonWettingSolidA3_);

        // dummys for free flow: no interface where there is only one phase
        aWettingNonWettingSurfaceParamsFreeFlow_.setA1(0.);
        aWettingNonWettingSurfaceParamsFreeFlow_.setA2(0.);
        aWettingNonWettingSurfaceParamsFreeFlow_.setA3(0.);
        aWettingNonWettingSurfaceParamsFreeFlow_.setPcMax(42.); // not needed because it is anyways zero; silencing valgrind

        // dummys for free flow: no interface where there is only one phase
        aNonWettingSolidSurfaceParamsFreeFlow_.setA1(0.);
        aNonWettingSolidSurfaceParamsFreeFlow_.setA2(0.);
        aNonWettingSolidSurfaceParamsFreeFlow_.setA3(0.);
    }

    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        const  auto & globalPos = scv.dofPosition();
        if (inFF_(globalPos))
            return intrinsicPermeabilityFF_ ;
        else if (inPM_(globalPos))
            return intrinsicPermeabilityPM_ ;
        else
            DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);
    }

    /*!
     * \brief Function for defining the porosity.
     *        That is possibly solution dependent.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return the porosity
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        const auto& globalPos =  scv.dofPosition();

        if (inFF_(globalPos))
            return porosityFF_ ;
        else if (inPM_(globalPos))
            return porosityPM_ ;
        else
            DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);
    }

    template<class ElementSolution>
    const MaterialLawParams& materialLawParams(const Element& element,
                                               const SubControlVolume& scv,
                                               const ElementSolution& elemSol) const
    { return materialLawParamsAtPos(scv.dofPosition()); }

    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        if (inFF_(globalPos))
            return materialParamsFF_;
        else if (inPM_(globalPos))
            return materialParamsPM_;
        else DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);
    }

    /*!\brief Return a reference to the container object for the
     *        parametrization of the surface between wetting and non-Wetting phase.
     *
     *        The position is determined based on the coordinate of
     *        the vertex belonging to the considered sub control volume.
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub-control volume */
    template<class ElementSolution>
    const AwnSurfaceParams& aWettingNonWettingSurfaceParams(const Element &element,
                                                            const SubControlVolume &scv,
                                                            const ElementSolution &elemSol) const
    {
        const auto& globalPos =  scv.dofPosition();
        if (inFF_(globalPos) )
            return aWettingNonWettingSurfaceParamsFreeFlow_  ;
        else if (inPM_(globalPos))
            return aWettingNonWettingSurfaceParams_ ;
        else DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);
     }

    /*!\brief Return a reference to the container object for the
     *        parametrization of the surface between non-Wetting and solid phase.
     *
     *        The position is determined based on the coordinate of
     *        the vertex belonging to the considered sub control volume.
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub-control volume */
    template<class ElementSolution>
    const AnsSurfaceParams& aNonWettingSolidSurfaceParams(const Element &element,
                                                          const SubControlVolume &scv,
                                                          const ElementSolution &elemSol) const
    {
        const auto& globalPos =  scv.dofPosition();
        if (inFF_(globalPos) )
            return aNonWettingSolidSurfaceParamsFreeFlow_  ;
        else if (inPM_(globalPos))
            return aNonWettingSolidSurfaceParams_ ;
        else DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);
     }

    /*!\brief Return the maximum capillary pressure for the given pc-Sw curve
     *
     *        Of course physically there is no such thing as a maximum capillary pressure.
     *        The parametrization (VG/BC) goes to infinity and physically there is only one pressure
     *        for single phase conditions.
     *        Here, this is used for fitting the interfacial area surface: the capillary pressure,
     *        where the interfacial area is zero.
     *        Technically this value is obtained as the capillary pressure of saturation zero.
     *        This value of course only exists for the case of a regularized pc-Sw relation.
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub-control volume */
    template<class ElementSolution>
    const Scalar pcMax(const Element &element,
                       const SubControlVolume &scv,
                       const ElementSolution &elemSol) const
    { return aWettingNonWettingSurfaceParams_.pcMax() ; }

    /*!\brief Return the characteristic length for the mass transfer.
     * \param globalPos The position in global coordinates.*/
    const Scalar characteristicLengthAtPos(const  GlobalPosition & globalPos) const
    {
        if (inFF_(globalPos) )
            return characteristicLengthFF_ ;
        else if (inPM_(globalPos))
            return characteristicLengthPM_ ;
        else DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);
    }

    /*!\brief Return the pre factor the the energy transfer
     * \param globalPos The position in global coordinates.*/
    const Scalar factorEnergyTransferAtPos(const  GlobalPosition & globalPos) const
    {
        if (inFF_(globalPos) )
            return factorEnergyTransfer_ ;
        else if (inPM_(globalPos))
            return factorEnergyTransfer_ ;
        else DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);
    }

    /*!\brief Return the pre factor the the mass transfer
     * \param globalPos The position in global coordinates.*/
    const Scalar factorMassTransferAtPos(const  GlobalPosition & globalPos) const
    {
        if (inFF_(globalPos) )
            return factorMassTransfer_ ;
        else if (inPM_(globalPos))
            return factorMassTransfer_ ;
        else DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The global position
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    {
        return FluidSystem::phase0Idx;
    }

    /*!\brief Give back whether the tested position (input) is a specific region (porous medium part) in the domain
     *
     * This setting ensures, that the boundary between the two domains has porous medium properties.
     * This is desirable, because I want to observe the boundary of the porous domain.
     * However, the position has to be the coordinate of the vertex and not the integration point
     * of the boundary flux. If the position is the ip of the neumann flux this leads to a situation
     * where the vertex belongs to porous medium and there is nonetheless injection over the boundary.
     * That does not work.
     * -> be careful with neumannAtPos
     */
    bool inPM_(const GlobalPosition & globalPos) const
    { return ( (globalPos[dimWorld-1] > 0. - eps_) and (globalPos[dimWorld-1] < (heightPM_ + eps_) ) );   }

    /*!
     * \brief Give back whether the tested position (input) is a specific region (above PM / "free flow") in the domain
     *
     * This setting ensures, that the boundary between the two domains has porous medium properties.
     * This is desirable, because I want to observe the boundary of the porous domain.
     * However, the position has to be the coordinate of the vertex and not the integration point
     * of the boundary flux. If the position is the ip of the neumann flux this leads to a situation
     * where the vertex belongs to porous medium and there is nonetheless injection over the boundary.
     * That does not work.
     * -> be careful with neumannAtPos
     */
    bool inFF_(const GlobalPosition & globalPos) const
    { return ( (globalPos[dimWorld-1] < heightDomain_ + eps_) and (globalPos[dimWorld-1] > (heightPM_ + eps_) ) );   }

    /*! \brief access function for the depth / height of the porous medium */
    const Scalar heightPM() const
    { return heightPM_; }

private:
    static constexpr Scalar eps_  = 1e-6;
    Scalar heightDomain_ ;

    AwnSurfaceParams aWettingNonWettingSurfaceParams_;
    AnsSurfaceParams aNonWettingSolidSurfaceParams_ ;
    AwnSurfaceParams aWettingNonWettingSurfaceParamsFreeFlow_;
    AnsSurfaceParams aNonWettingSolidSurfaceParamsFreeFlow_ ;

    Scalar pcMax_ ;

    // Porous Medium Domain
    Scalar intrinsicPermeabilityPM_ ;
    Scalar porosityPM_ ;
    Scalar heightPM_ ;
    Scalar factorEnergyTransfer_ ;
    Scalar factorMassTransfer_ ;
    Scalar characteristicLengthPM_ ;
    MaterialLawParams materialParamsPM_ ;

    // Free Flow Domain
    Scalar porosityFF_ ;
    Scalar intrinsicPermeabilityFF_ ;
    Scalar characteristicLengthFF_ ;
    MaterialLawParams materialParamsFF_ ;

    // interfacial area parameters
    Scalar aWettingNonWettingA1_ ;
    Scalar aWettingNonWettingA2_ ;
    Scalar aWettingNonWettingA3_ ;

    Scalar aNonWettingSolidA1_;
    Scalar aNonWettingSolidA2_;
    Scalar aNonWettingSolidA3_;

    // capillary pressures parameters
    Scalar BCPd_ ;
    Scalar BClambda_ ;
    Scalar Swr_ ;
    Scalar Snr_ ;
    std::vector<Scalar> gridVector_;
};

}

#endif // GUARDIAN
