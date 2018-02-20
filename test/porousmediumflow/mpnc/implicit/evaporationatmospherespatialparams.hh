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

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchtenoftemperature.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/mp/2padapter.hh>
#include <dumux/material/fluidmatrixinteractions/mp/2poftadapter.hh>

// material laws for interfacial area
#include <dumux/material/fluidmatrixinteractions/2pia/efftoabslawia.hh>
#include <dumux/material/fluidmatrixinteractions/2pia/awnsurfacepolynomial2ndorder.hh>
#include <dumux/material/fluidmatrixinteractions/2pia/awnsurfacepolynomialedgezero2ndorder.hh>
#include <dumux/material/fluidmatrixinteractions/2pia/awnsurfaceexpfct.hh>
#include <dumux/material/fluidmatrixinteractions/2pia/awnsurfacepcmaxfct.hh>
#include <dumux/material/fluidmatrixinteractions/2pia/awnsurfaceexpswpcto3.hh>

#include <dune/common/parametertreeparser.hh>

namespace Dumux
{
/*!
 * \ingroup MPNCTests
 * \brief spatialparameters for the kinetic test-case of the mpnc model. "Poor-mans" coupling of free-flow and porous medium.
 *
 */
//forward declaration
template<class TypeTag>
class EvaporationAtmosphereSpatialParams;

namespace Properties
{
// The spatial params TypeTag
NEW_TYPE_TAG(EvaporationAtmosphereSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(EvaporationAtmosphereSpatialParams, SpatialParams, EvaporationAtmosphereSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(EvaporationAtmosphereSpatialParams, MaterialLaw)
{
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    enum {wPhaseIdx   = FluidSystem::wPhaseIdx};

private:
    // define the material law which is parameterized by effective
    // saturations
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    // select appropriate parametrization
#define linear      0 // pc(S_w)
#define RegVG       0 // pc(S_w)
#define RegBC       1 // pc(S_w)
#define RegVGofT    0 // pc(S_w, T) // adapt also the setting of the fluidstate
#define VG          0 // pc(S_w)

    // define the material law
    #if linear
        using EffectiveLaw = LinearMaterial<Scalar>;
    #endif

    #if RegVG
        using EffectiveLaw = RegularizedVanGenuchten<Scalar>;
    #endif

    #if RegBC
        using EffectiveLaw = RegularizedBrooksCorey<Scalar>;
    #endif

    #if RegVGofT
        using EffectiveLaw =  RegularizedVanGenuchtenOfTemperature<Scalar>;
    #endif

    #if VG
        using EffectiveLaw =  VanGenuchten<Scalar> ;
    #endif

    using TwoPMaterialLaw = EffToAbsLaw<EffectiveLaw>;
    public:
        using type = TwoPAdapter<wPhaseIdx, TwoPMaterialLaw>;
};



// Set the interfacial area relation: wetting -- non-wetting
SET_PROP(EvaporationAtmosphereSpatialParams, AwnSurface)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;
    using EffectiveIALaw = AwnSurfacePcMaxFct<Scalar>;
public:
    using type = EffToAbsLawIA<EffectiveIALaw, MaterialLawParams>;
};


// Set the interfacial area relation: wetting -- solid
SET_PROP(EvaporationAtmosphereSpatialParams, AwsSurface)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;
    using EffectiveIALaw = AwnSurfacePolynomial2ndOrder<Scalar>;
public:
    using type = EffToAbsLawIA<EffectiveIALaw, MaterialLawParams>;
};

// Set the interfacial area relation: non-wetting -- solid
SET_PROP(EvaporationAtmosphereSpatialParams, AnsSurface)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;
    using EffectiveIALaw = AwnSurfaceExpSwPcTo3<Scalar>;
public:
    using type = EffToAbsLawIA<EffectiveIALaw, MaterialLawParams>;
};

} // end namespace properties

/**
 * \brief Definition of the spatial parameters for the evaporation atmosphere Problem (using a "poor man's coupling")
 */
template<class TypeTag>
class EvaporationAtmosphereSpatialParams : public FVSpatialParams<TypeTag>
{
    using ParentType = FVSpatialParams<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);

public:

    using AwnSurface = typename GET_PROP_TYPE(TypeTag, AwnSurface);
    using AwnSurfaceParams = typename  AwnSurface::Params;
    using AwsSurface = typename GET_PROP_TYPE(TypeTag, AwsSurface);
    using AwsSurfaceParams = typename  AwsSurface::Params;
    using AnsSurface = typename GET_PROP_TYPE(TypeTag, AnsSurface);
    using AnsSurfaceParams = typename  AnsSurface::Params;
    using PermeabilityType = Scalar;

    EvaporationAtmosphereSpatialParams(const Problem &problem)
        : ParentType(problem)
    {
        heightPM_               = getParam<std::vector<Scalar>>("Grid.Positions1")[1];
        heightDomain_           = getParam<std::vector<Scalar>>("Grid.Positions1")[2];

        porosityPM_                 = getParam<Scalar>("SpatialParams.PorousMedium.porosity");
        intrinsicPermeabilityPM_    = getParam<Scalar>("SpatialParams.PorousMedium.permeability");

        porosityFF_                 = getParam<Scalar>("SpatialParams.FreeFlow.porosity");
        intrinsicPermeabilityFF_    = getParam<Scalar>("SpatialParams.FreeFlow.permeability");

        solidDensity_               = getParam<Scalar>("SpatialParams.soil.density");
        solidThermalConductivity_    = getParam<Scalar>("SpatialParams.soil.thermalConductivity");
        solidHeatCapacity_               = getParam<Scalar>("SpatialParams.soil.heatCapacity");

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

        {//scope it
            // capillary pressure parameters
            FluidState fluidState ;
            Scalar S[numPhases] ;
            const auto &materialParams =  materialParamsPM_;
                    Scalar capPress[numPhases];
            //set saturation to inital values, this needs to be done in order for the fluidState to tell me pc
            for (int phaseIdx = 0; phaseIdx < numPhases ; ++phaseIdx) {
                // set saturation to zero for getting pcmax
                S[phaseIdx] = 0. ;
                Scalar TInitial  = getParam<Scalar>("InitialConditions.TInitial");
                fluidState.setSaturation(phaseIdx, S[phaseIdx]);
                fluidState.setTemperature(phaseIdx,TInitial);
            }

            //obtain pc according to saturation
            MaterialLaw::capillaryPressures(capPress, materialParams, fluidState);
            using std::abs;
            pcMax_ = abs(capPress[0]);

            // set pressures from capillary pressures
            aWettingNonWettingSurfaceParams_.setPcMax(pcMax_);
        }

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

    ~EvaporationAtmosphereSpatialParams()
    {}


     PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolutionVector& elemSol) const
    {
        const  auto & globalPos =  scv.dofPosition();
        if (inFF_(globalPos) )
            return intrinsicPermeabilityFF_ ;
        else if (inPM_(globalPos))
            return intrinsicPermeabilityPM_ ;
        else
            DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);
    }

    /*! \brief Return the porosity \f$[-]\f$ of the soil
     *
     *        The position is determined based on the coordinate of
     *        the vertex belonging to the considered sub control volume.
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume  */
    Scalar porosity(const Element &element,
                    const SubControlVolume &scv,
                    const ElementSolutionVector &elemSol) const
    {
        const auto& globalPos =  scv.dofPosition();

        if (inFF_(globalPos) )
            return porosityFF_ ;
        else if (inPM_(globalPos))
            return porosityPM_ ;
        else
            DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);
    }

    const MaterialLawParams& materialLawParams(const Element& element,
                                               const SubControlVolume& scv,
                                               const ElementSolutionVector& elemSol) const
    {
        const auto& globalPos =  scv.dofPosition();
        if (inFF_(globalPos) )
            return materialParamsFF_ ;
        else if (inPM_(globalPos))
            return materialParamsPM_ ;
        else             DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);
    }

    /*!
     * \brief Return a reference to the material parameters of the material law.
     * \param globalPos The position in global coordinates. */
    const MaterialLawParams & materialLawParamsAtPos(const GlobalPosition & globalPos) const
    {
        if (inFF_(globalPos) )
            return materialParamsFF_ ;
        else if (inPM_(globalPos))
            return materialParamsPM_ ;
        else             DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);
    }

    /*!\brief Return a reference to the container object for the
     *        parametrization of the surface between wetting and non-Wetting phase.
     *
     *        The position is determined based on the coordinate of
     *        the vertex belonging to the considered sub control volume.
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub-control volume */
    const AwnSurfaceParams & aWettingNonWettingSurfaceParams(const Element &element,
                                                             const SubControlVolume &scv,
                                                             const ElementSolutionVector &elemSol) const
    {
        const auto& globalPos =  scv.dofPosition();
        if (inFF_(globalPos) )
            return aWettingNonWettingSurfaceParamsFreeFlow_  ;
        else if (inPM_(globalPos))
            return aWettingNonWettingSurfaceParams_ ;
        else             DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);
     }

    /*!\brief Return a reference to the container object for the
     *        parametrization of the surface between non-Wetting and solid phase.
     *
     *        The position is determined based on the coordinate of
     *        the vertex belonging to the considered sub control volume.
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub-control volume */
    const AnsSurfaceParams & aNonWettingSolidSurfaceParams(const Element &element,
                                                             const SubControlVolume &scv,
                                                             const ElementSolutionVector &elemSol) const
    {
        const auto& globalPos =  scv.dofPosition();
        if (inFF_(globalPos) )
            return aNonWettingSolidSurfaceParamsFreeFlow_  ;
        else if (inPM_(globalPos))
            return aNonWettingSolidSurfaceParams_ ;
        else             DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);
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
    const Scalar pcMax(const Element &element,
                       const SubControlVolume &scv,
                       const ElementSolutionVector &elemSol) const
    { return aWettingNonWettingSurfaceParams_.pcMax() ; }


    /*!\brief Return the characteristic length for the mass transfer.
     *
     *        The position is determined based on the coordinate of
     *        the vertex belonging to the considered sub controle volume.
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub control volume */
    const Scalar characteristicLength(const Element & element,
                                      const SubControlVolume &scv,
                                      const ElementSolutionVector &elemSol) const

    {
        const auto& globalPos =  scv.dofPosition();
        return characteristicLengthAtPos(globalPos);
    }

    /*!\brief Return the characteristic length for the mass transfer.
     * \param globalPos The position in global coordinates.*/
    const Scalar characteristicLengthAtPos(const  GlobalPosition & globalPos) const
    {
        if (inFF_(globalPos) )
            return characteristicLengthFF_ ;
        else if (inPM_(globalPos))
            return characteristicLengthPM_ ;
        else             DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);
    }

    /*!\brief Return the pre factor the the energy transfer
     *
     *        The position is determined based on the coordinate of
     *        the vertex belonging to the considered sub controle volume.
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub control volume */
    const Scalar factorEnergyTransfer(const Element &element,
                                      const SubControlVolume &scv,
                                      const ElementSolutionVector &elemSol) const
    {
       const auto& globalPos =  scv.dofPosition();
       return factorEnergyTransferAtPos(globalPos);
    }

    /*!\brief Return the pre factor the the energy transfer
     * \param globalPos The position in global coordinates.*/
    const Scalar factorEnergyTransferAtPos(const  GlobalPosition & globalPos) const
    {
        if (inFF_(globalPos) )
            return factorEnergyTransfer_ ;
        else if (inPM_(globalPos))
            return factorEnergyTransfer_ ;
        else             DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);
    }

    /*!\brief Return the pre factor the the mass transfer
     *
     *        The position is determined based on the coordinate of
     *        the vertex belonging to the considered sub controle volume.
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub control volume */
    const Scalar factorMassTransfer(const Element &element,
                                      const SubControlVolume &scv,
                                      const ElementSolutionVector &elemSol) const
    {
       const auto& globalPos =  scv.dofPosition();
        return factorMassTransferAtPos(globalPos);

    }

    /*!\brief Return the pre factor the the mass transfer
     * \param globalPos The position in global coordinates.*/
    const Scalar factorMassTransferAtPos(const  GlobalPosition & globalPos) const
    {
        if (inFF_(globalPos) )
            return factorMassTransfer_ ;
        else if (inPM_(globalPos))
            return factorMassTransfer_ ;
        else             DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);
    }


    /*!
     * \brief Returns the heat capacity \f$[J / (kg K)]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The global position
     */
    Scalar solidHeatCapacityAtPos(const GlobalPosition& globalPos) const
    {  return solidHeatCapacity_ ;}  // specific heat capacity of solid  [J / (kg K)]

    /*!
     * \brief Returns the mass density \f$[kg / m^3]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The global position
     */
    Scalar solidDensityAtPos(const GlobalPosition& globalPos) const
    {return solidDensity_ ;} // density of solid [kg/m^3]

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the porous material.
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The global position
     */
    Scalar solidThermalConductivityAtPos(const GlobalPosition& globalPos) const
    { return solidThermalConductivity_ ;} // conductivity of solid  [W / (m K ) ]

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

    AwnSurfaceParams    aWettingNonWettingSurfaceParams_;
    AnsSurfaceParams    aNonWettingSolidSurfaceParams_ ;

    AwnSurfaceParams    aWettingNonWettingSurfaceParamsFreeFlow_;
    AnsSurfaceParams    aNonWettingSolidSurfaceParamsFreeFlow_ ;

    Scalar pcMax_ ;

    // Porous Medium Domain
    Scalar intrinsicPermeabilityPM_ ;
    Scalar porosityPM_ ;
    Scalar heightPM_ ;
    Scalar factorEnergyTransfer_ ;
    Scalar factorMassTransfer_ ;
    Scalar characteristicLengthPM_ ;
    MaterialLawParams   materialParamsPM_ ;

    // Free Flow Domain
    Scalar porosityFF_ ;
    Scalar intrinsicPermeabilityFF_ ;
    Scalar characteristicLengthFF_ ;
    MaterialLawParams   materialParamsFF_ ;

    // solid parameters
    Scalar solidDensity_ ;
    Scalar solidThermalConductivity_ ;
    Scalar solidHeatCapacity_ ;

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
