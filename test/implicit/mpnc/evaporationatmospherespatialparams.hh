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
 * \file evaporationatmospherespatialparams.hh
 * \brief spatialparameters for the kinetic test-case of the mpnc model. "Poor-mans" coupling of free-flow and porous medium.
 *
 */
#ifndef DUMUX_EVAPORATION_ATMOSPHERE_SPATIALPARAMS_HH
#define DUMUX_EVAPORATION_ATMOSPHERE_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/implicitspatialparams.hh>
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

#include <dumux/implicit/mpnc/mpncmodelkinetic.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class EvaporationAtmosphereSpatialParams;

namespace Properties
{
// The spatial params TypeTag
NEW_TYPE_TAG(EvaporationAtmosphereSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(EvaporationAtmosphereSpatialParams, SpatialParams, Dumux::EvaporationAtmosphereSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(EvaporationAtmosphereSpatialParams, MaterialLaw)
{
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum {wPhaseIdx   = FluidSystem::wPhaseIdx};

private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    // select appropriate parametrization
#define linear      0 // pc(S_w)
#define RegVG       0 // pc(S_w)
#define RegBC       1 // pc(S_w)
#define RegVGofT    0 // pc(S_w, T) // adapt also the setting of the fluidstate
#define VG          0 // pc(S_w)

    // define the material law
    #if linear
        typedef LinearMaterial<Scalar>          EffectiveLaw;
    #endif

    #if RegVG
        typedef RegularizedVanGenuchten<Scalar> EffectiveLaw;
    #endif

    #if RegBC
        typedef RegularizedBrooksCorey<Scalar> EffectiveLaw;
    #endif

    #if RegVGofT
        typedef RegularizedVanGenuchtenOfTemperature<Scalar> EffectiveLaw;
    #endif

    #if VG
        typedef VanGenuchten<Scalar> EffectiveLaw;
    #endif

    typedef EffToAbsLaw<EffectiveLaw>       TwoPMaterialLaw;
    public:
//        typedef TwoPOfTAdapter<wPhaseIdx, TwoPMaterialLaw> type; // adapter for incorporating temperature effects on pc-S
        typedef TwoPAdapter<wPhaseIdx, TwoPMaterialLaw> type;
};



// Set the interfacial area relation: wetting -- non-wetting
SET_PROP(EvaporationAtmosphereSpatialParams, AwnSurface)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef AwnSurfacePcMaxFct<Scalar>     EffectiveIALaw;
    //    typedef AwnSurfacePolynomial2ndOrder<Scalar>      EffectiveIALaw;
public:
    typedef EffToAbsLawIA<EffectiveIALaw, MaterialLawParams> type;
};


// Set the interfacial area relation: wetting -- solid
SET_PROP(EvaporationAtmosphereSpatialParams, AwsSurface)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef AwnSurfacePolynomial2ndOrder<Scalar>  EffectiveIALaw;
public:
    typedef EffToAbsLawIA<EffectiveIALaw, MaterialLawParams> type;
};

// Set the interfacial area relation: non-wetting -- solid
SET_PROP(EvaporationAtmosphereSpatialParams, AnsSurface)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef AwnSurfaceExpSwPcTo3<Scalar>      EffectiveIALaw;
public:
    typedef EffToAbsLawIA<EffectiveIALaw, MaterialLawParams> type;
};

} // end namespace properties

/**
 * \brief Definition of the spatial parameters for the evaporation atmosphere Problem (using a "poor man's coupling")
 */
template<class TypeTag>
class EvaporationAtmosphereSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum {dim=GridView::dimension };
    enum {dimWorld=GridView::dimensionworld};
    enum {wPhaseIdx = FluidSystem::wPhaseIdx};
    enum {nPhaseIdx = FluidSystem::nPhaseIdx};
    enum {sPhaseIdx = FluidSystem::sPhaseIdx};
    enum { numEnergyEquations  = GET_PROP_VALUE(TypeTag, NumEnergyEquations)};
    enum { numPhases       = GET_PROP_VALUE(TypeTag, NumPhases)};
    enum { enableEnergy         = GET_PROP_VALUE(TypeTag, EnableEnergy)};

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar,dim> LocalPosition;

    typedef Dune::FieldVector<Scalar,dim> DimVector;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

public:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    typedef typename GET_PROP_TYPE(TypeTag, AwnSurface) AwnSurface;
    typedef typename AwnSurface::Params AwnSurfaceParams;

    typedef typename GET_PROP_TYPE(TypeTag, AwsSurface) AwsSurface;
    typedef typename AwsSurface::Params AwsSurfaceParams;

    typedef typename GET_PROP_TYPE(TypeTag, AnsSurface) AnsSurface;
    typedef typename AnsSurface::Params AnsSurfaceParams;


    EvaporationAtmosphereSpatialParams(const GridView &gridView)
        : ParentType(gridView)
    {
    }

    ~EvaporationAtmosphereSpatialParams()
    {}

        // load parameters from input file and initialize parameter values
        void setInputInitialize()
        {
            try
            {
                eps_                    = 1e-6;
                heightPM_               = GET_RUNTIME_PARAM(TypeTag, Scalar, Grid.InterfacePos);
                heightDomain_           = GET_RUNTIME_PARAM(TypeTag, Scalar, Grid.YMax);
                // BEWARE! First the input values have to be set, than the material parameters can be set

                // this is the parameter value from file part
                porosityPM_                 = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.PorousMedium.porosity);
                intrinsicPermeabilityPM_    = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.PorousMedium.permeability);

                porosityFF_                 = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.FreeFlow.porosity);
                intrinsicPermeabilityFF_    = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.FreeFlow.permeability);

                solidDensity_               = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.soil.density);
                solidThermalConductivity_    = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.soil.thermalConductivity);
                solidHeatCapacity_               = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.soil.heatCapacity);

                aWettingNonWettingA1_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.soil.aWettingNonWettingA1);
                aWettingNonWettingA2_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.soil.aWettingNonWettingA2);
                aWettingNonWettingA3_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.soil.aWettingNonWettingA3);

                aNonWettingSolidA1_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.soil.aNonWettingSolidA1);
                aNonWettingSolidA2_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.soil.aNonWettingSolidA2);
                aNonWettingSolidA3_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.soil.aNonWettingSolidA3);

                BCPd_           = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.soil.BCPd);
                BClambda_       = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.soil.BClambda);
                Swr_            = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.soil.Swr);
                Snr_            = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.soil.Snr);

                characteristicLengthFF_         = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.FreeFlow.meanPoreSize);

                characteristicLengthPM_         = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.PorousMedium.meanPoreSize);
                factorEnergyTransfer_           = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.PorousMedium.factorEnergyTransfer);
                factorMassTransfer_             = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.PorousMedium.factorMassTransfer);
            }
            catch (Dumux::ParameterException &e) {
                std::cerr << e << ". Abort!\n";
                exit(1) ;
            }
            catch (...) {
                std::cerr << "Unknown exception thrown!\n";
                exit(1);
            }

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
                const MaterialLawParams &materialParams =  materialParamsPM_;
                        Scalar capPress[numPhases];
                //set saturation to inital values, this needs to be done in order for the fluidState to tell me pc
                for (int phaseIdx = 0; phaseIdx < numPhases ; ++phaseIdx) {
                    // set saturation to zero for getting pcmax
                    S[phaseIdx] = 0. ;
                    Scalar TInitial               = GET_RUNTIME_PARAM(TypeTag, Scalar, InitialConditions.TInitial);

                    fluidState.setSaturation(phaseIdx, S[phaseIdx]);
                    if(enableEnergy){
                        if (numEnergyEquations>1)
                            fluidState.setTemperature(phaseIdx,TInitial);
                        else
                            fluidState.setTemperature(TInitial);
                    }
                }

                //obtain pc according to saturation
                MaterialLaw::capillaryPressures(capPress, materialParams, fluidState);
                pcMax_ = std::abs(capPress[0]);

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

            // dummy for aws not necessary: aws = ans - as, ans=0, as=ans(0,pcmax)=0 -> aws=0
        }

    /*! \brief Update the spatial parameters with the flow solution
     *        after a timestep.
     */
    void update(const SolutionVector & globalSolutionFn)
    { }

    /*! \brief Returns the intrinsic permeability
     *
     *        The position is determined based on the coordinate of
     *        the vertex belonging to the considered sub control volume.
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index sub-control volume */
    const Scalar intrinsicPermeability(const Element & element,
                                       const FVElementGeometry & fvGeometry,
                                       const unsigned int scvIdx) const
    {
        const  GlobalPosition & globalPos =  fvGeometry.subContVol[scvIdx].global ;
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
    const Scalar porosity(const Element & element,
                          const FVElementGeometry & fvGeometry,
                          const unsigned int scvIdx) const
    {
        const  GlobalPosition & globalPos =  fvGeometry.subContVol[scvIdx].global ;
//        const  LocalPosition & localPos =  fvGeometry.subContVol[scvIdx].localCenter ;
//        const  GlobalPosition & globalPos =  element.geometry().global(localPos) ;

        if (inFF_(globalPos) )
            return porosityFF_ ;
        else if (inPM_(globalPos))
            return porosityPM_ ;
        else
            DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);
    }

    /*!
     * \brief Return a reference to the material parameters of the material law.
     *
     *        The position is determined based on the coordinate of
     *        the vertex belonging to the considered sub control volume.
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub-control volume */
    const MaterialLawParams & materialLawParams(const Element & element,
                                                const FVElementGeometry & fvGeometry,
                                                const unsigned int scvIdx) const
    {
        const GlobalPosition & globalPos =  fvGeometry.subContVol[scvIdx].global ;
        const MaterialLawParams & materialParams  = materialLawParamsAtPos(globalPos) ;
        return materialParams ;
    }

    /*!
     * \brief Return a reference to the material parameters of the material law.
     *
     * \param globalPos The position in global coordinates.
     */
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
    const AwnSurfaceParams & aWettingNonWettingSurfaceParams(const Element & element,
                                                             const FVElementGeometry & fvGeometry,
                                                             const unsigned int scvIdx) const
    {
        const  GlobalPosition & globalPos =  fvGeometry.subContVol[scvIdx].global ;
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
                                                           const FVElementGeometry &fvGeometry,
                                                           const unsigned int scvIdx) const
    {
        const  GlobalPosition & globalPos =  fvGeometry.subContVol[scvIdx].global ;
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
    const Scalar pcMax(const Element & element,
                       const FVElementGeometry & fvGeometry,
                       const unsigned int scvIdx) const
    { return aWettingNonWettingSurfaceParams_.pcMax() ; }

    /*!\brief Return the characteristic length for the mass transfer.
     *
     *        The position is determined based on the coordinate of
     *        the vertex belonging to the considered sub controle volume.
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub control volume */
    const Scalar characteristicLength(const Element & element,
                                      const FVElementGeometry & fvGeometry,
                                      const unsigned int scvIdx) const
    {   const  GlobalPosition & globalPos =  fvGeometry.subContVol[scvIdx].global ;
        return characteristicLengthAtPos(globalPos); }
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
    const Scalar factorEnergyTransfer(const Element & element,
                                      const FVElementGeometry & fvGeometry,
                                      const unsigned int scvIdx) const
    {   const  GlobalPosition & globalPos =  fvGeometry.subContVol[scvIdx].global ;
        return factorEnergyTransferAtPos(globalPos); }
    /*!\brief Return the pre factor the the energy transfer
     * \param globalPos The position in global coordinates.*/
    const Scalar factorEnergyTransferAtPos(const GlobalPosition & globalPos) const
    {
        if (inFF_(globalPos) )
            return factorEnergyTransfer_ ;
        else if (inPM_(globalPos))
            return factorEnergyTransfer_ ;
        else             DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);
    }

    /*!\brief Return the pre factor for the mass transfer
     *
     *        The position is determined based on the coordinate of
     *        the vertex belonging to the considered sub controle volume.
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub control volume */
    const Scalar factorMassTransfer(const Element & element,
                                    const FVElementGeometry & fvGeometry,
                                    const unsigned int scvIdx) const
    {   const  GlobalPosition & globalPos =  fvGeometry.subContVol[scvIdx].global ;
        return factorMassTransferAtPos(globalPos); }
    /*!\brief Return the pre factor the the mass transfer
     * \param globalPos The position in global coordinates.*/
    const Scalar factorMassTransferAtPos(const GlobalPosition & globalPos) const
    {
        if (inFF_(globalPos) )
            return factorMassTransfer_ ;
        else if (inPM_(globalPos))
            return factorMassTransfer_ ;
        else             DUNE_THROW(Dune::InvalidStateException, "You should not be here: x=" << globalPos[0] << " y= "<< globalPos[dimWorld-1]);
    }


    /*! \brief Returns the heat capacity \f$[J/m^3 K]\f$ of the rock matrix.
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub-control volume where
     *                    the heat capacity needs to be defined */
    const Scalar solidHeatCapacity(const Element & element,
                              const FVElementGeometry & fvGeometry,
                              const unsigned int scvIdx) const
    {  return solidHeatCapacity_ ;}  // specific heat capacity of solid  [J / (kg K)]

    /*!\brief Returns the density \f$[kg/m^3]\f$ of the rock matrix.
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub-control volume */
    const Scalar solidDensity(const Element & element,
                             const FVElementGeometry & fvGeometry,
                             const unsigned int scvIdx) const
    { return solidDensity_ ;} // density of solid [kg/m^3]

    /*!\brief Returns the thermal conductivity \f$[W/(m K)]\f$ of the rock matrix.
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub-control volume */
    const Scalar  solidThermalConductivity(const Element & element,
                                          const FVElementGeometry & fvGeometry,
                                          const unsigned int scvIdx)const
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
    const bool inPM_(const GlobalPosition & globalPos) const
    {       return ( (globalPos[dimWorld-1] > 0. - 1e-6) and (globalPos[dimWorld-1] < (heightPM_ + 1e-6 ) ) );   }

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
    const bool inFF_(const GlobalPosition & globalPos) const
    {       return ( (globalPos[dimWorld-1] < heightDomain_ + 1e-6) and (globalPos[dimWorld-1] > (heightPM_ + 1e-6) ) );   }

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
    const bool inInjection_(const GlobalPosition & globalPos) const
    {       return ( (globalPos[dimWorld-1] < heightDomain_ - 0.25*heightDomain_  + 1e-6) and (globalPos[dimWorld-1] > (heightPM_ + 1e-6) ) );   }

    /*! \brief access function for the depth / height of the porous medium */
    const Scalar heightPM() const
    { return heightPM_; }

private:
    Scalar eps_ ;
    Scalar heightDomain_ ;

    AwnSurfaceParams    aWettingNonWettingSurfaceParams_;
    AnsSurfaceParams 	aNonWettingSolidSurfaceParams_ ;

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
};

}

#endif // GUARDIAN
