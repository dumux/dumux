/*****************************************************************************
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file combustionspatialparams.hh
 *
 * \brief Spatialparameters for the combustionproblem1c. Parameters for the actual simulation domain and an outflow region are provided.
 */
#ifndef DUMUX_COMBUSTION_SPATIALPARAMS_HH
#define DUMUX_COMBUSTION_SPATIALPARAMS_HH

#include <dune/common/parametertreeparser.hh>

#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/heatpipelaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/mp/2padapter.hh>
#include <dumux/material/spatialparams/implicit.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class CombustionSpatialParams;

namespace Properties
{
// Some forward declarations
NEW_PROP_TAG(EnableEnergy);
NEW_PROP_TAG(FluidState);
NEW_PROP_TAG(FluidSystem);
NEW_PROP_TAG(NumEnergyEquations);
NEW_PROP_TAG(NumPhases);

// The spatial params TypeTag
NEW_TYPE_TAG(CombustionSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(CombustionSpatialParams, SpatialParams, CombustionSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(CombustionSpatialParams, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    enum {wPhaseIdx   = FluidSystem::wPhaseIdx};

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

//    actually other people call this Leverett
    typedef HeatPipeLaw<Scalar> EffectiveLaw;

    typedef EffToAbsLaw<EffectiveLaw>       TwoPMaterialLaw;
    public:
        typedef TwoPAdapter<wPhaseIdx, TwoPMaterialLaw> type;
};


}// end namespace properties

/**
 * \brief Definition of the spatial parameters for the one component combustion problem
 *
 */
template<class TypeTag>
class CombustionSpatialParams : public ImplicitSpatialParams<TypeTag>
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
    enum {numEnergyEquations    = GET_PROP_VALUE(TypeTag, NumEnergyEquations)};
    enum {numPhases = GET_PROP_VALUE(TypeTag, NumPhases)};
    enum {enableEnergy          = GET_PROP_VALUE(TypeTag, EnableEnergy)};

    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> DimVector;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

public:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    CombustionSpatialParams(const GridView &gv)
        : ParentType(gv)
    {
    }

    ~CombustionSpatialParams()
    {}

    //! load parameters from input file and initialize parameter values
        void setInputInitialize()
        {
                // BEWARE! First the input values have to be set, then the material parameters can be set

                // this is the parameter value from file part
                porosity_                       = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.PorousMedium.porosity);

                intrinsicPermeabilityOutFlow_   = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.OutFlow.permeabilityOutFlow);
                porosityOutFlow_                = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.OutFlow.porosityOutFlow);
                solidThermalConductivityOutflow_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.OutFlow.soilThermalConductivityOutFlow);

                solidDensity_                       = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.soil.density);
                solidThermalConductivity_           = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.soil.thermalConductivity);
                solidHeatCapacity_                      = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.soil.heatCapacity);

                interfacialTension_                 = GET_RUNTIME_PARAM(TypeTag, Scalar, Constants.interfacialTension);

                Swr_            = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.soil.Swr);
                Snr_            = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.soil.Snr);

                characteristicLength_   = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.PorousMedium.meanPoreSize);
                intrinsicPermeability_  =  (std::pow(characteristicLength_,2.0)  * std::pow(porosity_,3.0)) / (150.0 * std::pow((1.0-porosity_),2.0)); // 1.69e-10 ; //

                factorEnergyTransfer_         = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.PorousMedium.factorEnergyTransfer);
                factorMassTransfer_           = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.PorousMedium.factorMassTransfer);

                lengthPM_               = GET_RUNTIME_PARAM(TypeTag, Scalar,Grid.lengthPM);

            // residual saturations
            materialParams_.setSwr(Swr_) ;
            materialParams_.setSnr(Snr_) ;

            materialParams_.setP0(std::sqrt(porosity_/intrinsicPermeability_));
            materialParams_.setGamma(interfacialTension_); // interfacial tension of water-air at 100°C
        }

    /*! \brief Update the spatial parameters with the flow solution
     *        after a timestep. */
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
        if ( inOutFlow(globalPos) )
            return intrinsicPermeabilityOutFlow_ ;
        else
            return intrinsicPermeability_ ;
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
        if ( inOutFlow(globalPos) )
            return porosityOutFlow_ ;
        else
            return porosity_ ;
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
        const  GlobalPosition & globalPos =  fvGeometry.subContVol[scvIdx].global ;
        const MaterialLawParams & materialParams  = materialLawParamsAtPos(globalPos) ;
        return materialParams ;
    }

    /*!
     * \brief Return a reference to the material parameters of the material law.
     * \param globalPos The position in global coordinates. */
    const MaterialLawParams & materialLawParamsAtPos(const GlobalPosition & globalPos) const
    {
            return materialParams_ ;
    }

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
    { return characteristicLength_ ; }

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
    const Scalar factorEnergyTransferAtPos(const  GlobalPosition & globalPos) const
    {
        return factorEnergyTransfer_ ;
    }


    /*!\brief Return the pre factor the the mass transfer
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
    const Scalar factorMassTransferAtPos(const  GlobalPosition & globalPos) const
    {
        return factorMassTransfer_ ;
    }


    /*! \brief Returns the heat capacity \f$[J/m^3 K]\f$ of the rock matrix.
     * \param element     The finite element
     * \param fvGeometry  The finite volume geometry
     * \param scvIdx      The local index of the sub-control volume where
     *                    the heat capacity needs to be defined */
    const Scalar solidHeatCapacity(const Element &element,
                               const FVElementGeometry &fvGeometry,
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
    const Scalar  solidThermalConductivity(const Element &element,
                                    const FVElementGeometry &fvGeometry,
                                    const unsigned int scvIdx)const
    {
        const  GlobalPosition & globalPos =  fvGeometry.subContVol[scvIdx].global ;
        if ( inOutFlow(globalPos) )
            return solidThermalConductivityOutflow_ ;
        else
            return solidThermalConductivity_ ;
    } // conductivity of solid  [W / (m K ) ]

    /*!
     * \brief Give back whether the tested position (input) is a specific region (right end of porous medium) in the domain
     */
    bool inOutFlow(const GlobalPosition & globalPos) const
    {        return globalPos[0] > (lengthPM_ - eps_) ;    }

    /*!
     * \brief Give back how long the porous medium domain is.
     */
    Scalar lengthPM() const
    {
        return lengthPM_ ;
    }

    /*!
     * \brief Give back the itnerfacial tension
     */
    Scalar interfacialTension() const
    {
        return interfacialTension_ ;
    }

private:
    static constexpr Scalar eps_ = 1e-6;;

    // Porous Medium Domain
    Scalar intrinsicPermeability_ ;
    Scalar porosity_ ;
    Scalar factorEnergyTransfer_ ;
    Scalar factorMassTransfer_ ;
    Scalar characteristicLength_ ;
    MaterialLawParams   materialParams_ ;

    // Outflow Domain
    Scalar intrinsicPermeabilityOutFlow_ ;
    Scalar porosityOutFlow_ ;

    // solid parameters
    Scalar solidDensity_ ;
    Scalar solidThermalConductivity_ ;
    Scalar solidThermalConductivityOutflow_ ;
    Scalar solidHeatCapacity_ ;
    Scalar interfacialTension_ ;


    // capillary pressures parameters
    Scalar Swr_ ;
    Scalar Snr_ ;

    // grid
    Scalar lengthPM_ ;
};

}

#endif // GUARDIAN
