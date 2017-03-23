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
 * \brief This class contains the volume variables required for the
 *        modules which require the specific interfacial area between
 *        fluid phases.
 *
 * This files contains all specializations which use 'real'
 * interfacial areas.
 */
#ifndef DUMUX_MPNC_VOLUME_VARIABLES_IA_KINETIC_HH
#define DUMUX_MPNC_VOLUME_VARIABLES_IA_KINETIC_HH

#include <dumux/common/dimensionlessnumbers.hh>
#include "volumevariablesia.hh"
#include "propertieskinetic.hh"

namespace Dumux
{


////////////////////////////////////////////////////////////////////////////////////////////////////
// specialization for the case of kinetic mass AND energy transfer
////////////////////////////////////////////////////////////////////////////////////////////////////
/*!
 * \brief Volume variables required for the modules which require the specific interfacial
 *        area between fluid phases. This is the specialization for the case of kinetic
 *         mass AND energy transfer
 */
template <class TypeTag, bool enableKinetic >
class MPNCVolumeVariablesIA<TypeTag, enableKinetic, /*numEnergyEqs=*/3>
{
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename FluidSystem::ParameterCache ParameterCache;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;


    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { wPhaseIdx = FluidSystem::wPhaseIdx };
    enum { nPhaseIdx = FluidSystem::nPhaseIdx };
    enum { sPhaseIdx = FluidSystem::sPhaseIdx };
    enum { nCompIdx = FluidSystem::nCompIdx } ;
    enum { wCompIdx = FluidSystem::wCompIdx } ;
    enum { dim = GridView::dimension};
    enum { dimWorld = GridView::dimensionworld};
    enum { numEnergyEqs = Indices::numPrimaryEnergyVars};
    enum { nusseltFormulation = GET_PROP_VALUE(TypeTag, NusseltFormulation)} ;
    enum { sherwoodFormulation = GET_PROP_VALUE(TypeTag, SherwoodFormulation)} ;
    enum { enableDiffusion = GET_PROP_VALUE(TypeTag, EnableDiffusion)} ;


    typedef DimensionlessNumbers<Scalar> DimLessNum ;

    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;


    typedef typename GET_PROP_TYPE(TypeTag, AwnSurface) AwnSurface;
    typedef typename AwnSurface::Params AwnSurfaceParams;

    typedef typename GET_PROP_TYPE(TypeTag, AwsSurface) AwsSurface;
    typedef typename AwsSurface::Params AwsSurfaceParams;

    typedef typename GET_PROP_TYPE(TypeTag, AnsSurface) AnsSurface;
    typedef typename AnsSurface::Params AnsSurfaceParams;


public:
    /*!
     * \brief Updates the volume specific interfacial area [m^2 / m^3] between the phases.
     */
    void update(const VolumeVariables & volVars,
                const FluidState & fluidState,
                const ParameterCache &paramCache,
                const PrimaryVariables &priVars,
                const Problem &problem,
                const Element & element,
                const FVElementGeometry & fvGeometry,
                const unsigned int scvIdx)
    {
        // obtain (standard) material parameters (needed for the residual saturations)
        const MaterialLawParams & materialParams  = problem.spatialParams().materialLawParams(element,fvGeometry,scvIdx);

        //obtain parameters for interfacial area constitutive relations
        const AwnSurfaceParams & aWettingNonWettingSurfaceParams
               = problem.spatialParams().aWettingNonWettingSurfaceParams(element,fvGeometry,scvIdx);
        const AnsSurfaceParams & aNonWettingSolidSurfaceParams
               = problem.spatialParams().aNonWettingSolidSurfaceParams(element,fvGeometry,scvIdx);

        Valgrind::CheckDefined(aWettingNonWettingSurfaceParams);
        Valgrind::CheckDefined(aNonWettingSolidSurfaceParams);

        const Scalar pc = fluidState.pressure(nPhaseIdx) - fluidState.pressure(wPhaseIdx);
        const Scalar Sw = fluidState.saturation(wPhaseIdx);
        Valgrind::CheckDefined(Sw);
        Valgrind::CheckDefined(pc);
        Scalar awn;

#define AwnRegul 0
        // This regularizes the interfacial area between the fluid phases.
        // This makes sure, that
        // a) some saturation cannot be lost: Never leave two phase region.
        // b) We cannot leave the fit region: no crazy (e.g. negative) values possible
//        const Scalar Swr =  aWettingNonWettingSurfaceParams.Swr() ;
//        const Scalar Snr =  aWettingNonWettingSurfaceParams.Snr() ;

        // this just leads to a stalling newton error as soon as this kicks in.
        // May be a spline or sth like this would help, but I do not which derivatives
        // to specify.
#if AwnRegul
        if(Sw < 5e-3 ) // or Sw > (1.-1e-5 )
        {
            awn = 0. ; // 10.; //
        }
        else
#endif
        awn = AwnSurface::interfacialArea(aWettingNonWettingSurfaceParams, materialParams, Sw, pc ); // 10.; //

        interfacialArea_[wPhaseIdx][nPhaseIdx] = awn ; //10. ;//
        interfacialArea_[nPhaseIdx][wPhaseIdx] = interfacialArea_[wPhaseIdx][nPhaseIdx];
        interfacialArea_[wPhaseIdx][wPhaseIdx] = 0. ;

        Scalar ans = AnsSurface::interfacialArea(aNonWettingSolidSurfaceParams, materialParams,Sw, pc ); // 10.; //
//        if (ans <0 )
//            ans = 0 ;

// Switch for using a a_{wn} relations that has some "maximum capillary pressure" as parameter.
// That value is obtained by regularization of the pc(Sw) function.
#if USE_PCMAX
       const Scalar pcMax = problem.spatialParams().pcMax(element,
                                                            fvGeometry,
                                                            scvIdx);
        // I know the solid surface from the pore network. But it is more consistent to use the fit value.
        // solidSurface_   = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.soil.specificSolidsurface);
        solidSurface_   = AnsSurface::interfacialArea(aNonWettingSolidSurfaceParams, materialParams, /*Sw=*/0., pcMax );
        Valgrind::CheckDefined(solidSurface_);

//        if (ans > solidSurface_){
//            const  GlobalPosition & globalPos =  fvGeometry.subContVol[scvIdx].global ;
//            std::stringstream positionString ;
//            positionString << " Here:";
//            for(int i=0; i<dim; i++)
//                positionString << " x"<< (i+1) << "="  << globalPos[i] << " "   ;
//            positionString << "\n";
//
//              std::cout<<"a_{ns} > a_s, set a_{ns}=" << ans <<" to a_{ns}=a_s="<<solidSurface_ << "
//                                    with S_w="<< Sw << " p_c= "<< pc <<  positionString.str() ;
//
//              ans = solidSurface_ ;
//        }

#endif


        interfacialArea_[nPhaseIdx][sPhaseIdx] = ans ; //10. ; //
        interfacialArea_[sPhaseIdx][nPhaseIdx] = interfacialArea_[nPhaseIdx][sPhaseIdx];
        interfacialArea_[nPhaseIdx][nPhaseIdx] = 0. ;

#if USE_PCMAX
        const Scalar aws = solidSurface_ - ans ;
        interfacialArea_[wPhaseIdx][sPhaseIdx] = aws ; //10. ; //
        interfacialArea_[sPhaseIdx][wPhaseIdx] = interfacialArea_[wPhaseIdx][sPhaseIdx];
        interfacialArea_[sPhaseIdx][sPhaseIdx] = 0. ;
#else
        const AwsSurfaceParams & aWettingSolidSurfaceParams
               = problem.spatialParams().aWettingSolidSurfaceParams();
        Valgrind::CheckDefined(aWettingSolidSurfaceParams);
        const Scalar aws = AwsSurface::interfacialArea(aWettingSolidSurfaceParams,materialParams, Sw, pc ); // 10.; //
        interfacialArea_[wPhaseIdx][sPhaseIdx] = aws ;
        interfacialArea_[sPhaseIdx][wPhaseIdx] = interfacialArea_[wPhaseIdx][sPhaseIdx];
        interfacialArea_[sPhaseIdx][sPhaseIdx] = 0. ;
#endif

        Valgrind::CheckDefined(interfacialArea_);

        factorMassTransfer_   = problem.spatialParams().factorMassTransfer(element,
                                                                           fvGeometry,
                                                                           scvIdx);

        factorEnergyTransfer_   = problem.spatialParams().factorEnergyTransfer(element,
                                                                               fvGeometry,
                                                                               scvIdx);

        characteristicLength_   = problem.spatialParams().characteristicLength(element,
                                                                               fvGeometry,
                                                                               scvIdx);

        // setting the dimensionless numbers.
        // obtaining the respective quantities.
        const unsigned int globalVertexIdx = problem.vertexMapper().subIndex(element, scvIdx, dim);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            const Scalar darcyMagVelocity     = problem.model().volumeDarcyMagVelocity(phaseIdx, globalVertexIdx);
            const Scalar dynamicViscosity     = fluidState.viscosity(phaseIdx);
            const Scalar density              = fluidState.density(phaseIdx);
            const Scalar kinematicViscosity   = dynamicViscosity / density;
            const Scalar heatCapacity         = FluidSystem::heatCapacity(fluidState, paramCache, phaseIdx);
            const Scalar thermalConductivity  = FluidSystem::thermalConductivity(fluidState, paramCache, phaseIdx);


            // diffusion coefficient of non-wetting component in wetting phase
            const Scalar diffCoeff = volVars.diffCoeff(phaseIdx, wCompIdx, nCompIdx) ;
            const Scalar porosity = problem.spatialParams().porosity(element,
                                                                   fvGeometry,
                                                                   scvIdx);

            reynoldsNumber_[phaseIdx]   = DimLessNum::reynoldsNumber(darcyMagVelocity,
                                                                     characteristicLength_,
                                                                     kinematicViscosity);

            prandtlNumber_[phaseIdx]    = DimLessNum::prandtlNumber(dynamicViscosity,
                                                                    heatCapacity,
                                                                    thermalConductivity);

            nusseltNumber_[phaseIdx]    = DimLessNum::nusseltNumberForced(reynoldsNumber_[phaseIdx],
                                                                          prandtlNumber_[phaseIdx],
                                                                          porosity,
                                                                          nusseltFormulation);

            schmidtNumber_[phaseIdx]    = DimLessNum::schmidtNumber(dynamicViscosity,
                                                                    density,
                                                                    diffCoeff);

            // If Diffusion is not enabled, Sherwood is divided by zero
            sherwoodNumber_[phaseIdx]   = DimLessNum::sherwoodNumber(reynoldsNumber_[phaseIdx],
                                                                      schmidtNumber_[phaseIdx],
                                                                      sherwoodFormulation);
        }
    }

    /*!
     * \brief The specific interfacial area between two fluid phases [m^2 / m^3]
     *
     * This is _only_ required by the kinetic mass/energy modules
     *
     */
    const Scalar interfacialArea(const unsigned int phaseIIdx, const unsigned int phaseJIdx) const
    {
        // there is no interfacial area between a phase and itself
        assert(phaseIIdx not_eq phaseJIdx);
        return interfacialArea_[phaseIIdx][phaseJIdx];
    }

    //! access function Reynolds Number
    const Scalar reynoldsNumber(const unsigned int phaseIdx) const
    { return reynoldsNumber_[phaseIdx]; }

    //! access function Prandtl Number
    const Scalar prandtlNumber(const unsigned int phaseIdx) const
    { return prandtlNumber_[phaseIdx]; }

    //! access function Nusselt Number
    const Scalar nusseltNumber(const unsigned int phaseIdx) const
    { return nusseltNumber_[phaseIdx]; }

    //! access function Schmidt Number
    const Scalar schmidtNumber(const unsigned int phaseIdx) const
    { return schmidtNumber_[phaseIdx]; }

    //! access function Sherwood Number
    const Scalar sherwoodNumber(const unsigned int phaseIdx) const
    { return sherwoodNumber_[phaseIdx]; }

    //! access function characteristic length
    const Scalar characteristicLength() const
    { return characteristicLength_; }

    //! access function pre factor energy transfer
    const Scalar factorEnergyTransfer() const
    { return factorEnergyTransfer_; }

    //! access function pre factor mass transfer
    const Scalar factorMassTransfer() const
    { return factorMassTransfer_; }

    /*!
     * \brief If running in valgrind this makes sure that all
     *        quantities in the volume variables are defined.
     */
    void checkDefined() const
    {
#if !defined NDEBUG && HAVE_VALGRIND
        Valgrind::CheckDefined(reynoldsNumber_);
        Valgrind::CheckDefined(prandtlNumber_);
        Valgrind::CheckDefined(nusseltNumber_);
        Valgrind::CheckDefined(schmidtNumber_);
        Valgrind::CheckDefined(sherwoodNumber_);
        Valgrind::CheckDefined(characteristicLength_);
        Valgrind::CheckDefined(factorEnergyTransfer_);
        Valgrind::CheckDefined(factorMassTransfer_);
        Valgrind::CheckDefined(interfacialArea_);
#endif
    }

private:
    //! dimensionless numbers
    Scalar reynoldsNumber_[numPhases];
    Scalar prandtlNumber_[numPhases];
    Scalar nusseltNumber_[numPhases];
    Scalar schmidtNumber_[numPhases];
    Scalar sherwoodNumber_[numPhases];
    Scalar characteristicLength_;
    Scalar factorEnergyTransfer_;
    Scalar factorMassTransfer_;
    Scalar solidSurface_ ;
    Scalar interfacialArea_[numPhases+1][numPhases+1];
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// specialization for the case of NO kinetic mass but kinteic energy transfer of a fluid mixture and solid
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*!
 * \brief Volume variables required for the modules which require the specific interfacial
 *        area between fluid phases. This is the specialization for the case of NO kinetic
 *        mass but kinteic energy transfer of a fluid mixture and solid
 */
template <class TypeTag>
class MPNCVolumeVariablesIA<TypeTag, /*enableKinetic=*/false, /*numEnergyEqs=*/2>
{
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename FluidSystem::ParameterCache ParameterCache;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { dim = GridView::dimension};
    enum { dimWorld = GridView::dimensionworld};
    enum { numEnergyEqs = Indices::numPrimaryEnergyVars};
    enum { nusseltFormulation = GET_PROP_VALUE(TypeTag, NusseltFormulation)} ;

    typedef DimensionlessNumbers<Scalar> DimLessNum ;
public:
    /*!
     * \brief Updates the volume specific interfacial area [m^2 / m^3] between the phases.
     */
    void update(const VolumeVariables & volVars,
                const FluidState & fluidState,
                const ParameterCache &paramCache,
                const PrimaryVariables &priVars,
                const Problem &problem,
                const Element & element,
                const FVElementGeometry & fvGeometry,
                const unsigned int scvIdx)
    {
        factorMassTransfer_   = problem.spatialParams().factorMassTransfer(element,
                                                                           fvGeometry,
                                                                           scvIdx);

        factorEnergyTransfer_   = problem.spatialParams().factorEnergyTransfer(element,
                                                                               fvGeometry,
                                                                               scvIdx);

        characteristicLength_   = problem.spatialParams().characteristicLength(element,
                                                                               fvGeometry,
                                                                               scvIdx);

        // setting the dimensionless numbers.
        // obtaining the respective quantities.
        const unsigned int globalVertexIdx = problem.vertexMapper().subIndex(element, scvIdx, dim);

        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            const Scalar darcyMagVelocity     = problem.model().volumeDarcyMagVelocity(phaseIdx, globalVertexIdx);
            const Scalar dynamicViscosity     = fluidState.viscosity(phaseIdx);
            const Scalar density              = fluidState.density(phaseIdx);
            const Scalar kinematicViscosity   = dynamicViscosity / density;
            const Scalar heatCapacity         = FluidSystem::heatCapacity(fluidState,
                                                                          paramCache,
                                                                          phaseIdx);
            const Scalar thermalConductivity  = FluidSystem::thermalConductivity(fluidState,
                                                                           paramCache,
                                                                           phaseIdx);

            const Scalar porosity = problem.spatialParams().porosity(element,
                                                                   fvGeometry,
                                                                   scvIdx);

            reynoldsNumber_[phaseIdx]   = DimLessNum::reynoldsNumber(darcyMagVelocity,
                                                                     characteristicLength_,
                                                                     kinematicViscosity);

            prandtlNumber_[phaseIdx]    = DimLessNum::prandtlNumber(dynamicViscosity,
                                                                    heatCapacity,
                                                                    thermalConductivity);


            nusseltNumber_[phaseIdx]    = DimLessNum::nusseltNumberForced(reynoldsNumber_[phaseIdx],
                                                                          prandtlNumber_[phaseIdx],
                                                                          porosity,
                                                                          nusseltFormulation);
        }
    }


    //! access function Reynolds Number
    const Scalar reynoldsNumber(const unsigned int phaseIdx) const
    { return reynoldsNumber_[phaseIdx]; }

    //! access function Prandtl Number
    const Scalar prandtlNumber(const unsigned int phaseIdx) const
    { return prandtlNumber_[phaseIdx]; }

    //! access function Nusselt Number
    const Scalar nusseltNumber(const unsigned int phaseIdx) const
    { return nusseltNumber_[phaseIdx]; }

    //! access function characteristic length
    const Scalar characteristicLength() const
    { return characteristicLength_; }

    //! access function pre factor energy transfer
    const Scalar factorEnergyTransfer() const
    { return factorEnergyTransfer_; }

    //! access function pre factor mass transfer
    const Scalar factorMassTransfer() const
    { return factorMassTransfer_; }

    /*!
     * \brief If running in valgrind this makes sure that all
     *        quantities in the volume variables are defined.
     */
    void checkDefined() const
    {
#if !defined NDEBUG && HAVE_VALGRIND
        Valgrind::CheckDefined(reynoldsNumber_);
        Valgrind::CheckDefined(prandtlNumber_);
        Valgrind::CheckDefined(nusseltNumber_);
        Valgrind::CheckDefined(characteristicLength_);
        Valgrind::CheckDefined(factorEnergyTransfer_);
        Valgrind::CheckDefined(factorMassTransfer_);
#endif
    }

private:
    //! dimensionless numbers
    Scalar reynoldsNumber_[numPhases];
    Scalar prandtlNumber_[numPhases];
    Scalar nusseltNumber_[numPhases];
    Scalar characteristicLength_;
    Scalar factorEnergyTransfer_;
    Scalar factorMassTransfer_;
    Scalar solidSurface_ ;
};


////////////////////////////////////////////////////////////////////////////////////////////////////
// specialization for the case of (only) kinetic mass transfer
////////////////////////////////////////////////////////////////////////////////////////////////////
/*!
 * \brief Volume variables required for the modules which require the specific interfacial
 *        area between fluid phases. This is the specialization specialization for the case
 *        of (only) kinetic mass transfer
 */
template <class TypeTag>
class MPNCVolumeVariablesIA<TypeTag, /*enableKinetic=*/true, /*numEnergyEqs=*/0>
{
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename FluidSystem::ParameterCache ParameterCache;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef DimensionlessNumbers<Scalar> DimLessNum;

    typedef typename GET_PROP_TYPE(TypeTag, AwnSurface) AwnSurface;
    typedef typename AwnSurface::Params AwnSurfaceParams;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { wPhaseIdx = FluidSystem::wPhaseIdx };
    enum { nPhaseIdx = FluidSystem::nPhaseIdx };
    enum { wCompIdx  = FluidSystem::wCompIdx };
    enum { nCompIdx  = FluidSystem::nCompIdx };
    enum { dim       = GridView::dimension};
    enum { dimWorld  = GridView::dimensionworld};
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { sherwoodFormulation = GET_PROP_VALUE(TypeTag, SherwoodFormulation)} ;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;


public:
    /*!
     * \brief Updates the volume specific interfacial area [m^2 / m^3] between the phases.
     */
    void update(const VolumeVariables & volVars,
                const FluidState & fluidState,
                const ParameterCache & paramCache,
                const PrimaryVariables &priVars,
                const Problem &problem,
                const Element & element,
                const FVElementGeometry & fvGeometry,
                const unsigned int scvIdx)
    {
        //obtain parameters for awnsurface
        const AwnSurfaceParams & awnSurfaceParams = problem.spatialParams().aWettingNonWettingSurfaceParams(element,fvGeometry,scvIdx) ;

        // obtain (standard) material parameters (needed for the residual saturations)
        const MaterialLawParams & materialParams  = problem.spatialParams().materialLawParams(element,fvGeometry,scvIdx) ;

        Valgrind::CheckDefined(awnSurfaceParams);
        const Scalar Sw = fluidState.saturation(wPhaseIdx) ;
        const Scalar pc = fluidState.pressure(nPhaseIdx) - fluidState.pressure(wPhaseIdx);

        // so far there is only a model for kinetic mass transfer between fluid phases
        interfacialArea_ = AwnSurface::interfacialArea(awnSurfaceParams, materialParams, Sw, pc );

        Valgrind::CheckDefined(interfacialArea_);

        factorMassTransfer_   = problem.spatialParams().factorMassTransfer(element,
                                                                           fvGeometry,
                                                                           scvIdx);

        characteristicLength_   = problem.spatialParams().characteristicLength(element,
                                                                               fvGeometry,
                                                                               scvIdx);
        // setting the dimensionless numbers.
        // obtaining the respective quantities.
        int globalVertexIdx = problem.vertexMapper().subIndex(element, scvIdx, dim);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            const Scalar darcyMagVelocity     = problem.model().volumeDarcyMagVelocity(phaseIdx, globalVertexIdx);
            const Scalar dynamicViscosity     = fluidState.viscosity(phaseIdx);
            const Scalar density              = fluidState.density(phaseIdx);
            const Scalar kinematicViscosity   = dynamicViscosity / density;

            // diffusion coefficient of non-wetting component in wetting phase
            const Scalar diffCoeff = volVars.diffCoeff(phaseIdx, wCompIdx, nCompIdx) ;

            reynoldsNumber_[phaseIdx]   = DimLessNum::reynoldsNumber(darcyMagVelocity,
                                                                     characteristicLength_,
                                                                     kinematicViscosity);

            schmidtNumber_[phaseIdx]    = DimLessNum::schmidtNumber(dynamicViscosity,
                                                                    density,
                                                                    diffCoeff);

            sherwoodNumber_[phaseIdx]   = DimLessNum::sherwoodNumber(reynoldsNumber_[phaseIdx],
                                                                      schmidtNumber_[phaseIdx],
                                                                      sherwoodFormulation);
        }
    }

    /*!
     * \brief The specific interfacial area between two fluid phases [m^2 / m^3]
     */
    const Scalar interfacialArea(const unsigned int phaseIIdx, const unsigned int phaseJIdx) const
    {
        // so far there is only a model for kinetic mass transfer between fluid phases
        assert((phaseIIdx == nPhaseIdx && phaseJIdx == wPhaseIdx)
              || (phaseIIdx == wPhaseIdx && phaseJIdx == nPhaseIdx));
        return interfacialArea_;
    }

    //! access function Reynolds Number
    const Scalar reynoldsNumber(const unsigned int phaseIdx) const
    { return reynoldsNumber_[phaseIdx]; }

    //! access function Schmidt Number
    const Scalar schmidtNumber(const unsigned int phaseIdx) const
    { return schmidtNumber_[phaseIdx]; }

    //! access function Sherwood Number
    const Scalar sherwoodNumber(const unsigned int phaseIdx) const
    { return sherwoodNumber_[phaseIdx]; }

    //! access function characteristic length
    const Scalar characteristicLength() const
    { return characteristicLength_; }

    //! access function pre factor mass transfer
    const Scalar factorMassTransfer() const
    { return factorMassTransfer_; }

    /*!
     * \brief If running in valgrind this makes sure that all
     *        quantities in the volume variables are defined.
     */
    void checkDefined() const
    {
        Valgrind::CheckDefined(interfacialArea_);
        Valgrind::CheckDefined(characteristicLength_);
        Valgrind::CheckDefined(factorMassTransfer_);
        Valgrind::CheckDefined(reynoldsNumber_);
        Valgrind::CheckDefined(schmidtNumber_);
        Valgrind::CheckDefined(sherwoodNumber_);
    }

private:
    Scalar characteristicLength_;
    Scalar factorMassTransfer_;
    Scalar solidSurface_ ;
    Scalar interfacialArea_ ;
    Scalar sherwoodNumber_[numPhases] ;
    Scalar schmidtNumber_[numPhases] ;
    Scalar reynoldsNumber_[numPhases] ;
};



} // namespace Dumux

#endif
