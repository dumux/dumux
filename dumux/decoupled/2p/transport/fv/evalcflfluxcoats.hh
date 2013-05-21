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
#ifndef DUMUX_EVALCFLFLUX_COATS_HH
#define DUMUX_EVALCFLFLUX_COATS_HH

/**
 * @file
 * @brief  CFL-flux-function to evaluate a CFL-Condition after Coats 2003
 */

#include <dumux/decoupled/common/impetproperties.hh> 
#include "evalcflflux.hh"

namespace Dumux
{
/*!\ingroup IMPES
 * @brief  CFL-flux-function to evaluate a CFL-Condition after Coats 2003
 *
 * tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class EvalCflFluxCoats: public EvalCflFlux<TypeTag>
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename SpatialParams::MaterialLaw MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;

    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;

    enum
        {
            dim = GridView::dimension, dimWorld = GridView::dimensionworld
        };
    enum
        {
            wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
            eqIdxPress = Indices::pressureEqIdx,
            eqIdxSat = Indices::satEqIdx,
            numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
        };

    enum
        {
            pw = Indices::pressureW,
            pn = Indices::pressureNW,
            pglobal = Indices::pressureGlobal,
            vw = Indices::velocityW,
            vn = Indices::velocityNW,
            vt = Indices::velocityTotal,
            Sw = Indices::saturationW,
            Sn = Indices::saturationNW
        };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;

public:
    //! \brief Initializes the cfl-flux-model
    void initialize()
    {
        ElementIterator element = problem_.gridView().template begin<0>();
        FluidState fluidState;
        fluidState.setPressure(wPhaseIdx, problem_.referencePressure(*element));
        fluidState.setPressure(nPhaseIdx, problem_.referencePressure(*element));
        fluidState.setTemperature(problem_.temperature(*element));
        fluidState.setSaturation(wPhaseIdx, 1.);
        fluidState.setSaturation(nPhaseIdx, 0.);
        density_[wPhaseIdx] = FluidSystem::density(fluidState, wPhaseIdx);
        density_[nPhaseIdx] = FluidSystem::density(fluidState, nPhaseIdx);
    }

    /*! \brief adds a flux to the cfl-criterion evaluation
     *
     * \copydetails EvalCflFlux::addFlux(Scalar&,Scalar&,Scalar&,Scalar&,Scalar,const Element&,int)
     */
    void addFlux(Scalar& lambdaW, Scalar& lambdaNW, Scalar& viscosityW, Scalar& viscosityNW, Scalar flux,
                 const Element& element, int phaseIdx = -1)
    {
        addDefaultFlux(flux, phaseIdx);
    }

    /*! \brief adds a flux to the cfl-criterion evaluation
     *
     * \copydetails EvalCflFlux::addFlux(Scalar&,Scalar&,Scalar&,Scalar&,Scalar,const Intersection&,int)
     */
    void addFlux(Scalar& lambdaW, Scalar& lambdaNW, Scalar& viscosityW, Scalar& viscosityNW, Scalar flux,
                 const Intersection& intersection, int phaseIdx = -1)
    {
        addDefaultFlux(flux, phaseIdx);
        addCoatsFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, flux, intersection, phaseIdx);
    }

    /*! \brief Returns the CFL flux-function
     *
     * \copydetails EvalCflFlux::getCFLFluxFunction(const Element&)
     */
    Scalar getCFLFluxFunction(const Element& element)
    {
    	if (rejectForTimeStepping_)
    		return 1e100;

        Scalar cflFluxDefault = getCFLFluxFunctionDefault();

        if (std::isnan(cflFluxFunctionCoats_) || std::isinf(cflFluxFunctionCoats_))
        {
            return 0.99 / cflFluxDefault;
        }
        else if (cflFluxFunctionCoats_ <= 0)
        {
            return 0.99 / cflFluxDefault;
        }
        else if (cflFluxDefault > cflFluxFunctionCoats_)
        {
        	return 0.99 / cflFluxDefault;
        }
        else
        {
            return 1.0 / cflFluxFunctionCoats_;
        }
    }

    /*! \brief  Returns the CFL time-step
     *
     * \copydetails EvalCflFlux::getDt(const Element&)
     */
    Scalar getDt(const Element& element)
    {
    	if (rejectForTimeStepping_)
    		return 1e100;

        Scalar porosity = std::max(problem_.spatialParams().porosity(element), porosityThreshold_);

        return (getCFLFluxFunction(element) * porosity * element.geometry().volume());
    }

    //! \brief  Resets the Timestep-estimator
    void reset()
    {
        cflFluxFunctionCoats_ = 0;
        rejectForTimeStepping_ = false;
        fluxWettingOut_ = 0;
        fluxNonwettingOut_ = 0;
        fluxIn_ = 0;
        fluxOut_ = 0;
    }

    /*! \brief Constructs an EvalCflFluxDefault object
     *
     * \param problem A problem type object
     */
    EvalCflFluxCoats(Problem& problem) :
        problem_(problem), epsDerivative_(5e-3), threshold_(1e-12)
    {
        reset();
        density_[wPhaseIdx] = 0.;
        density_[nPhaseIdx] = 0.;

        porosityThreshold_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Impet, PorosityThreshold);
    }

private:
    Scalar getCFLFluxFunctionDefault()
    {
        if (std::isnan(fluxIn_) || std::isinf(fluxIn_))
        {
            fluxIn_ = 1e-100;
        }

        Scalar cFLFluxIn = fluxIn_;


        Scalar cFLFluxOut = 0;

        if (velocityType_ == vt)
        {
            if (std::isnan(fluxOut_) || std::isinf(fluxOut_))
            {
                fluxOut_ = 1e-100;
            }

            cFLFluxOut = fluxOut_;
        }
        else
        {
            if (std::isnan(fluxWettingOut_) || std::isinf(fluxWettingOut_))
            {
                fluxWettingOut_ = 1e-100;
            }
            if (std::isnan(fluxNonwettingOut_) || std::isinf(fluxNonwettingOut_))
            {
                fluxNonwettingOut_ = 1e-100;
            }

            cFLFluxOut = std::max(fluxWettingOut_, fluxNonwettingOut_);
        }


        //determine timestep
        Scalar cFLFluxFunction = std::max(cFLFluxIn, cFLFluxOut);

        return cFLFluxFunction;
    }

    void addDefaultFlux(Scalar flux,int phaseIdx);

    void addCoatsFlux(Scalar& lambdaW, Scalar& lambdaNW, Scalar& viscosityW, Scalar& viscosityNW, Scalar flux,
                      const Intersection& intersection, int phaseIdx);

    Problem& problem_;//problem data
    Scalar cflFluxFunctionCoats_;
    Scalar fluxWettingOut_;
    Scalar fluxNonwettingOut_;
    Scalar fluxOut_;
    Scalar fluxIn_;
    bool rejectForTimeStepping_;
    Scalar density_[numPhases];
    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PressureFormulation);
    static const int velocityType_ = GET_PROP_VALUE(TypeTag, VelocityFormulation);
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, SaturationFormulation);
    const Scalar epsDerivative_;
    const Scalar threshold_;
    Scalar porosityThreshold_;
};

template<class TypeTag>
void EvalCflFluxCoats<TypeTag>::addDefaultFlux(Scalar flux, int phaseIdx)
{
    switch (phaseIdx)
    {
    case wPhaseIdx:
        {
            //for time step criterion
            if (flux >= 0)
            {
                fluxWettingOut_ += flux;
            }
            if (flux < 0)
            {
                fluxIn_ -= flux;
            }

            break;
        }

        //for time step criterion if the non-wetting phase velocity is used
    case nPhaseIdx:
        {
            if (flux >= 0)
            {
                fluxNonwettingOut_ += flux;
            }
            if (flux < 0)
            {
                fluxIn_ -= flux;
            }

            break;
        }
    default:
        {
            if (flux >= 0)
            {
                fluxOut_ += flux;
            }
            if (flux < 0)
            {
                fluxIn_ -= flux;
            }

            break;
        }
    }
}

/*! \brief adds a flux to the cfl-criterion evaluation
 *
 * \copydetails EvalCflFlux::addFlux(Scalar&,Scalar&,Scalar&,Scalar&,Scalar,const Intersection&,int)
 */
template<class TypeTag>
void EvalCflFluxCoats<TypeTag>::addCoatsFlux(Scalar& lambdaW, Scalar& lambdaNW, Scalar& viscosityW, Scalar& viscosityNW, Scalar flux,
                                             const Intersection& intersection, int phaseIdx)
{
	if (rejectForTimeStepping_)
		return;

    Scalar lambdaT = (lambdaW + lambdaNW);

    ElementPointer element = intersection.inside();

    //coordinates of cell center
    const GlobalPosition& globalPos = element->geometry().center();

    // cell index
    int globalIdxI = problem_.variables().index(*element);

    CellData& cellDataI = problem_.variables().cellData(globalIdxI);

    if (cellDataI.pressure(wPhaseIdx) < 0.0 || cellDataI.pressure(nPhaseIdx) < 0.0 )
    {
    	rejectForTimeStepping_ = true;
    	cflFluxFunctionCoats_ = 0;
    	return;
    }

    int indexInInside = intersection.indexInInside();

    Scalar satI = cellDataI.saturation(wPhaseIdx);
    Scalar lambdaWI = cellDataI.mobility(wPhaseIdx);
    Scalar lambdaNWI = cellDataI.mobility(nPhaseIdx);

    Scalar dPcdSI = MaterialLaw::dpC_dSw(problem_.spatialParams().materialLawParams(*element), satI);

    const GlobalPosition& unitOuterNormal = intersection.centerUnitOuterNormal();

    if (intersection.neighbor())
    {
        ElementPointer neighbor = intersection.outside();

        const GlobalPosition& globalPosNeighbor = neighbor->geometry().center();

        // distance vector between barycenters
        Dune::FieldVector < Scalar, dimWorld > distVec = globalPosNeighbor - globalPos;

        // compute distance between cell centers
        Scalar dist = distVec.two_norm();

        int globalIdxJ = problem_.variables().index(*neighbor);

        CellData& cellDataJ = problem_.variables().cellData(globalIdxJ);

        if (cellDataJ.pressure(wPhaseIdx) < 0.0 || cellDataJ.pressure(nPhaseIdx) < 0.0 )
        {
        	rejectForTimeStepping_ = true;
        	cflFluxFunctionCoats_ = 0;
        	return;
        }

        //calculate potential gradients
        if (element.level() != neighbor.level())
        {
            rejectForTimeStepping_ = true;
            cflFluxFunctionCoats_ = 0;
            return;
        }

        if (lambdaT <= 0.0)
        {
            return;
        }

        Scalar satJ = cellDataJ.saturation(wPhaseIdx);
        Scalar lambdaWJ = cellDataI.mobility(wPhaseIdx);
        Scalar lambdaNWJ = cellDataI.mobility(nPhaseIdx);

        Scalar dPcdSJ = MaterialLaw::dpC_dSw(problem_.spatialParams().materialLawParams(*neighbor), satJ);

        // compute vectorized permeabilities
        DimMatrix meanPermeability(0);

        problem_.spatialParams().meanK(meanPermeability,
                                       problem_.spatialParams().intrinsicPermeability(*element),
                                       problem_.spatialParams().intrinsicPermeability(*neighbor));

        Dune::FieldVector < Scalar, dim > permeability(0);
        meanPermeability.mv(unitOuterNormal, permeability);

        Scalar transmissibility = (unitOuterNormal * permeability) * intersection.geometry().volume() / dist;

        Scalar satUpw = 0;
        if (cellDataI.fluxData().isUpwindCell(wPhaseIdx, indexInInside))
        {
            satUpw = std::max(satI, 0.0);
        }
        else
        {
            satUpw = std::max(satJ, 0.0);
        }

        Scalar dS = epsDerivative_;

        Scalar satPlus = satUpw + epsDerivative_;
        Scalar satMinus = satUpw;
        if (satMinus - epsDerivative_ > 0.0)
        {
            satMinus -= epsDerivative_;
            dS += epsDerivative_;
        }

        Scalar dLambdaWdS = MaterialLaw::krw(problem_.spatialParams().materialLawParams(*neighbor), std::abs(satPlus)) / viscosityW;
        dLambdaWdS -= MaterialLaw::krw(problem_.spatialParams().materialLawParams(*neighbor), std::abs(satMinus)) / viscosityW;
        dLambdaWdS /= (dS);

        if (cellDataI.fluxData().isUpwindCell(nPhaseIdx, indexInInside))
        {
            satUpw = std::max(1 - satI, 0.0);
        }
        else
        {
            satUpw = std::max(1 - satJ, 0.0);
        }

        dS = epsDerivative_;

        satPlus = satUpw + epsDerivative_;
        satMinus = satUpw;
        if (satMinus - epsDerivative_ > 0.0)
        {
            satMinus -= epsDerivative_;
            dS += epsDerivative_;
        }

        Scalar dLambdaNWdS = MaterialLaw::krn(problem_.spatialParams().materialLawParams(*neighbor), satPlus) / viscosityNW;
        dLambdaNWdS -= MaterialLaw::krn(problem_.spatialParams().materialLawParams(*neighbor), satMinus) / viscosityNW;
        dLambdaNWdS /= (dS);

        Scalar lambdaWCap = 0.5 * (lambdaWI + lambdaWJ);
        Scalar lambdaNWCap = 0.5 * (lambdaNWI + lambdaNWJ);

        if (phaseIdx == wPhaseIdx)
        {
            Scalar pressDiff = cellDataI.pressure(wPhaseIdx) - cellDataJ.pressure(wPhaseIdx);
            cflFluxFunctionCoats_+= transmissibility * lambdaNW * dLambdaWdS * std::abs(pressDiff) / lambdaT;

            cflFluxFunctionCoats_ -= transmissibility * lambdaWCap * lambdaNWCap * (dPcdSI + dPcdSJ) / lambdaT;
        }
        else if (phaseIdx == nPhaseIdx)
        {
            Scalar pressDiff = cellDataI.pressure(nPhaseIdx) - cellDataJ.pressure(nPhaseIdx);
            cflFluxFunctionCoats_ -= transmissibility * lambdaW * dLambdaNWdS * std::abs(pressDiff) / lambdaT;
        }
        else
        {
            if (flux != 0)
            {
                switch (saturationType_)
                {
                case Sw:
                    {
                        cflFluxFunctionCoats_ += dLambdaWdS / (dLambdaWdS + dLambdaNWdS) * std::abs(flux);
                        break;
                    }
                case Sn:
                    {
                        cflFluxFunctionCoats_ +=  dLambdaNWdS / (dLambdaWdS + dLambdaNWdS) * std::abs(flux);
                        break;
                    }
                }
            }
        }
    }
    else
    {
        if (lambdaT <= 0.0)
        {
            return;
        }

        // center of face in global coordinates
        GlobalPosition globalPosFace = intersection.geometry().center();

        //get boundary type
        BoundaryTypes bcType;
        problem_.boundaryTypes(bcType, intersection);

        // distance vector between barycenters
        Dune::FieldVector < Scalar, dimWorld > distVec = globalPosFace - globalPos;

        // compute distance between cell centers
        Scalar dist = distVec.two_norm();

        //permeability vector at boundary
        // compute vectorized permeabilities
        DimMatrix meanPermeability(0);

        problem_.spatialParams().meanK(meanPermeability,
                                       problem_.spatialParams().intrinsicPermeability(*element));

        Dune::FieldVector<Scalar, dim> permeability(0);
        meanPermeability.mv(unitOuterNormal, permeability);

        Scalar transmissibility = (unitOuterNormal * permeability) * intersection.geometry().volume() / dist;

        Scalar satWBound =  cellDataI.saturation(wPhaseIdx);
        if (bcType.isDirichlet(eqIdxSat))
        {
            PrimaryVariables bcValues;
            problem_.dirichlet(bcValues, intersection);
            switch (saturationType_)
            {
            case Sw:
                {
                    satWBound = bcValues[eqIdxSat];
                    break;
                }
            case Sn:
                {
                    satWBound = 1 - bcValues[eqIdxSat];
                    break;
                }
            default:
                {
                    DUNE_THROW(Dune::RangeError, "saturation type not implemented");
                    break;
                }
            }

        }

        Scalar potWBound =  cellDataI.pressure(wPhaseIdx);
        Scalar potNWBound =  cellDataI.pressure(nPhaseIdx);
        Scalar gdeltaZ = (problem_.bboxMax()-globalPosFace) * problem_.gravity();
        if (bcType.isDirichlet(eqIdxPress))
        {
            PrimaryVariables bcValues;
            problem_.dirichlet(bcValues, intersection);
            switch (pressureType_)
            {
            case pw:
                {
                    potWBound = bcValues[eqIdxPress] + density_[wPhaseIdx] * gdeltaZ;
                    potNWBound = bcValues[eqIdxPress] + MaterialLaw::pC(problem_.spatialParams().materialLawParams(*element), satWBound) + density_[nPhaseIdx] * gdeltaZ;
                    break;
                }
            case pn:
                {
                    potWBound = bcValues[eqIdxPress] - MaterialLaw::pC(problem_.spatialParams().materialLawParams(*element), satWBound) + density_[wPhaseIdx] * gdeltaZ;
                    potNWBound = bcValues[eqIdxPress] + density_[nPhaseIdx] * gdeltaZ;
                    break;
                }
            default:
                {
                    DUNE_THROW(Dune::RangeError, "pressure type not implemented");
                    break;
                }
            }
        }
        else if (bcType.isNeumann(eqIdxPress))
        {
            PrimaryVariables bcValues;
            problem_.neumann(bcValues, intersection);

            if (lambdaW != 0 && bcValues[wPhaseIdx] != 0)
            {
                potWBound += bcValues[wPhaseIdx] / (transmissibility * lambdaW);
            }
            if (lambdaNW != 0 && bcValues[nPhaseIdx] != 0)
            {
                potNWBound += bcValues[nPhaseIdx] / (transmissibility * lambdaNW);
            }
        }

        Scalar dPcdSBound = MaterialLaw::dpC_dSw(problem_.spatialParams().materialLawParams(*element), satWBound);

        Scalar lambdaWBound = 0;
        Scalar lambdaNWBound = 0;

        Scalar temperature = problem_.temperature(*element);
        Scalar referencePressure = problem_.referencePressure(*element);
        FluidState fluidState;
        fluidState.setPressure(wPhaseIdx, referencePressure);
        fluidState.setPressure(nPhaseIdx, referencePressure);
        fluidState.setTemperature(temperature);

        Scalar viscosityWBound = FluidSystem::viscosity(fluidState, wPhaseIdx);
        Scalar viscosityNWBound =
            FluidSystem::viscosity(fluidState, nPhaseIdx);
        lambdaWBound = MaterialLaw::krw(problem_.spatialParams().materialLawParams(*element), satWBound) / viscosityWBound;
        lambdaNWBound = MaterialLaw::krn(problem_.spatialParams().materialLawParams(*element), satWBound) / viscosityNWBound;

        Scalar satUpw = 0;
        if (cellDataI.fluxData().isUpwindCell(wPhaseIdx, indexInInside))
        {
            satUpw = std::max(satI, 0.0);
        }
        else
        {
            satUpw = std::max(satWBound, 0.0);
        }

        Scalar dS = epsDerivative_;

        Scalar satPlus = satUpw + epsDerivative_;
        Scalar satMinus = satUpw;
        if (satMinus - epsDerivative_ > 0.0)
        {
            satMinus -= epsDerivative_;
            dS += epsDerivative_;
        }

        Scalar dLambdaWdS = MaterialLaw::krw(problem_.spatialParams().materialLawParams(*element), satPlus) / viscosityW;
        dLambdaWdS -= MaterialLaw::krw(problem_.spatialParams().materialLawParams(*element), satMinus) / viscosityW;
        dLambdaWdS /= (dS);

        if (cellDataI.fluxData().isUpwindCell(nPhaseIdx, indexInInside))
        {
            satUpw = std::max(1 - satI, 0.0);
        }
        else
        {
            satUpw = std::max(1 - satWBound, 0.0);
        }

        dS = epsDerivative_;

        satPlus = satUpw + epsDerivative_;
        satMinus = satUpw;
        if (satMinus - epsDerivative_ > 0.0)
        {
            satMinus -= epsDerivative_;
            dS += epsDerivative_;
        }

        Scalar dLambdaNWdS = MaterialLaw::krn(problem_.spatialParams().materialLawParams(*element), satPlus) / viscosityNW;
        dLambdaNWdS -= MaterialLaw::krn(problem_.spatialParams().materialLawParams(*element), satMinus) / viscosityNW;
        dLambdaNWdS /= (dS);

        Scalar lambdaWCap = 0.5 * (lambdaWI + lambdaWBound);
        Scalar lambdaNWCap = 0.5 * (lambdaNWI + lambdaNWBound);

        if (phaseIdx == wPhaseIdx)
        {
            Scalar potDiff = cellDataI.pressure(wPhaseIdx) - potWBound;
            cflFluxFunctionCoats_ += transmissibility * lambdaNW * dLambdaWdS * std::abs(potDiff) / lambdaT;

            cflFluxFunctionCoats_ -= transmissibility * lambdaWCap * lambdaNWCap * (dPcdSI + dPcdSBound) / lambdaT;
        }
        else if (phaseIdx == nPhaseIdx)
        {
            Scalar potDiff = cellDataI.pressure(nPhaseIdx) - potNWBound;
            cflFluxFunctionCoats_ -= transmissibility * lambdaW * dLambdaNWdS * std::abs(potDiff) / lambdaT;
        }
        else
        {
            if (flux != 0)
            {
                switch (saturationType_)
                {
                case Sw:
                    {
                        cflFluxFunctionCoats_ += dLambdaWdS / (dLambdaWdS + dLambdaNWdS) * std::abs(flux);
                        break;
                    }
                case Sn:
                    {
                        cflFluxFunctionCoats_ +=  dLambdaNWdS / (dLambdaWdS + dLambdaNWdS) * std::abs(flux);
                        break;
                    }
                }
            }
        }
    }
}

}

#endif
