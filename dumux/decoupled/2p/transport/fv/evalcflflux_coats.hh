// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                 *
 *   Institute of Hydraulic Engineering                                      *
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
#ifndef DUMUX_EVALCFLFLUX_COATS_HH
#define DUMUX_EVALCFLFLUX_COATS_HH

/**
 * @file
 * @brief  CFL-flux-function to evaluate a CFL-Condition after Coats 2003
 * @author Markus Wolff
 */

#include "evalcflflux_default.hh"

namespace Dumux
{
/*!\ingroup Saturation2p
 * @brief  CFL-flux-function to evaluate a CFL-Condition after Coats 2003
 *
 *
 * Template parameters are:

 - TypeTag PropertyTag of the problem implementation
 */
template<class TypeTag>
class EvalCflFluxCoats: public EvalCflFluxDefault<TypeTag>
{
private:
    typedef EvalCflFluxDefault<TypeTag> ParentType;typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
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
        Sn = Indices::saturationNW,
    };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;

public:

    void addFlux(Scalar& lambdaW, Scalar& lambdaNW, Scalar& viscosityW, Scalar& viscosityNW, Scalar flux,
            const Element& element, int phaseIdx = -1)
    {
        ParentType::addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, flux, element, phaseIdx);
    }

    void addFlux(Scalar& lambdaW, Scalar& lambdaNW, Scalar& viscosityW, Scalar& viscosityNW, Scalar flux,
            const Intersection& intersection, int phaseIdx = -1)
    {
        Scalar parentMob = 0.5;
        Scalar parentVis = 1;
        //        ParentType::addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, flux, intersection, phaseIdx);
        ParentType::addFlux(parentMob, parentMob, parentVis, parentVis, flux, intersection, phaseIdx);

        //coordinates of cell center
        const GlobalPosition& globalPos = intersection.inside()->geometry().center();

        // cell index
        int globalIdxI = problem_.variables().index(*(intersection.inside()));

        int indexInInside = intersection.indexInInside();

        // get absolute permeability
        const FieldMatrix& permeabilityI(problem_.spatialParameters().intrinsicPermeability(globalPos,
                *(intersection.inside())));

        Scalar satI = problem_.variables().saturation()[globalIdxI];
        Scalar lambdaWI = problem_.variables().mobilityWetting(globalIdxI);
        Scalar lambdaNWI = problem_.variables().mobilityNonwetting(globalIdxI);

        Scalar dPcdSI = MaterialLaw::dpC_dSw(problem_.spatialParameters().materialLawParams(globalPos,
                *(intersection.inside())), satI);

        const Dune::FieldVector<Scalar, dimWorld>& unitOuterNormal = intersection.centerUnitOuterNormal();

        if (intersection.neighbor())
        {
            const GlobalPosition& globalPosNeighbor = intersection.outside()->geometry().center();

            // distance vector between barycenters
            Dune::FieldVector < Scalar, dimWorld > distVec = globalPosNeighbor - globalPos;

            // compute distance between cell centers
            Scalar dist = distVec.two_norm();

            int globalIdxJ = problem_.variables().index(*(intersection.outside()));

            //calculate potential gradients
            Scalar potentialW = 0;
            Scalar potentialNW = 0;

            potentialW = problem_.variables().potentialWetting(globalIdxI, indexInInside);
            potentialNW = problem_.variables().potentialNonwetting(globalIdxI, indexInInside);

            Scalar satJ = problem_.variables().saturation()[globalIdxI];
            Scalar lambdaWJ = problem_.variables().mobilityWetting(globalIdxJ);
            Scalar lambdaNWJ = problem_.variables().mobilityNonwetting(globalIdxJ);

            Scalar dPcdSJ = MaterialLaw::dpC_dSw(problem_.spatialParameters().materialLawParams(globalPosNeighbor,
                    (*intersection.outside())), satJ);

            // get absolute permeability
            const FieldMatrix& permeabilityJ(problem_.spatialParameters().intrinsicPermeability(globalPos,
                    *(intersection.outside())));

            // compute vectorized permeabilities
            FieldMatrix meanPermeability(0);

            //                // harmonic mean of permeability
            for (int x = 0; x < dim; x++)
            {
                meanPermeability[x][x] = 2 * permeabilityI[x][x] * permeabilityJ[x][x] / (permeabilityI[x][x]
                        + permeabilityJ[x][x]);
                for (int y = 0; y < dim; y++)
                {
                    if (x != y)
                    {//use arithmetic mean for the off-diagonal entries to keep the tensor property!
                        meanPermeability[x][y] = 0.5 * (permeabilityI[x][y] + permeabilityJ[x][y]);
                    }
                }
            }

            Dune::FieldVector < Scalar, dim > permeability(0);
            meanPermeability.mv(unitOuterNormal, permeability);

            Scalar transmissibility = (unitOuterNormal * permeability) * intersection.geometry().volume() / dist;

            Scalar satUpw = 0;
            if (potentialW >= 0)
            {
                satUpw = std::max(satI, 0.0);
            }
            else
            {
                satUpw = std::max(satJ, 0.0);
            }

            Scalar dS = eps_;

            Scalar satPlus = satUpw + eps_;
            Scalar satMinus = satUpw;
            if (satMinus - eps_ > 0.0)
            {
                satMinus -= eps_;
                dS += eps_;
            }

            Scalar dLambdaWdS = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(globalPos,
                    (*intersection.outside())), std::abs(satPlus)) / viscosityW;
            dLambdaWdS -= MaterialLaw::krw(problem_.spatialParameters().materialLawParams(globalPos,
                    (*intersection.outside())), std::abs(satMinus)) / viscosityW;
            dLambdaWdS /= (dS);

            if (potentialNW >= 0)
            {
                satUpw = std::max(1 - satI, 0.0);
            }
            else
            {
                satUpw = std::max(1 - satJ, 0.0);
            }

            dS = eps_;

            satPlus = satUpw + eps_;
            satMinus = satUpw;
            if (satMinus - eps_ > 0.0)
            {
                satMinus -= eps_;
                dS += eps_;
            }

            Scalar dLambdaNWdS = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(globalPos,
                    (*intersection.outside())), satPlus) / viscosityNW;
            dLambdaNWdS -= MaterialLaw::krn(problem_.spatialParameters().materialLawParams(globalPos,
                    (*intersection.outside())), satMinus) / viscosityNW;
            dLambdaNWdS /= (dS);

            Scalar lambdaWCap = 0.5 * (lambdaWI + lambdaWJ);
            Scalar lambdaNWCap = 0.5 * (lambdaNWI + lambdaNWJ);

            if (phaseIdx == wPhaseIdx)
            {
                if (flux != 0)
                {
                    cflFluxFunction_ += lambdaNW * dLambdaWdS * std::abs(flux) / (lambdaW * (lambdaW + lambdaNW));
                }

                cflFluxFunction_ -= transmissibility * lambdaWCap * lambdaNWCap * (dPcdSI + dPcdSJ) / (lambdaW
                        + lambdaNW);
            }
            else if (phaseIdx == nPhaseIdx)
            {
                if (flux != 0)
                {
                    cflFluxFunction_ -= lambdaW * dLambdaNWdS * std::abs(flux) / (lambdaNW * (lambdaW + lambdaNW));
                }
            }
            else
            {
                if (flux != 0)
                {
                    switch (saturationType_)
                    {
                    case Sw:
                    {
                    cflFluxFunction_ += dLambdaWdS / (dLambdaWdS + dLambdaNWdS) * std::abs(flux);
                    break;
                    }
                    case Sn:
                    {
                    cflFluxFunction_ +=  dLambdaNWdS / (dLambdaWdS + dLambdaNWdS) * std::abs(flux);
                    break;
                    }
                    }
                }
            }
        }
        else
        {
            // center of face in global coordinates
            GlobalPosition globalPosFace = intersection.geometry().center();

            //get boundary type
            BoundaryConditions::Flags bcTypeSat = problem_.bctypeSat(globalPosFace, intersection);

            // distance vector between barycenters
            Dune::FieldVector < Scalar, dimWorld > distVec = globalPosFace - globalPos;

            // compute distance between cell centers
            Scalar dist = distVec.two_norm();

            //multiply with normal vector at the boundary
            Dune::FieldVector < Scalar, dim > permeability(0);
            permeabilityI.mv(unitOuterNormal, permeability);

            Scalar satBound = 0;
            if (bcTypeSat == BoundaryConditions::dirichlet)
            {
                satBound = problem_.dirichletSat(globalPosFace, intersection);
            }
            else
            {
                satBound = problem_.variables().saturation()[globalIdxI];
            }

            //determine phase saturations from primary saturation variable
            Scalar satW;
            Scalar satNW;
            switch (saturationType_)
            {
            case Sw:
            {
                satW = satBound;
                satNW = 1 - satBound;
                break;
            }
            case Sn:
            {
                satW = 1 - satBound;
                satNW = satBound;
                break;
            }
            default:
            {
                DUNE_THROW(Dune::RangeError, "saturation type not implemented");
            }
            }

            Scalar dPcdSBound = MaterialLaw::dpC_dSw(problem_.spatialParameters().materialLawParams(globalPos,
                    *(intersection.inside())), satBound);

            Scalar lambdaWBound = 0;
            Scalar lambdaNWBound = 0;

            Scalar temperature = problem_.temperature(globalPosFace, (*intersection.inside()));
            Scalar referencePressure = problem_.referencePressure(globalPos, (*intersection.inside()));
            FluidState fluidState;
            Scalar viscosityWBound = FluidSystem::phaseViscosity(wPhaseIdx, temperature, referencePressure, fluidState);
            Scalar viscosityNWBound =
                    FluidSystem::phaseViscosity(nPhaseIdx, temperature, referencePressure, fluidState);
            lambdaWBound = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(globalPos,
                    (*intersection.inside())), satW) / viscosityWBound;
            lambdaNWBound = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(globalPos,
                    (*intersection.inside())), satW) / viscosityNWBound;

            Scalar potentialW = 0;
            Scalar potentialNW = 0;

            potentialW = problem_.variables().potentialWetting(globalIdxI, indexInInside);
            potentialNW = problem_.variables().potentialNonwetting(globalIdxI, indexInInside);

            Scalar transmissibility = (unitOuterNormal * permeability) * intersection.geometry().volume() / dist;

            Scalar satUpw = 0;
            if (potentialW >= 0)
            {
                satUpw = std::max(satI, 0.0);
            }
            else
            {
                satUpw = std::max(satBound, 0.0);
            }

            Scalar dS = eps_;

            Scalar satPlus = satUpw + eps_;
            Scalar satMinus = satUpw;
            if (satMinus - eps_ > 0.0)
            {
                satMinus -= eps_;
                dS += eps_;
            }

            Scalar dLambdaWdS = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(globalPos,
                    (*intersection.inside())), satPlus) / viscosityW;
            dLambdaWdS -= MaterialLaw::krw(problem_.spatialParameters().materialLawParams(globalPos,
                    (*intersection.inside())), satMinus) / viscosityW;
            dLambdaWdS /= (dS);

            if (potentialNW >= 0)
            {
                satUpw = std::max(1 - satI, 0.0);
            }
            else
            {
                satUpw = std::max(1 - satBound, 0.0);
            }

            dS = eps_;

            satPlus = satUpw + eps_;
            satMinus = satUpw;
            if (satMinus - eps_ > 0.0)
            {
                satMinus -= eps_;
                dS += eps_;
            }

            Scalar dLambdaNWdS = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(globalPos,
                    (*intersection.inside())), satPlus) / viscosityNW;
            dLambdaNWdS -= MaterialLaw::krn(problem_.spatialParameters().materialLawParams(globalPos,
                    (*intersection.inside())), satMinus) / viscosityNW;
            dLambdaNWdS /= (dS);

            Scalar lambdaWCap = 0.5 * (lambdaWI + lambdaWBound);
            Scalar lambdaNWCap = 0.5 * (lambdaNWI + lambdaNWBound);

            if (phaseIdx == wPhaseIdx)
            {
                if (flux != 0)
                {
                    cflFluxFunction_ += lambdaNW * dLambdaWdS * std::abs(flux) / (lambdaW * (lambdaW + lambdaNW));
                }

                cflFluxFunction_ -= transmissibility * lambdaWCap * lambdaNWCap * (dPcdSI + dPcdSBound) / (lambdaW
                        + lambdaNW);
            }
            else if (phaseIdx == nPhaseIdx)
            {
                if (flux != 0)
                {
                    cflFluxFunction_ -= lambdaW * dLambdaNWdS * std::abs(flux) / (lambdaNW * (lambdaW + lambdaNW));
                }
            }
            else
            {
                if (flux != 0)
                {
                    switch (saturationType_)
                    {
                    case Sw:
                    {
                    cflFluxFunction_ += dLambdaWdS / (dLambdaWdS + dLambdaNWdS) * std::abs(flux);
                    break;
                    }
                    case Sn:
                    {
                    cflFluxFunction_ +=  dLambdaNWdS / (dLambdaWdS + dLambdaNWdS) * std::abs(flux);
                    break;
                    }
                    }
                }
            }
        }
    }

    Scalar getCflFluxFunction(const GlobalPosition& globalPos, const Element& element)
    {
        Scalar cflFluxDefault = 1 / ParentType::getCflFluxFunction(globalPos, element);

        if (std::isnan(cflFluxFunction_) || std::isinf(cflFluxFunction_) || cflFluxFunction_ > 100 * cflFluxDefault)
        {
            cflFluxFunction_ = 0;
        }

        if (cflFluxDefault > 1)
            std::cout << "cflFluxDefault = " << cflFluxDefault << "\n";

        if (cflFluxFunction_ > 1)
            std::cout << "cflFluxFunction_ = " << cflFluxFunction_ << "\n";

        Scalar returnValue = std::max(cflFluxFunction_, cflFluxDefault);
        reset();

        if (returnValue > 0)
        {
            return (returnValue == cflFluxDefault) ? 0.95 / returnValue : 1 / returnValue;
        }
        else
            return 1e100;
    }

    EvalCflFluxCoats(Problem& problem) :
        ParentType(problem), problem_(problem)
    {
        reset();
    }

private:
    void reset()
    {
        ParentType::reset();
        cflFluxFunction_ = 0;
    }

    Problem& problem_;//problem data
    Scalar cflFluxFunction_;
    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PTAG(PressureFormulation));
    static const int velocityType_ = GET_PROP_VALUE(TypeTag, PTAG(VelocityFormulation));
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, PTAG(SaturationFormulation));
    static const Scalar eps_ = 1e-2;
};

}

#endif
