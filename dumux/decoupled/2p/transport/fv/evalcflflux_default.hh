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
#ifndef DUMUX_EVALCFLFLUX_DEFAULT_HH
#define DUMUX_EVALCFLFLUX_DEFAULT_HH

/**
 * @file
 * @brief  Fluxes to evaluate a CFL-Condition
 * @author Markus Wolff
 */

#include "evalcflflux.hh"

namespace Dumux
{
/*!\ingroup Saturation2p
 * @brief  Default implementation of cfl-fluxes to evaluate a CFL-Condition
 *
 * Compares the maximum of inflow and outflow to the element volume weighted by relative permeability and viscosity ratios.
 *
 * Template parameters are:

 - TypeTag PropertyTag of the problem implementation
 */
template<class TypeTag>
class EvalCflFluxDefault: public EvalCflFlux<TypeTag>
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
      typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
      typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

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

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::Intersection Intersection;

public:

    void addFlux(Scalar& lambdaW, Scalar& lambdaNW, Scalar& viscosityW, Scalar& viscosityNW, Scalar flux, const Element& element, int phaseIdx = -1)
    {
        addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, flux, phaseIdx);
    }

    void addFlux(Scalar& lambdaW, Scalar& lambdaNW, Scalar& viscosityW, Scalar& viscosityNW, Scalar flux, const Intersection& intersection, int phaseIdx = -1)
    {
        addFlux(lambdaW, lambdaNW, viscosityW, viscosityNW, flux, phaseIdx);
    }

    Scalar getCflFluxFunction(const GlobalPosition& globalPos, const Element& element);

    EvalCflFluxDefault (Problem& problem)
    : problem_(problem)
    {
        reset();
    }

private:

    void addFlux(Scalar& lambdaW, Scalar& lambdaNW, Scalar& viscosityW, Scalar& viscosityNW, Scalar flux, int phaseIdx = -1)
    {
        Scalar krSum = lambdaW * viscosityW + lambdaNW * viscosityNW;
        Scalar viscosityRatio = 1 - fabs(0.5 - viscosityNW / (viscosityW + viscosityNW));//1 - fabs(viscosityWI-viscosityNWI)/(viscosityWI+viscosityNWI);

        switch (phaseIdx)
         {
         case wPhaseIdx:
         {
             //for time step criterion
             if (flux >= 0)
             {
                 fluxWettingOut_ += flux / (krSum * viscosityRatio);
             }
             if (flux < 0)
             {
                 fluxIn_ -= flux / (krSum * viscosityRatio);
             }

             break;
         }

             //for time step criterion if the non-wetting phase velocity is used
         case nPhaseIdx:
         {
             if (flux >= 0)
             {
                 fluxNonwettingOut_ += flux / (krSum * viscosityRatio);
             }
             if (flux < 0)
             {
                 fluxIn_ -= flux / (krSum * viscosityRatio);
             }

             break;
         }
         default:
         {
             if (flux >= 0)
             {
                 fluxOut_ += flux / (krSum * viscosityRatio);
             }
             if (flux < 0)
             {
                 fluxIn_ -= flux / (krSum * viscosityRatio);
             }

             break;
         }
         }
    }

    Scalar getCFLFluxIn(int phaseIdx = 0)
    {
        if (std::isnan(fluxIn_) || std::isinf(fluxIn_))
        {
            fluxIn_ = 1e-100;
        }

            return fluxIn_;
    }

    Scalar getCFLFluxOut(int phaseIdx = 0)
    {
        if (std::isnan(fluxWettingOut_) || std::isinf(fluxWettingOut_))
        {
            fluxWettingOut_ = 1e-100;
        }

        if (std::isnan(fluxNonwettingOut_) || std::isinf(fluxNonwettingOut_))
        {
            fluxNonwettingOut_ = 1e-100;
        }
        if (std::isnan(fluxOut_) || std::isinf(fluxOut_))
        {
            fluxOut_ = 1e-100;
        }

        if (phaseIdx == wPhaseIdx)
            return fluxWettingOut_;
        else if (phaseIdx == nPhaseIdx)
            return fluxNonwettingOut_;
        else
            return fluxOut_;
    }


protected:

    //! resets the accumulated CFL-fluxes to zero
    void reset()
    {
        fluxWettingOut_ = 0;
        fluxNonwettingOut_ = 0;
        fluxIn_ = 0;
        fluxOut_ = 0;
    }

private:
    Problem& problem_;//problem data
    Scalar fluxWettingOut_;
    Scalar fluxNonwettingOut_;
    Scalar fluxOut_;
    Scalar fluxIn_;
    static const int velocityType_ = GET_PROP_VALUE(TypeTag, PTAG(VelocityFormulation));
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, PTAG(SaturationFormulation));
};

template<class TypeTag>
typename EvalCflFluxDefault<TypeTag>::Scalar EvalCflFluxDefault<TypeTag>::getCflFluxFunction(const GlobalPosition& globalPos, const Element& element)
{
    Scalar residualSatW = problem_.spatialParameters().materialLawParams(globalPos, element).Swr();
    Scalar residualSatNW = problem_.spatialParameters().materialLawParams(globalPos, element).Snr();

    // compute dt restriction
    Scalar volumeCorrectionFactor = 1 - residualSatW - residualSatNW;
    Scalar volumeCorrectionFactorOutW = 0;
    Scalar volumeCorrectionFactorOutNW = 0;

    Scalar sat = problem_.variables().saturation()[problem_.variables().index(element)];

    if (saturationType_ == Sw)
    {
        volumeCorrectionFactorOutW = std::max((sat - residualSatW), 1e-2);
        volumeCorrectionFactorOutNW = std::max((1 - sat - residualSatNW), 1e-2);
    }
    if (saturationType_ == Sn)
    {
        volumeCorrectionFactorOutNW = std::max((1 - sat - residualSatNW), 1e-2);
        volumeCorrectionFactorOutNW = std::max((sat - residualSatNW), 1e-2);
    }

    //make sure correction is in the right range. If not: force dt to be not min-dt!
    if (volumeCorrectionFactorOutW <= 0)
    {
        volumeCorrectionFactorOutW = 1e100;
    }
    if (volumeCorrectionFactorOutNW <= 0)
    {
        volumeCorrectionFactorOutNW = 1e100;
    }

    //correct volume
    Scalar cFLFluxIn = volumeCorrectionFactor / getCFLFluxIn();
    Scalar cFLFluxOut = 0;

    if (velocityType_ == vt)
    {
        cFLFluxOut = volumeCorrectionFactor / getCFLFluxOut(-1);
    }
    else
    {
        cFLFluxOut = std::min(volumeCorrectionFactorOutW / getCFLFluxOut(wPhaseIdx), volumeCorrectionFactorOutNW / getCFLFluxOut(nPhaseIdx));
    }

    //determine timestep
    Scalar cFLFluxFunction = std::min(cFLFluxIn, cFLFluxOut);

    reset();

    return cFLFluxFunction;
}

}

#endif
