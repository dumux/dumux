/*****************************************************************************
 *   Copyright (C) 2008,2009 by Andreas Lauser                               *
 *   Copyright (C) 2008,2009 by Melanie Darcis                               *
 *                                                                           *
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
/*!
 * \file
 *
 * \brief Contains the quantities to calculate the energy flux in the
 *        MpNc box model.
 */
#ifndef DUMUX_MPNC_ENERGY_FLUX_VARIABLES_HH
#define DUMUX_MPNC_ENERGY_FLUX_VARIABLES_HH

namespace Dumux
{

template <class TypeTag, bool enableEnergy/*=false*/, bool kineticEnergyTransfer/*=false*/>
class MPNCFluxVariablesEnergy
{
    static_assert(!(kineticEnergyTransfer && !enableEnergy),
                  "No kinetic energy transfer may only be enabled "
                  "if energy is enabled in general.");
    static_assert(!kineticEnergyTransfer,
                  "No kinetic energy transfer module included, "
                  "but kinetic energy transfer enabled.");

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;

public:
    MPNCFluxVariablesEnergy()
    {
    }

    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvElemGeom,
                int scvfIdx,
                const FluxVariables &fluxDat,
                const ElementVolumeVariables &elemVolVars)
    {};
};

template <class TypeTag>
class MPNCFluxVariablesEnergy<TypeTag, /*enableEnergy=*/true,  /*kineticEnergyTransfer=*/false>
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        gPhaseIdx = FluidSystem::gPhaseIdx,
        lPhaseIdx = FluidSystem::lPhaseIdx,
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
    };

    typedef Dune::FieldVector<CoordScalar, dimWorld>  Vector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolume SCV;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

public:
    MPNCFluxVariablesEnergy()
    {}

    void update(const  Problem & problem,
                const Element &element,
                const FVElementGeometry &fvElemGeom,
                int scvfIdx,
                const FluxVariables &fluxDat,
                const ElementVolumeVariables &elemVolVars)
    {
        // calculate temperature gradient using finite element
        // gradients
        Vector tmp(0.0);
        Vector temperatureGradient(0.);
        for (int scvIdx = 0; scvIdx < fvElemGeom.numVertices; scvIdx++)
        {
            tmp = fvElemGeom.subContVolFace[scvfIdx].grad[scvIdx];
            tmp *= elemVolVars[scvIdx].fluidState().temperature(/*phaseIdx=*/0);
            temperatureGradient += tmp;
        }

        // project the heat flux vector on the face's normal vector
        temperatureGradientNormal_ = temperatureGradient * fvElemGeom.subContVolFace[scvfIdx].normal;


        lambdaPm_ = lumpedLambdaPm(problem,
                                   element,
                                    fvElemGeom,
                                    scvfIdx,
                                    elemVolVars) ;

    }

    Scalar lumpedLambdaPm(const Problem &problem,
                          const Element &element,
                           const FVElementGeometry & fvElemGeom,
                           int faceIdx,
                           const ElementVolumeVariables & elemVolVars)
    {
         // arithmetic mean of the liquid saturation and the porosity
         const int i = fvElemGeom.subContVolFace[faceIdx].i;
         const int j = fvElemGeom.subContVolFace[faceIdx].j;

         typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables))::FluidState FluidState;
         const FluidState &fsI = elemVolVars[i].fluidState();
         const FluidState &fsJ = elemVolVars[j].fluidState();
         const Scalar Sli = fsI.saturation(lPhaseIdx);
         const Scalar Slj = fsJ.saturation(lPhaseIdx);

         typename FluidSystem::ParameterCache paramCacheI, paramCacheJ;
         paramCacheI.updateAll(fsI);
         paramCacheJ.updateAll(fsJ);

         const Scalar Sl = std::max<Scalar>(0.0, 0.5*(Sli + Slj));

         //        const Scalar lambdaDry = 0.583; // W / (K m) // works, orig
         //        const Scalar lambdaWet = 1.13; // W / (K m) // works, orig
         
         Scalar lambdaSoilI = problem.spatialParameters().soilThermalConductivity(element, fvElemGeom, i);
         Scalar lambdaSoilJ = problem.spatialParameters().soilThermalConductivity(element, fvElemGeom, i);
         const Scalar lambdaDry = 0.5 * (lambdaSoilI + FluidSystem::thermalConductivity(fsI, paramCacheI, gPhaseIdx)); // W / (K m)
         const Scalar lambdaWet = 0.5 * (lambdaSoilJ + FluidSystem::thermalConductivity(fsJ, paramCacheJ, lPhaseIdx)) ; // W / (K m)
         
         // the heat conductivity of the matrix. in general this is a
         // tensorial value, but we assume isotropic heat conductivity.
         // This is the Sommerton approach with lambdaDry =
         // lambdaSn100%.  Taken from: H. Class: "Theorie und
         // numerische Modellierung nichtisothermer Mehrphasenprozesse
         // in NAPL-kontaminierten poroesen Medien", PhD Thesis, University of
         // Stuttgart, Institute of Hydraulic Engineering, p. 57
         
         Scalar result;
         if (Sl < 0.1) {
             // regularization
             Dumux::Spline<Scalar> sp(0, 0.1, // x1, x2
                                      0, sqrt(0.1), // y1, y2
                                      5*0.5/sqrt(0.1), 0.5/sqrt(0.1)); // m1, m2
             result = lambdaDry + sp.eval(Sl)*(lambdaWet - lambdaDry);
         }
         else
             result = lambdaDry + std::sqrt(Sl)*(lambdaWet - lambdaDry);
         
         return result;
    }

    /*!
     * \brief The lumped / average conductivity of solid plus phases \f$[W/mK]\f$.
     */
    Scalar lambdaPm() const
    { return lambdaPm_; }

    /*!
     * \brief The normal of the gradient of temperature .
     */
    Scalar temperatureGradientNormal() const
    {
        return temperatureGradientNormal_;
    }

private:
    Scalar lambdaPm_;
    Scalar temperatureGradientNormal_;
};

} // end namepace

#endif
