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
 * \brief Contains the quantities to calculate the energy flux in the
 *        MpNc fully implicit model.
 */
#ifndef DUMUX_MPNC_ENERGY_FLUX_VARIABLES_HH
#define DUMUX_MPNC_ENERGY_FLUX_VARIABLES_HH

#include <dune/common/fvector.hh>

#include <dumux/implicit/mpnc/mpncproperties.hh>
#include <dumux/common/spline.hh>

namespace Dumux
{

/*!
 * \ingroup MPNCModel
 * \ingroup ImplicitFluxVariables
 * \brief Variables for the enthalpy fluxes in the MpNc model
 */
template <class TypeTag, bool enableEnergy/*=false*/, bool kineticEnergyTransfer/*=false*/>
class MPNCFluxVariablesEnergy
{
    static_assert(!(kineticEnergyTransfer && !enableEnergy),
                  "No kinetic energy transfer may only be enabled "
                  "if energy is enabled in general.");
    static_assert(!kineticEnergyTransfer,
                  "No kinetic energy transfer module included, "
                  "but kinetic energy transfer enabled.");

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

public:
    /*!
     * \brief The constructor
     */
    MPNCFluxVariablesEnergy()
    {
    }
    /*!
     * \brief update
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param face The SCV (sub-control-volume) face
     * \param fluxVars The flux variables
     * \param elemVolVars The volume variables of the current element
     */
    void update(const Problem & problem,
                const Element & element,
                const FVElementGeometry & fvGeometry,
                const SCVFace & face,
                const FluxVariables & fluxVars,
                const ElementVolumeVariables & elemVolVars)
    {}
};

template <class TypeTag>
class MPNCFluxVariablesEnergy<TypeTag, /*enableEnergy=*/true,  /*kineticEnergyTransfer=*/false>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename FluidSystem::ParameterCache ParameterCache;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum{dim = GridView::dimension};
    enum{dimWorld = GridView::dimensionworld};
    enum{nPhaseIdx = FluidSystem::nPhaseIdx};
    enum{wPhaseIdx = FluidSystem::wPhaseIdx};
    enum{numPhases = GET_PROP_VALUE(TypeTag, NumPhases)};

    typedef Dune::FieldVector<CoordScalar, dimWorld>  DimVector;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

public:
    /*!
     * \brief The constructor
     */
    MPNCFluxVariablesEnergy()
    {}
    /*!
     * \brief update
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param face The SCV (sub-control-volume) face
     * \param fluxVars The flux variables
     * \param elemVolVars The volume variables of the current element
     */
    void update(const Problem & problem,
                const Element & element,
                const FVElementGeometry & fvGeometry,
                const SCVFace & face,
                const FluxVariables & fluxVars,
                const ElementVolumeVariables & elemVolVars)
    {
        // calculate temperature gradient using finite element
        // gradients
        DimVector tmp(0.0);
        DimVector temperatureGradient(0.);
        for (int idx = 0; idx < face.numFap; idx++)
        {
            tmp = face.grad[idx];

            // index for the element volume variables 
            int volVarsIdx = face.fapIndices[idx];
            
            tmp *= elemVolVars[volVarsIdx].fluidState().temperature(/*phaseIdx=*/0);
            temperatureGradient += tmp;
        }

        // project the heat flux vector on the face's normal vector
        temperatureGradientNormal_ = temperatureGradient * face.normal;


        lambdaPm_ = lumpedLambdaPm(problem,
                                   element,
                                   fvGeometry,
                                   face,
                                   elemVolVars) ;

    }
    /*!
     * \brief The lumped / average conductivity of solid plus phases \f$[W/mK]\f$.
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param face The SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     */
    Scalar lumpedLambdaPm(const Problem &problem,
                          const Element &element,
                          const FVElementGeometry & fvGeometry,
                          const SCVFace & face,
                          const ElementVolumeVariables & elemVolVars)
    {
         // arithmetic mean of the liquid saturation and the porosity
         const unsigned int i = face.i;
         const unsigned int j = face.j;

         const FluidState &fsI = elemVolVars[i].fluidState();
         const FluidState &fsJ = elemVolVars[j].fluidState();
         const Scalar swi = fsI.saturation(wPhaseIdx);
         const Scalar swj = fsJ.saturation(wPhaseIdx);

         typename FluidSystem::ParameterCache paramCacheI, paramCacheJ;
         paramCacheI.updateAll(fsI);
         paramCacheJ.updateAll(fsJ);

         const Scalar sw = std::max<Scalar>(0.0, 0.5*(swi + swj));

         //        const Scalar lambdaDry = 0.583; // W / (K m) // works, orig
         //        const Scalar lambdaWet = 1.13; // W / (K m) // works, orig

         const Scalar lambdaSoilI = problem.spatialParams().soilThermalConductivity(element, fvGeometry, i);
         const Scalar lambdaSoilJ = problem.spatialParams().soilThermalConductivity(element, fvGeometry, i);
         const Scalar lambdaDry = 0.5 * (lambdaSoilI + FluidSystem::thermalConductivity(fsI, paramCacheI, nPhaseIdx)); // W / (K m)
         const Scalar lambdaWet = 0.5 * (lambdaSoilJ + FluidSystem::thermalConductivity(fsJ, paramCacheJ, wPhaseIdx)) ; // W / (K m)

         // the heat conductivity of the matrix. in general this is a
         // tensorial value, but we assume isotropic heat conductivity.
         // This is the Sommerton approach with lambdaDry =
         // lambdaSn100%.  Taken from: H. Class: "Theorie und
         // numerische Modellierung nichtisothermer Mehrphasenprozesse
         // in NAPL-kontaminierten poroesen Medien", PhD Thesis, University of
         // Stuttgart, Institute of Hydraulic Engineering, p. 57

         Scalar result;
         if (sw < 0.1) {
             // regularization
             Dumux::Spline<Scalar> sp(0, 0.1, // x1, x2
                                      0, sqrt(0.1), // y1, y2
                                      5*0.5/sqrt(0.1), 0.5/sqrt(0.1)); // m1, m2
             result = lambdaDry + sp.eval(sw)*(lambdaWet - lambdaDry);
         }
         else
             result = lambdaDry + std::sqrt(sw)*(lambdaWet - lambdaDry);

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
