/*****************************************************************************
 *   Copyright (C) 2009-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
 * \brief This file contains the data which is required to calculate
 *        all fluxes of components over a face of a finite volume.
 *
 * This means pressure, concentration and temperature gradients, phase
 * densities at the integration point, etc.
 */
#ifndef DUMUX_MPNC_FLUX_VARIABLES_HH
#define DUMUX_MPNC_FLUX_VARIABLES_HH

#include <dumux/common/spline.hh>

#include "diffusion/fluxvariables.hh"
#include "energy/MpNcfluxvariablesenergy.hh"

namespace Dumux
{

/*!
 * \ingroup MPNCModel
 * \ingroup BoxFluxVariables
 * \brief This template class contains the data which is required to
 *        calculate all fluxes of components over a face of a finite
 *        volume for the two-phase, three-component model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the intergration point, etc.
 */
template <class TypeTag>
class MPNCFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;

    enum {
        dimWorld = GridView::dimensionworld,
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
        numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)),

        enableDiffusion = GET_PROP_VALUE(TypeTag, PTAG(EnableDiffusion)),
        enableEnergy = GET_PROP_VALUE(TypeTag, PTAG(EnableEnergy)),
        enableKinetic = GET_PROP_VALUE(TypeTag, PTAG(EnableKinetic)),
        enableKineticEnergy = GET_PROP_VALUE(TypeTag, PTAG(EnableKineticEnergy)),
        enableGravity = GET_PROP_VALUE(TypeTag, PTAG(EnableGravity)),
    };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> Vector;
    typedef Dune::FieldMatrix<CoordScalar, dimWorld, dimWorld> Tensor;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolume SCV;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;

    typedef MPNCFluxVariablesDiffusion<TypeTag, enableDiffusion>                  FluxVariablesDiffusion;
    typedef MPNCFluxVariablesEnergy<TypeTag, enableEnergy, enableKineticEnergy>    FluxVariablesEnergy;

public:
    MPNCFluxVariables(const Problem &problem,
                   const Element &element,
                   const FVElementGeometry &elemGeom,
                   int scvfIdx,
                   const ElementVolumeVariables &elemVolVars)
        : fvElemGeom_(elemGeom), volVars_(elemVolVars)
    {
        scvfIdx_ = scvfIdx;

        // update the base module (i.e. advection)
        calculateGradients_(problem, element, elemVolVars);
        calculateVelocities_(problem, element, elemVolVars);

        // update the flux data of the energy module (i.e. isothermal
        // or non-isothermal)
        energyDat_.update(problem, element, elemGeom, scvfIdx, *this, elemVolVars);

        // update the flux data of the diffusion module (i.e. with or
        // without diffusion)
        diffusionDat_.update(problem, element, elemGeom, scvfIdx, elemVolVars);

        extrusionFactor_ =
            (elemVolVars[face().i].extrusionFactor() 
             + elemVolVars[face().j].extrusionFactor()) / 2;
    }


    /*!
     * \brief Calculate a phase's darcy velocity [m/s] at a
     *        sub-control volume face.
     *
     * So far, this method only exists in the Mp-Nc model!
     *
     *  Of course, in the setting of a finite volume scheme, the velocities are
     *  on the faces rather than in the volume. Therefore, the velocity
     *
     * \param vDarcy the resulting Darcy velocity
     * \param elemVolVars element volume variables
     * \param phaseIdx phase index
     */
    void computeDarcy(Vector & vDarcy,
                      const ElementVolumeVariables &elemVolVars,
                      int phaseIdx) const
    {
        intrinsicPermeability().mv(potentialGrad(phaseIdx),
                                            vDarcy);
        // darcy velocity is along *negative* potential gradient
        // (i.e. from high to low pressures), this means that we need
        // to negate the product of the intrinsic permeability and the
        // potential gradient!
        vDarcy *= -1;

        // JUST for upstream decision
        Scalar normalFlux = vDarcy * face().normal;
        // data attached to upstream and the downstream vertices
        // of the current phase
        int upIdx = face().i;
        int dnIdx = face().j;

        if (!std::isfinite(normalFlux))
            DUNE_THROW(NumericalProblem, "Calculated non-finite normal flux");

        if (normalFlux < 0)
            std::swap(upIdx, dnIdx);

        const VolumeVariables &up = elemVolVars[upIdx];

        ////////
        // Jipie this is a velocity now, finally deserves the name
        ////////
        vDarcy *= up.mobility(phaseIdx);
    }

    const SCVFace &face() const
    { return fvElemGeom_.subContVolFace[scvfIdx_]; }

    const VolumeVariables &volVars(int idx) const
    { return volVars_[idx]; }

    /*!
     * \brief Returns th extrusion factor for the sub-control volume face
     */
    Scalar extrusionFactor() const
    { return extrusionFactor_; }

    /*!
     * \brief Return the intrinsic permeability.
     */
    const Tensor &intrinsicPermeability() const
    { return K_; }

    /*!
     * \brief Return the pressure potential gradient.
     */
    const Vector &potentialGrad(int phaseIdx) const
    { return potentialGrad_[phaseIdx]; }

    ////////////////////////////////////////////////
    // forward calls to the diffusion module
    Scalar porousDiffCoeffL(int compIdx) const
    { return diffusionDat_.porousDiffCoeffL(compIdx); };

    Scalar porousDiffCoeffG(int compIIdx, int compJIdx) const
    { return diffusionDat_.porousDiffCoeffG(compIIdx, compJIdx); };

    const Scalar moleFrac(int phaseIdx, int compIdx) const
    { return diffusionDat_.moleFrac(phaseIdx, compIdx); };

    const Vector &moleFracGrad(int phaseIdx,
                               int compIdx) const
    { return diffusionDat_.moleFracGrad(phaseIdx, compIdx); };
    // end of forward calls to the diffusion module
    ////////////////////////////////////////////////

    ////////////////////////////////////////////////
    // forward calls to the temperature module
    const Vector &temperatureGrad() const
    { return energyDat_.temperatureGrad(); };

    const FluxVariablesEnergy &energyData() const
    { return energyDat_; }
    // end of forward calls to the temperature module
    ////////////////////////////////////////////////

private:
    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemVolVars)
    {
        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++) {
            potentialGrad_[phaseIdx] = Scalar(0);
        }

        // calculate pressure gradients using finite element gradients
        Vector tmp(0.0);
        for (int idx = 0;
             idx < fvElemGeom_.numVertices;
             idx++) // loop over adjacent vertices
        {
            // FE gradient at vertex idx
            const Vector &feGrad = face().grad[idx];

            // TODO: only calculate the gradients for the present
            // phases.
            //
            // compute sum of pressure gradients for each phase
            for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
            {
                // the pressure gradient
                tmp = feGrad;
                tmp *= elemVolVars[idx].fluidState().pressure(phaseIdx);
                potentialGrad_[phaseIdx] += tmp;
            }
        }

        ///////////////
        // correct the pressure gradients by the gravitational acceleration
        ///////////////
        if (enableGravity) {
            // estimate the gravitational acceleration at a given SCV face
            // using the arithmetic mean
            Vector g(problem.boxGravity(element, fvElemGeom_, face().i));
            g += problem.boxGravity(element, fvElemGeom_, face().j);
            g /= 2;

            for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
            {
                // calculate the phase density at the integration point. we
                // only do this if the wetting phase is present in both cells
                Scalar SI = elemVolVars[face().i].fluidState().saturation(phaseIdx);
                Scalar SJ = elemVolVars[face().j].fluidState().saturation(phaseIdx);
                Scalar rhoI = elemVolVars[face().i].fluidState().density(phaseIdx);
                Scalar rhoJ = elemVolVars[face().j].fluidState().density(phaseIdx);
                Scalar fI = std::max(0.0, std::min(SI/1e-5, 0.5));
                Scalar fJ = std::max(0.0, std::min(SJ/1e-5, 0.5));
                if (fI + fJ == 0)
                    // doesn't matter because no wetting phase is present in
                    // both cells!
                    fI = fJ = 0.5;
                Scalar density = (fI*rhoI + fJ*rhoJ)/(fI + fJ);
                
                // make gravity acceleration a force
                Vector f(g);
                f *= density;
        
                // calculate the final potential gradient
                potentialGrad_[phaseIdx] -= f;
            }
        }
    }

    void calculateVelocities_(const Problem &problem,
                              const Element &element,
                              const ElementVolumeVariables &elemVolVars)
    {
        // multiply the pressure potential with the intrinsic
        // permeability
        const SpatialParameters &sp = problem.spatialParameters();
        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
        {
            sp.meanK(K_,
                     sp.intrinsicPermeability(element,
                                              fvElemGeom_,
                                              face().i),
                     sp.intrinsicPermeability(element,
                                              fvElemGeom_,
                                              face().j));
        }
    }




    const FVElementGeometry &fvElemGeom_;
    int scvfIdx_;

    const ElementVolumeVariables &volVars_;

    // The extrusion factor for the sub-control volume face
    Scalar extrusionFactor_;

    // pressure potential gradients
    Vector potentialGrad_[numPhases];

    // intrinsic permeability
    Tensor K_;

    FluxVariablesDiffusion     diffusionDat_;
    FluxVariablesEnergy     energyDat_;
};

} // end namepace

#endif
