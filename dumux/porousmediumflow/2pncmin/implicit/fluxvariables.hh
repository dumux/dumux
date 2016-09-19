// -**- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
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
 * \brief Contains the data which is required to calculate
 *        all fluxes of components over a face of a finite volume for
 *        the two-phase two-component mineralization model fully implicit model.
 */
#ifndef DUMUX_2PNCMIN_FLUX_VARIABLES_HH
#define DUMUX_2PNCMIN_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/common/spline.hh>
#include <dumux/porousmediumflow/2pnc/implicit/fluxvariables.hh>
#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup TwoPNCMinModel
 * \ingroup ImplicitFluxVariables
 * \brief Contains the data which is required to calculate
 *        all fluxes of components over a face of a finite volume for
 *        the two-phase n-component mineralization fully implicit model.
 *
 * This means pressure and concentration gradients, phase densities at
 * the integration point, etc.
 */

template <class TypeTag>
class TwoPNCMinFluxVariables : public TwoPNCFluxVariables<TypeTag>
{
    friend typename GET_PROP_TYPE(TypeTag, BaseFluxVariables); // be friends with base class
    friend class TwoPNCFluxVariables<TypeTag>; // be friends with parent class

    typedef TwoPNCFluxVariables<TypeTag> ParentType;
    typedef TwoPNCMinFluxVariables<TypeTag> ThisType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),
    };

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef Dune::FieldMatrix<CoordScalar, dimWorld, dimWorld> DimWorldMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,
        wCompIdx  = FluidSystem::wCompIdx,
    };

protected:
    /*!
     * \brief Actual calculation of the normal Darcy velocities.
     * \note this overloads the darcy flux variables velocity calculation
     *       This is only necessary because of the permeability factor.
     * \todo Remove this once we have solDependent spatialParams!
     *
     * \param problem The problem
     * \param element The finite element
     * \param elemVolVars The volume variables of the current element
     */
    void calculateNormalVelocity_(const Problem &problem,
                                  const Element &element,
                                  const ElementVolumeVariables &elemVolVars)
    {
        // calculate the mean intrinsic permeability
        const auto& spatialParams = problem.spatialParams();
        DimWorldMatrix K(0.0);
        const auto& volVarsI = elemVolVars[this->face().i];
        const auto& volVarsJ = elemVolVars[this->face().j];

        if (GET_PROP_VALUE(TypeTag, ImplicitIsBox))
        {
            auto Ki = spatialParams.intrinsicPermeability(element, this->fvGeometry_(), this->face().i);
            Ki *= volVarsI.permeabilityFactor();

            auto Kj = spatialParams.intrinsicPermeability(element, this->fvGeometry_(), this->face().j);
            Kj *= volVarsJ.permeabilityFactor();

            spatialParams.meanK(K, Ki, Kj);
        }
        else
        {
            const Element& elementI = this->fvGeometry_().neighbors[this->face().i];
            FVElementGeometry fvGeometryI;
            fvGeometryI.subContVol[0].global = elementI.geometry().center();

            const Element& elementJ = this->fvGeometry_().neighbors[this->face().j];
            FVElementGeometry fvGeometryJ;
            fvGeometryJ.subContVol[0].global = elementJ.geometry().center();

            auto Ki = spatialParams.intrinsicPermeability(elementI, fvGeometryI, 0);
            Ki *= volVarsI.permeabilityFactor();

            auto Kj = spatialParams.intrinsicPermeability(elementJ, fvGeometryJ, 0);
            Kj *= volVarsJ.permeabilityFactor();

            spatialParams.meanK(K, Ki, Kj);
        }

        // loop over all phases
        for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
        {
            // calculate the flux in the normal direction of the
            // current sub control volume face:
            //
            // v = - (K_f grad phi) * n
            // with K_f = rho g / mu K
            //
            // Mind, that the normal has the length of it's area.
            // This means that we are actually calculating
            //  Q = - (K grad phi) dot n /|n| * A


            K.mv(this->potentialGrad_[phaseIdx], this->kGradP_[phaseIdx]);
            this->kGradPNormal_[phaseIdx] = this->kGradP_[phaseIdx]*this->face().normal;

            // determine the upwind direction
            if (this->kGradPNormal_[phaseIdx] < 0)
            {
                this->upstreamIdx_[phaseIdx] = this->face().i;
                this->downstreamIdx_[phaseIdx] = this->face().j;
            }
            else
            {
                this->upstreamIdx_[phaseIdx] = this->face().j;
                this->downstreamIdx_[phaseIdx] = this->face().i;
            }

            // obtain the upwind volume variables
            const auto& upVolVars = elemVolVars[ this->upstreamIdx(phaseIdx) ];
            const auto& downVolVars = elemVolVars[ this->downstreamIdx(phaseIdx) ];

            // the minus comes from the Darcy relation which states that
            // the flux is from high to low potentials.
            // set the velocity
            this->velocity_[phaseIdx] = this->kGradP_[phaseIdx];
            this->velocity_[phaseIdx] *= - ( this->mobilityUpwindWeight_*upVolVars.mobility(phaseIdx)
                    + (1.0 - this->mobilityUpwindWeight_)*downVolVars.mobility(phaseIdx)) ;

            // set the volume flux
            this->volumeFlux_[phaseIdx] = this->velocity_[phaseIdx] * this->face().normal;
        }// loop all phases
    }
};

} // end namespace

#endif
