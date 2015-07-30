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
#include <dumux/implicit/2pnc/2pncfluxvariables.hh>
#include "2pncminproperties.hh"

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
    typedef TwoPNCFluxVariables<TypeTag> ParentType;
    typedef TwoPNCMinFluxVariables<TypeTag> ThisType;
    
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;

    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),
    };

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename FVElementGeometry::SubControlVolume SCV;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    typedef Dune::FieldVector<CoordScalar, dimWorld> DimVector;
    typedef Dune::FieldMatrix<CoordScalar, dim, dim> DimMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,
        wCompIdx  = FluidSystem::wCompIdx,
    };

public:
    /*!
     * \brief The constructor
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the fully implicit scheme
     * \param fIdx The local index of the sub-control-volume face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary Evaluate flux at inner sub-control-volume face or on a boundary face
     */
    TwoPNCMinFluxVariables(const Problem &problem,
                     const Element &element,
                     const FVElementGeometry &fvGeometry,
                     const int fIdx,
                     const ElementVolumeVariables &elemVolVars,
                     const bool onBoundary = false)
    : ParentType(problem, element, fvGeometry, fIdx, elemVolVars, onBoundary) 
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            this->density_[phaseIdx] = Scalar(0);
            this->molarDensity_[phaseIdx] = Scalar(0);
            this->potentialGrad_[phaseIdx] = Scalar(0);
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
            	this->massFractionGrad_[phaseIdx][compIdx] = Scalar(0);
            	this->moleFractionGrad_[phaseIdx][compIdx] = Scalar(0);
            }
        }
        this->calculateGradients_(problem, element, elemVolVars);
        this->calculateVelocities_(problem, element, elemVolVars);
        this->calculateporousDiffCoeff_(problem, element, elemVolVars);
    };

protected:    
    void calculateVelocities_(const Problem &problem,
                              const Element &element,
                              const ElementVolumeVariables &elemVolVars)
    {
        const SpatialParams &spatialParams = problem.spatialParams();
        // multiply the pressure potential with the intrinsic permeability
        DimMatrix K(0.0);

        for (int phaseIdx=0; phaseIdx < numPhases; phaseIdx++)
        {
            const VolumeVariables &volVarsI = elemVolVars[this->face().i];
            const VolumeVariables &volVarsJ = elemVolVars[this->face().j];

            auto K_i = spatialParams.intrinsicPermeability(element,this->fvGeometry_,this->face().i);
            K_i *= volVarsI.permeabilityFactor();

            auto K_j = spatialParams.intrinsicPermeability(element,this->fvGeometry_,this->face().j);
            K_j *= volVarsJ.permeabilityFactor();

            spatialParams.meanK(K,K_i,K_j);

            K.mv(this->potentialGrad_[phaseIdx], this->Kmvp_[phaseIdx]);
            this->KmvpNormal_[phaseIdx] = - (this->Kmvp_[phaseIdx] * this->face().normal);
        }

        // set the upstream and downstream vertices
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            this->upstreamIdx_[phaseIdx] = this->face().i;
            this->downstreamIdx_[phaseIdx] = this->face().j;

            if (this->KmvpNormal_[phaseIdx] < 0) {
                std::swap(this->upstreamIdx_[phaseIdx],
                          this->downstreamIdx_[phaseIdx]);
            }
        }
    }
};

} // end namespace

#endif
