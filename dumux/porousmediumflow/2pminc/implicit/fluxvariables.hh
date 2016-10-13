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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 * \brief This file contains the data which is required to calculate
 *        all fluxes of fluid phases over a face of a finite volume.
 *
 * This means pressure and temperature gradients, phase densities at
 * the integration point, etc.
 */
#ifndef DUMUX_TWOPMINC_FLUX_VARIABLES_HH
#define DUMUX_TWOPMINC_FLUX_VARIABLES_HH

#include <dune/common/exceptions.hh>

#include "properties.hh"

#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>

namespace Dumux
{

namespace Properties
{
// forward declaration of properties
NEW_PROP_TAG(ImplicitMobilityUpwindWeight);
NEW_PROP_TAG(SpatialParams);
NEW_PROP_TAG(NumPhases);
NEW_PROP_TAG(ProblemEnableGravity);
}

/*!
 * \ingroup ImplicitFluxVariables
 * \brief Evaluates the normal component of the Darcy velocity
 * on a (sub)control volume face.
 */
template <class TypeTag>
class TwoPMincFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum { dimWorld = GridView::dimensionworld} ;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases)} ;
    enum { fractureIdx = 0} ; /*!< index of the fracture (continua are indexed nC = 0,1,..numContinua-1)*/

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

public:
    /*!
     * \brief Compute / update the flux variables
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param fIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary A boolean variable to specify whether the flux variables
     * are calculated for interior SCV faces or boundary faces, default=false
     * \todo The fvGeometry should be better initialized, passed and stored as an std::shared_ptr
     */
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int fIdx,
                const ElementVolumeVariables &elemVolVars,
                const bool onBoundary = false)
    {
        fvGeometryPtr_ = &fvGeometry;
        onBoundary_ = onBoundary;
        faceIdx_ = fIdx;

        mobilityUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MobilityUpwindWeight);
        calculateGradients_(problem, element, elemVolVars);
        calculateNormalVelocity_(problem, element, elemVolVars);
    }

public:
    /*!
     * \brief Return the volumetric flux over a face of a given phase.
     *
     *        This is the calculated velocity multiplied by the unit normal
     *        and the area of the face.
     *        face().normal
     *        has already the magnitude of the area.
     *
     * \param phaseIdx index of the phase
     */
    Scalar volumeFlux(const unsigned int phaseIdx) const
    { return volumeFlux_[phaseIdx]; }

    /*!
     * \brief Return the velocity of a given phase.
     *
     *        This is the full velocity vector on the
     *        face (without being multiplied with normal).
     *
     * \param phaseIdx index of the phase
     */
    DimVector velocity(const unsigned int phaseIdx) const
    { return velocity_[phaseIdx] ; }

    /*!
     * \brief Return intrinsic permeability multiplied with potential
     *        gradient multiplied with normal.
     *        I.e. everything that does not need upwind decisions.
     *
     * \param phaseIdx index of the phase
     */
    Scalar kGradPNormal(const unsigned int phaseIdx) const
    { return kGradPNormal_[phaseIdx] ; }

    /*!
     * \brief Return the local index of the downstream control volume
     *        for a given phase.
     *
     * \param phaseIdx index of the phase
     */
    unsigned int downstreamIdx(const unsigned phaseIdx) const
    { return downstreamIdx_[phaseIdx]; }

    /*!
     * \brief Return the local index of the upstream control volume
     *        for a given phase.
     *
     * \param phaseIdx index of the phase
     */
    unsigned int upstreamIdx(const unsigned phaseIdx) const
    { return upstreamIdx_[phaseIdx]; }

    /*!
     * \brief Return the SCV (sub-control-volume) face. This may be either
     *        a face within the element or a face on the element boundary,
     *        depending on the value of onBoundary_.
     */
    const SCVFace &face() const
    {
        if (onBoundary_)
            return fvGeometry_().boundaryFace[faceIdx_];
        else
            return fvGeometry_().subContVolFace[faceIdx_];
    }

protected:

    /*!
     * \brief Calculation of the potential gradients between neighbored fracture elements
     * Inter-cell flux only occurs through the fractures,
     * i.e., continua 0.
     *
     * \param problem The problem
     * \param element The grid element
     * \param elemVolVars The volume variables of the current element
     * are calculated for interior SCV faces or boundary faces, default=false
     */
    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemVolVars)
    {
        // loop over all phases
        for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
        {
            gradPotential_[phaseIdx]= 0.0 ;
            for (int idx = 0;
                 idx < face().numFap;
                 idx++) // loop over adjacent vertices
            {


                // index for the element volume variables
                int volVarsIdx = face().fapIndices[idx];

                // the pressure gradient
                DimVector gradPressure = face().grad[idx];
                gradPressure *= elemVolVars[volVarsIdx].pressure(phaseIdx, fractureIdx);
                gradPotential_[phaseIdx] += gradPressure;
            }

            // correct the pressure gradient by the gravitational acceleration
            if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
            {

                // calculate the phase density at the integration point. we
                // only do this if the wetting phase is present in both cells
                Scalar si = elemVolVars[face().i].saturation(phaseIdx, fractureIdx);
                Scalar sj = elemVolVars[face().j].saturation(phaseIdx, fractureIdx);
                Scalar rhoi = elemVolVars[face().i].density(phaseIdx, fractureIdx);
                Scalar rhoj = elemVolVars[face().j].density(phaseIdx, fractureIdx);
                Scalar fi = std::max(0.0, std::min(si/1e-5, 0.5));
                Scalar fj = std::max(0.0, std::min(sj/1e-5, 0.5));
                if (fi + fj == 0)
                    // doesn't matter because no wetting phase is present in
                    // both cells!
                    fi = fj = 0.5;
                const Scalar density = (fi*rhoi + fj*rhoj)/(fi + fj);

                // ask for the gravitational acceleration at the given SCV face
                DimVector f(problem.gravityAtPos(face().ipGlobal));

                // calculate the final potential gradient
                f *= density;
                gradPotential_[phaseIdx] -= f;
            } // gravity
        } // loop over all phases
     }

    /*!
     * \brief Actual calculation of the normal Darcy velocities.
     *
     * \param problem The problem
     * \param element The grid element
     * \param elemVolVars The volume variables of the current element
     */
    void calculateNormalVelocity_(const Problem &problem,
                                  const Element &element,
                                  const ElementVolumeVariables &elemVolVars)
    {
        // calculate the mean intrinsic permeability
        const SpatialParams &spatialParams = problem.spatialParams();
        Tensor K;
        if (GET_PROP_VALUE(TypeTag, ImplicitIsBox))
        {
            spatialParams.meanK(K,
                                spatialParams.intrinsicPermeability(element,
                                                                    fvGeometry_(),
                                                                    face().i, fractureIdx),
                                spatialParams.intrinsicPermeability(element,
                                                                    fvGeometry_(),
                                                                    face().j, fractureIdx));
        }
        else
        {
            const Element& elementi = fvGeometry_().neighbors[face().i];
            FVElementGeometry fvGeometryi;
            fvGeometryi.subContVol[0].global = elementi.geometry().center();

            const Element& elementj = fvGeometry_().neighbors[face().j];
            FVElementGeometry fvGeometryj;
            fvGeometryj.subContVol[0].global = elementj.geometry().center();

            spatialParams.meanK(K,
                                spatialParams.intrinsicPermeability(elementi, fvGeometryi, 0,fractureIdx),
                                spatialParams.intrinsicPermeability(elementj, fvGeometryj, 0,fractureIdx));
        }

        // loop over all phases
        for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
        {
            /*
             *  calculate the flux in the normal direction of the
             *  current sub control volume face:
             *  v = - (K_f grad phi) * n
             *  with K_f = rho g / mu K
             *
             *  Mind, that the normal has the length of it's area.
             *  This means that we are actually calculating
             *  Q = - (K grad phi) dot n /|n| * A
             */

            K.mv(gradPotential_[phaseIdx], kGradP_[phaseIdx]);
            kGradPNormal_[phaseIdx] = kGradP_[phaseIdx]*face().normal;

            // determine the upwind direction
            if (kGradPNormal_[phaseIdx] < 0)
            {
                upstreamIdx_[phaseIdx] = face().i;
                downstreamIdx_[phaseIdx] = face().j;
            }
            else
            {
                upstreamIdx_[phaseIdx] = face().j;
                downstreamIdx_[phaseIdx] = face().i;
            }

            // obtain the upwind volume variables
            const VolumeVariables& upVolVars = elemVolVars[ upstreamIdx(phaseIdx) ];
            const VolumeVariables& downVolVars = elemVolVars[ downstreamIdx(phaseIdx) ];

            // the minus comes from the Darcy relation which states that
            // the flux is from high to low potentials.
            // set the velocity
            velocity_[phaseIdx] = kGradP_[phaseIdx];
            velocity_[phaseIdx] *= - ( mobilityUpwindWeight_*upVolVars.mobility(phaseIdx, fractureIdx)
                    + (1.0 - mobilityUpwindWeight_)*downVolVars.mobility(phaseIdx, fractureIdx)) ;

            // set the volume flux
            volumeFlux_[phaseIdx] = velocity_[phaseIdx] * face().normal ;
        }// loop all phases
    }

    // return const reference to fvGeometry
    const FVElementGeometry& fvGeometry_() const
    { return *fvGeometryPtr_; }

    unsigned int faceIdx_;                      //!< The index of the sub control volume face
    bool onBoundary_;                           //!< Specifying whether we are currently on the boundary of the simulation domain
    unsigned int    upstreamIdx_[numPhases] , downstreamIdx_[numPhases]; //!< local index of the upstream / downstream vertex
    Scalar          volumeFlux_[numPhases] ;    //!< Velocity multiplied with normal (magnitude=area)
    DimVector       velocity_[numPhases] ;      //!< The velocity as determined by Darcy's law or by the Forchheimer relation
    Scalar          kGradPNormal_[numPhases] ;  //!< Permeability multiplied with gradient in potential, multiplied with normal (magnitude=area)
    DimVector       kGradP_[numPhases] ;        //!< Permeability multiplied with gradient in potential
    DimVector       gradPotential_[numPhases] ; //!< Gradient of potential, which drives flow
    Scalar          mobilityUpwindWeight_;      //!< Upwind weight for mobility. Set to one for full upstream weighting

private:
    const FVElementGeometry* fvGeometryPtr_; //!< Information about the geometry of discretization
};

} // end namespace

#endif
