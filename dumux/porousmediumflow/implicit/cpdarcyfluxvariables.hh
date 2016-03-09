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
 * \brief This file contains the data which is required to calculate
 *        volume fluxes of fluid phases over a face of a finite volume by means
 *        of the Darcy approximation.
 *
 */
#ifndef DUMUX_CP_DARCY_FLUX_VARIABLES_HH
#define DUMUX_CP_DARCY_FLUX_VARIABLES_HH

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>

#include <dumux/implicit/properties.hh>


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
class CpDarcyFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases)} ;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimWorldMatrix;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> DimVector;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

public:
    /*!
     * \brief The constructor
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param fIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary A boolean variable to specify whether the flux variables
     * are calculated for interior SCV faces or boundary faces, default=false
     */
    CpDarcyFluxVariables(const Problem &problem,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const int fIdx,
                 const ElementVolumeVariables &elemVolVars,
                 const bool onBoundary = false)
    : fvGeometry_(fvGeometry), faceIdx_(fIdx), onBoundary_(onBoundary)
    {
        mobilityUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MobilityUpwindWeight);
        calculateVolumeFlux_(problem, element, elemVolVars);
    }

public:
    /*!
     * \brief Return the volumetric flux over a face of a given phase.
     *
     *        This is the calculated velocity multiplied by the unit normal
     *        and the area of the face. face().normal has already the
     *        magnitude of the area.
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
    GlobalPosition velocity(const unsigned int phaseIdx) const
    { return velocity_[phaseIdx] ; }

    /*!
     * \brief Return the local index of the downstream control volume
     *        for a given phase.
     *
     * \param phaseIdx index of the phase
     */
    const unsigned int downstreamIdx(const unsigned phaseIdx) const
    { return downstreamIdx_[phaseIdx]; }

    /*!
     * \brief Return the local index of the upstream control volume
     *        for a given phase.
     *
     * \param phaseIdx index of the phase
     */
    const unsigned int upstreamIdx(const unsigned phaseIdx) const
    { return upstreamIdx_[phaseIdx]; }

    /*!
     * \brief Return the SCV (sub-control-volume) face. This may be either
     *        a face within the element or a face on the element boundary,
     *        depending on the value of onBoundary_.
     */
    const SCVFace &face() const
    {
        if (onBoundary_)
            return fvGeometry_.boundaryFace[faceIdx_];
        else
            return fvGeometry_.subContVolFace[faceIdx_];
    }

protected:
    /*!
     * \brief Actual calculation of the volume flux.
     *
     * \param problem The problem
     * \param element The finite element
     * \param elemVolVars The volume variables of the current element
     */
    void calculateVolumeFlux_(const Problem &problem,
                              const Element &element,
                              const ElementVolumeVariables &elemVolVars)
    {
        // calculate the transmissibilities
        const SpatialParams &spatialParams = problem.spatialParams();

        const Element& elementI = fvGeometry_.neighbors[face().i];
        FVElementGeometry fvGeometryI;
        fvGeometryI.subContVol[0].global = elementI.geometry().center();
        auto ki = spatialParams.intrinsicPermeability(elementI, fvGeometryI, 0);
        Dune::FieldVector<Scalar, dimWorld> kin;
        ki.mv(face().normal, kin);
        kin /= face().area;
        auto di = face().ipGlobal;
        di -= elementI.geometry().center();
        auto ti = std::abs(di*kin*face().area/(2*di.two_norm2()));

        auto tij = 2*ti;
        if (!onBoundary_)
        {
            const Element& elementJ = fvGeometry_.neighbors[face().j];
            FVElementGeometry fvGeometryJ;
            fvGeometryJ.subContVol[0].global = elementJ.geometry().center();
            auto kj = spatialParams.intrinsicPermeability(elementJ, fvGeometryJ, 0);
            Dune::FieldVector<Scalar, dimWorld> kjn;
            kj.mv(face().normal, kjn);
            kjn /= face().area;
            auto dj = face().ipGlobal;
            dj -= elementJ.geometry().center();
            auto tj = std::abs(dj*kjn*face().area/(2*dj.two_norm2()));
            tij = harmonicMean(ti, tj);
        }

        // loop over all phases
        for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
        {
            const auto& volVarsI = elemVolVars[face().i];
            auto potentialI = volVarsI.pressure(phaseIdx);

            const auto& volVarsJ = elemVolVars[face().j];
            auto potentialJ = volVarsJ.pressure(phaseIdx);

            if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
            {
                // calculate the phase density at the integration point. we
                // only do this if the phase is present in both cells
                Scalar SI = volVarsI.fluidState().saturation(phaseIdx);
                Scalar SJ = volVarsJ.fluidState().saturation(phaseIdx);
                Scalar rhoI = volVarsI.fluidState().density(phaseIdx);
                Scalar rhoJ = volVarsJ.fluidState().density(phaseIdx);
                Scalar fI = std::max(0.0, std::min(SI/1e-5, 0.5));
                Scalar fJ = std::max(0.0, std::min(SJ/1e-5, 0.5));
                if (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(fI + fJ, 0.0, 1.0e-30))
                    // doesn't matter because no wetting phase is present in
                    // both cells!
                    fI = fJ = 0.5;
                const Scalar density = (fI*rhoI + fJ*rhoJ)/(fI + fJ);

                auto globalPosI = elementI.geometry().center();
                potentialI -= density*(problem.gravityAtPos(globalPosI)*globalPosI);

                if (onBoundary_)
                {
                    potentialJ -= density*(problem.gravityAtPos(face().ipGlobal)*face().ipGlobal);
                }
                else
                {
                    const Element& elementJ = fvGeometry_.neighbors[face().j];
                    auto globalPosJ = elementJ.geometry().center();
                    potentialJ -= density*(problem.gravityAtPos(globalPosJ)*globalPosJ);
                }
            }

            auto potentialDiff = potentialI - potentialJ;

            // determine the upwind direction
            if (potentialDiff > 0)
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

            // set the volume flux
            volumeFlux_[phaseIdx] = tij*potentialDiff;
            volumeFlux_[phaseIdx] *= mobilityUpwindWeight_*upVolVars.mobility(phaseIdx)
                    + (1.0 - mobilityUpwindWeight_)*downVolVars.mobility(phaseIdx);

            velocity_[phaseIdx] = face().normal;
            velocity_[phaseIdx] *= volumeFlux_[phaseIdx]/face().area;
        } // over loop all phases
    }

    const FVElementGeometry &fvGeometry_;       //!< Information about the geometry of discretization
    const unsigned int faceIdx_;                //!< The index of the sub control volume face
    const bool      onBoundary_;                //!< Specifying whether we are currently on the boundary of the simulation domain
    unsigned int    upstreamIdx_[numPhases] , downstreamIdx_[numPhases]; //!< local index of the upstream / downstream vertex
    Scalar          volumeFlux_[numPhases] ;    //!< Velocity multiplied with normal (magnitude=area)
    GlobalPosition  velocity_[numPhases] ;      //!< The velocity as determined by Darcy's law or by the Forchheimer relation
    Scalar          mobilityUpwindWeight_;      //!< Upwind weight for mobility. Set to one for full upstream weighting
};

} // end namespace

#endif // DUMUX_CP_DARCY_FLUX_VARIABLES_HH
