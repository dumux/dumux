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
#ifndef DUMUX_IMPLICIT_DARCY_FLUX_VARIABLES_HH
#define DUMUX_IMPLICIT_DARCY_FLUX_VARIABLES_HH

#include <dune/common/float_cmp.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/deprecated.hh>

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
class ImplicitDarcyFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) Implementation;
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
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimWorldMatrix;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

public:
    /*!
     * \brief The old constructor
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param fIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary A boolean variable to specify whether the flux variables
     * are calculated for interior SCV faces or boundary faces, default=false
     */
    DUNE_DEPRECATED_MSG("FluxVariables now have to be default constructed and updated.")
    ImplicitDarcyFluxVariables(const Problem &problem,
                               const Element &element,
                               const FVElementGeometry &fvGeometry,
                               const int fIdx,
                               const ElementVolumeVariables &elemVolVars,
                               const bool onBoundary = false)
    {
        DUNE_THROW(Dune::InvalidStateException, "The FluxVariables now have to be default contructed. "
                                                << "In case you have your own FluxVariables you have to make them default "
                                                << " constructable too. All calls to the old constructor will throw this error. "
                                                << "Everywhere you instantiate FluxVariables do this now by default constructing "
                                                << "a FluxVariables object (FluxVariables fluxVars;) and then updating it where "
                                                << "the update method has the same signature as the old constructor (fluxVars.update(...).)");
    }

    /*!
     * \brief Default constructor
     * \note This can be removed when the deprecated constructor is removed.
     */
    ImplicitDarcyFluxVariables() = default;

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
        asImp_().calculateGradients_(problem, element, elemVolVars);
        asImp_().calculateNormalVelocity_(problem, element, elemVolVars);
    }


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
            return fvGeometry_().boundaryFace[faceIdx_];
        else
            return fvGeometry_().subContVolFace[faceIdx_];
    }

protected:

    //! Returns the implementation of the flux variables (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    /*!
     * \brief Calculation of the potential gradients
     *
     * \param problem The problem
     * \param element The finite element
     * \param elemVolVars The volume variables of the current element
     */
    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemVolVars)
    {
        // loop over all phases
        for (int phaseIdx = 0; phaseIdx < numPhases; phaseIdx++)
        {
            potentialGrad_[phaseIdx]= 0.0;

            for (unsigned int idx = 0;
                 idx < face().numFap;
                 idx++) // loop over adjacent vertices
            {
                // FE gradient at vertex idx
                const GlobalPosition &feGrad = face().grad[idx];

                // index for the element volume variables
                int volVarsIdx = face().fapIndices[idx];

                // the pressure gradient
                GlobalPosition tmp(feGrad);
                tmp *= elemVolVars[volVarsIdx].fluidState().pressure(phaseIdx);
                potentialGrad_[phaseIdx] += tmp;
            }

            // correct the pressure gradient by the gravitational acceleration
            if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
            {
                // ask for the gravitational acceleration at the given SCV face
                GlobalPosition g(problem.gravityAtPos(face().ipGlobal));

                // calculate the phase density at the integration point. we
                // only do this if the wetting phase is present in both cells
                Scalar SI = elemVolVars[face().i].fluidState().saturation(phaseIdx);
                Scalar SJ = elemVolVars[face().j].fluidState().saturation(phaseIdx);
                Scalar rhoI = elemVolVars[face().i].fluidState().density(phaseIdx);
                Scalar rhoJ = elemVolVars[face().j].fluidState().density(phaseIdx);
                Scalar fI = std::max(0.0, std::min(SI/1e-5, 0.5));
                Scalar fJ = std::max(0.0, std::min(SJ/1e-5, 0.5));
                if (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(fI + fJ, 0.0, 1.0e-30))
                    // doesn't matter because no wetting phase is present in
                    // both cells!
                    fI = fJ = 0.5;
                const Scalar density = (fI*rhoI + fJ*rhoJ)/(fI + fJ);

                // make gravity acceleration a force
                GlobalPosition f(g);
                f *= density;

                // calculate the final potential gradient
                potentialGrad_[phaseIdx] -= f;
            } // gravity
        } // loop over all phases
     }

    /*!
     * \brief Actual calculation of the normal Darcy velocities.
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
        const SpatialParams &spatialParams = problem.spatialParams();
        DimWorldMatrix K;
        if (GET_PROP_VALUE(TypeTag, ImplicitIsBox))
        {
            spatialParams.meanK(K,
                                spatialParams.intrinsicPermeability(element,
                                                                    fvGeometry_(),
                                                                    face().i),
                                spatialParams.intrinsicPermeability(element,
                                                                    fvGeometry_(),
                                                                    face().j));
        }
        else
        {
            const Element& elementI = fvGeometry_().neighbors[face().i];
            FVElementGeometry fvGeometryI;
            fvGeometryI.subContVol[0].global = elementI.geometry().center();

            const Element& elementJ = fvGeometry_().neighbors[face().j];
            FVElementGeometry fvGeometryJ;
            fvGeometryJ.subContVol[0].global = elementJ.geometry().center();

            spatialParams.meanK(K,
                                spatialParams.intrinsicPermeability(elementI, fvGeometryI, 0),
                                spatialParams.intrinsicPermeability(elementJ, fvGeometryJ, 0));
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


            K.mv(potentialGrad_[phaseIdx], kGradP_[phaseIdx]);
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
            velocity_[phaseIdx] *= - ( mobilityUpwindWeight_*upVolVars.mobility(phaseIdx)
                    + (1.0 - mobilityUpwindWeight_)*downVolVars.mobility(phaseIdx)) ;

            // set the volume flux
            volumeFlux_[phaseIdx] = velocity_[phaseIdx] * face().normal;
        }// loop all phases
    }

    // return const reference to the fvGeometry
    const FVElementGeometry& fvGeometry_() const
    { return *fvGeometryPtr_; }

    unsigned int faceIdx_;                //!< The index of the sub control volume face
    bool onBoundary_;                //!< Specifying whether we are currently on the boundary of the simulation domain
    unsigned int    upstreamIdx_[numPhases] , downstreamIdx_[numPhases]; //!< local index of the upstream / downstream vertex
    Scalar          volumeFlux_[numPhases] ;    //!< Velocity multiplied with normal (magnitude=area)
    GlobalPosition  velocity_[numPhases] ;      //!< The velocity as determined by Darcy's law or by the Forchheimer relation
    Scalar          kGradPNormal_[numPhases] ;  //!< Permeability multiplied with gradient in potential, multiplied with normal (magnitude=area)
    GlobalPosition  kGradP_[numPhases] ; //!< Permeability multiplied with gradient in potential
    GlobalPosition  potentialGrad_[numPhases] ; //!< Gradient of potential, which drives flow
    Scalar          mobilityUpwindWeight_;      //!< Upwind weight for mobility. Set to one for full upstream weighting

private:
    const FVElementGeometry* fvGeometryPtr_; //!< Information about the geometry of discretization

};

} // end namespace

#endif // DUMUX_IMPLICIT_DARCY_FLUX_VARIABLES_HH
