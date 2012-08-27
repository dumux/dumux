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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Data which is required to calculate the flux of fluid over a
 *        face of a finite volume
 */
#ifndef DUMUX_RICHARDS_FLUX_VARIABLES_HH
#define DUMUX_RICHARDS_FLUX_VARIABLES_HH

#warning This file is deprecated. Use boxfluxvariables instead. 

#include <dumux/boxmodels/common/boxdarcyfluxvariables.hh>
#include <dumux/common/math.hh>
#include "richardsproperties.hh"

namespace Dumux
{

/*!
 * \ingroup RichardsModel
 * \ingroup BoxFluxVariables
 * \brief Calculates and stores the data which is required to
 *        calculate the flux of fluid over a face of a finite volume.
 */
template <class TypeTag>
class RichardsFluxVariables : public BoxDarcyFluxVariables<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        wPhaseIdx = Indices::wPhaseIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum { dim = GridView::dimension};

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;
public:
    /*!
     * \brief Constructor
     *
     * \param problem The representation of the physical problem
     * \param element The DUNE Codim<0> entity which contains the face of
     *                the finite volume
     * \param fvGeometry The finite volume geometry of the element
     * \param faceIdx The local index of the sub-control volume face in the
     *                element's finite volume geometry.
     * \param elemVolVars An array containing the volume variables for all
     *                    sub-control volumes of the element.
     * \param onBoundary A boolean variable to specify whether the flux variables
     * are calculated for interior SCV faces or boundary faces, default=false
     */
    RichardsFluxVariables(const Problem &problem,
                          const Element &element,
                          const FVElementGeometry &fvGeometry,
                          const int faceIdx,
                          const ElementVolumeVariables &elemVolVars,
                          const bool onBoundary = false)
        : BoxDarcyFluxVariables<TypeTag>(problem, element, fvGeometry, faceIdx, elemVolVars, onBoundary), 
          fvGeometry_(fvGeometry), faceIdx_(faceIdx), onBoundary_(onBoundary)
    {


        calculateGradients_(problem, element, elemVolVars);
        calculateK_(problem, element, elemVolVars);
    };

    /*!
     * \brief Return the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     */
    const DimMatrix &intrinsicPermeability() const
    { return K_; }

    /*!
     * \brief Return the pressure potential gradient \f$\mathrm{[Pa/m]}\f$
     */
    const DimVector &potentialGradW() const
    { return potentialGrad_; }

    /*!
     * \brief Given the intrinisc permeability times the pressure
     *        potential gradient and SCV face normal for a phase,
     *        return the local index of the downstream control volume
     *        for a given phase.
     *
     * \param normalFlux The mass flux over the face multiplied with
     *                   the face's normal.
     */
    int downstreamIdx(Scalar normalFlux) const
    { return (normalFlux >= 0)?face().j:face().i; }

    /*!
     * \brief Given the intrinisc permeability times the pressure
     *        potential gradient and SCV face normal for a phase,
     *        return the local index of the upstream control volume
     *        for a given phase.
     *
     * \param normalFlux The mass flux over the face multiplied with
     *                   the face's normal.
     */
    int upstreamIdx(Scalar normalFlux) const
    { return (normalFlux > 0)?face().i:face().j; }



    /*!
     * \brief The face of the current sub-control volume. This may be either
     *        an inner sub-control-volume face or a face on the boundary.
     */
    const SCVFace &face() const
    {
        if (onBoundary_)
            return fvGeometry_.boundaryFace[faceIdx_];
        else
            return fvGeometry_.subContVolFace[faceIdx_];
    }



protected:
    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemVolVars)
    {
        potentialGrad_ = 0.0;
        // calculate gradients
        for (int idx = 0;
             idx < fvGeometry_.numFAP;
             idx++) // loop over adjacent vertices
        {
            // FE gradient at vertex index
            const DimVector &feGrad = face().grad[idx];

            // index for the element volume variables 
            int volVarsIdx = face().fapIndices[idx];

            // the pressure gradient
            DimVector tmp(feGrad);
            tmp *= elemVolVars[volVarsIdx].pressure(wPhaseIdx);

            potentialGrad_ += tmp;
        }

        ///////////////
        // correct the pressure gradients by the gravitational acceleration
        ///////////////
        if (GET_PARAM(TypeTag, bool, EnableGravity)) {
            // calculate the phase density at the integration point. we
            // only do this if the wetting phase is present in both cells
            Scalar SI = elemVolVars[face().i].saturation(wPhaseIdx);
            Scalar SJ = elemVolVars[face().j].saturation(wPhaseIdx);
            Scalar rhoI = elemVolVars[face().i].density(wPhaseIdx);
            Scalar rhoJ = elemVolVars[face().j].density(wPhaseIdx);
            Scalar fI = std::max(0.0, std::min(SI/1e-5, 0.5));
            Scalar fJ = std::max(0.0, std::min(SJ/1e-5, 0.5));
            if (fI + fJ == 0)
                // doesn't matter because no wetting phase is present in
                // both cells!
                fI = fJ = 0.5;
            Scalar density = (fI*rhoI + fJ*rhoJ)/(fI + fJ);

            // estimate the gravitational acceleration at a given SCV face
            // using the arithmetic mean
            DimVector f(problem.boxGravity(element, fvGeometry_, face().i));
            f += problem.boxGravity(element, fvGeometry_, face().j);
            f /= 2;

            // make it a force
            f *= density;

            // calculate the final potential gradient
            potentialGrad_ -= f;
        }
    }

    void calculateK_(const Problem &problem,
                     const Element &element,
                     const ElementVolumeVariables &elemDat)
    {
        const SpatialParams &spatialParams = problem.spatialParams();
        // calculate the intrinsic permeability
        spatialParams.meanK(K_,
                 spatialParams.intrinsicPermeability(element,
                                          fvGeometry_,
                                          face().i),
                 spatialParams.intrinsicPermeability(element,
                                          fvGeometry_,
                                          face().j));
    }

    const FVElementGeometry &fvGeometry_;
    const int faceIdx_;
    const bool onBoundary_;

    // gradients
    DimVector potentialGrad_;

    // intrinsic permeability
    DimMatrix K_;
};

} // end namepace

#endif
