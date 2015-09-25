// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Onur Dogan                                        *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 * \brief Data which is required to calculate the flux of fluid over a
 *        face of a finite volume
 */
#ifndef DUMUX_RICHARDS_FLUX_VARIABLES_HH
#define DUMUX_RICHARDS_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>

namespace Dumux
{

/*!
 * \ingroup RichardsModel
 * \brief Calculates and stores the data which is required to
 *        calculate the flux of fluid over a face of a finite volume.
 */
template <class TypeTag>
class RichardsFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(RichardsIndices)) Indices;
    enum {
        dimWorld = GridView::dimensionworld,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
    };

    typedef typename FVElementGeometry::SubControlVolume SCV;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

    typedef Dune::FieldVector<Scalar, dimWorld> Vector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;

    typedef typename GridView::template Codim<0>::Entity Element;

public:
    /*!
     * \brief Constructor
     *
     * \param problem The representation of the physical problem
     * \param element The DUNE Codim<0> entity which contains the face of
     *                the finite volume
     * \param fvElemGeom The finite volume geometry of the element
     * \param scvfIdx The local index of the sub-control volume face in the
     *                element's finite volume geometry.
     * \param elemVolVars An array containing the volume variables for all
     *                    sub-control volumes of the element.
     */
    RichardsFluxVariables(const Problem &problem,
                          const Element &element,
                          const FVElementGeometry &fvElemGeom,
                          int scvfIdx,
                          const ElementVolumeVariables &elemVolVars)
        : fvElemGeom_(fvElemGeom)
    {
        scvfIdx_ = scvfIdx;

        calculateGradients_(problem, element, elemVolVars);
        calculateK_(problem, element, elemVolVars);
    };

    /*!
     * \brief Return the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     */
    const Tensor &intrinsicPermeability() const
    { return K_; }

    /*!
     * \brief Return the pressure potential gradient \f$\mathrm{[Pa/m]}\f$
     */
    const Vector &potentialGradW() const
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
     * \brief Returns the face of the element's finite volume geometry
     *        which the flux variables object looks at
     */
    const SCVFace &face() const
    { return fvElemGeom_.subContVolFace[scvfIdx_]; }

protected:
    void calculateGradients_(const Problem &problem,
                             const Element &element,
                             const ElementVolumeVariables &elemVolVars)
    {
        potentialGrad_ = 0.0;
        // calculate gradients
        for (int idx = 0;
             idx < fvElemGeom_.numVertices;
             idx++) // loop over adjacent vertices
        {
            // FE gradient at vertex index
            const Vector &feGrad = face().grad[idx];

            // the pressure gradient
            Vector tmp(feGrad);
            tmp *= elemVolVars[idx].pressure(wPhaseIdx);
            potentialGrad_ += tmp;
        }

        ///////////////
        // correct the pressure gradients by the hydrostatic
        // pressure due to gravity
        ///////////////

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

        Vector tmp(problem.gravity());
        tmp *= density;
        potentialGrad_ -= tmp;
    }

    void calculateK_(const Problem &problem,
                     const Element &element,
                     const ElementVolumeVariables &elemDat)
    {
        const SpatialParameters &sp = problem.spatialParameters();
        // calculate the intrinsic permeability
        sp.meanK(K_,
                 sp.intrinsicPermeability(element,
                                          fvElemGeom_,
                                          face().i),
                 sp.intrinsicPermeability(element,
                                          fvElemGeom_,
                                          face().j));
    }

    const FVElementGeometry &fvElemGeom_;
    int scvfIdx_;

    // gradients
    Vector potentialGrad_;

    // intrinsic permeability
    Tensor K_;
};

} // end namepace

#endif
