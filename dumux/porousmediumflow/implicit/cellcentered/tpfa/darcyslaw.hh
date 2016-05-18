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
 *        volume and mass fluxes of fluid phases over a face of a finite volume by means
 *        of the Darcy approximation.
 */
#ifndef DUMUX_CC_TPFA_DARCYS_LAW_HH
#define DUMUX_CC_TPFA_DARCYS_LAW_HH

#include <memory>

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>

#include <dumux/implicit/properties.hh>


namespace Dumux
{

namespace Properties
{
// forward declaration of properties
NEW_PROP_TAG(NumPhases);
NEW_PROP_TAG(ProblemEnableGravity);
}

/*!
 * \ingroup CCTpfaDarcysLaw
 * \brief Evaluates the normal component of the Darcy velocity
 * on a (sub)control volume face.
 */
template <class TypeTag>
class CCTpfaDarcysLaw
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace) SubControlVolumeFace;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GridView::IndexSet::IndexType IndexType;
    typedef std::vector<Scalar> TransmissivityVector;
    typedef std::vector<IndexType> Stencil;

    using Element = typename GridView::template Codim<0>::Entity;
    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases)} ;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimWorldMatrix;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:

    static Scalar flux(const Problem& problem,
                       const SubControlVolumeFace& scvFace,
                       const IndexType phaseIdx,
                       std::shared_ptr<VolumeVariables> boundaryVolVars)
    {
        const auto& tij = getTransmissibilities(problem, scvFace);

        // Get the inside volume variables
        const auto insideScvIdx = scvFace.insideScvIdx();
        const auto& insideScv = problem.model().fvGeometries().subControlVolume(insideScvIdx);
        const auto* insideVolVars = &problem.model().curVolVars(insideScv);

        // and the outside volume variables
        const VolumeVariables* outsideVolVars;
        if (!scvFace.boundary())
            outsideVolVars = &problem.model().curVolVars(scvFace.outsideScvIdx());
        else
        {
            if (!boundaryVolVars)
                DUNE_THROW(Dune::InvalidStateException, "Trying to access invalid boundary volume variables.");
            outsideVolVars = boundaryVolVars.get();
        }

        auto hInside = insideVolVars->pressure(phaseIdx);
        auto hOutside = outsideVolVars->pressure(phaseIdx);

        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
        {
            // do averaging for the density
            const auto rhoInside = insideVolVars->density(phaseIdx);
            const auto rhoOutide = outsideVolVars->density(phaseIdx);
            const auto rho = (rhoInside + rhoOutide)*0.5;


            // ask for the gravitational acceleration in the inside neighbor
            const auto xInside = insideScv.center();
            const auto gInside = problem.gravityAtPos(xInside);

            hInside -= rho*(gInside*xInside);

            // and the outside neighbor
            if (scvFace.boundary())
            {
                const auto xOutside = scvFace.center();
                const auto gOutside = problem.gravityAtPos(xOutside);
                hOutside -= rho*(gOutside*xOutside);
            }
            else
            {
                const auto outsideScvIdx = scvFace.outsideScvIdx();
                const auto& outsideScv = problem.model().fvGeometries().subControlVolume(outsideScvIdx);
                const auto xOutside = outsideScv.center();
                const auto gOutside = problem.gravityAtPos(xOutside);
                hOutside -= rho*(gOutside*xOutside);
            }
        }

        return tij[0]*(hInside - hOutside);
    }

    static Stencil stencil(const Problem& problem, const SubControlVolumeFace& scvFace)
    {
        std::vector<IndexType> stencil;
        stencil.clear();
        if (!scvFace.boundary())
        {
            stencil.push_back(scvFace.insideScvIdx());
            stencil.push_back(scvFace.outsideScvIdx());
        }
        else
            stencil.push_back(scvFace.insideScvIdx());

        return stencil;
    }

    template <typename T = TypeTag>
    static const typename std::enable_if<GET_PROP_VALUE(T, EnableFluxVariablesCache), TransmissivityVector>::type& getTransmissibilities(const Problem& problem, const SubControlVolumeFace& scvFace)
    {
        return problem.model().fluxVarsCache(scvFace).tij();
    }

    template <typename T = TypeTag>
    static const typename std::enable_if<!GET_PROP_VALUE(T, EnableFluxVariablesCache), TransmissivityVector>::type getTransmissibilities(const Problem& problem, const SubControlVolumeFace& scvFace)
    {
        return calculateTransmissibilities(problem, scvFace);
    }

    static TransmissivityVector calculateTransmissibilities(const Problem& problem, const SubControlVolumeFace& scvFace)
    {
        Scalar tij;

        const auto insideScvIdx = scvFace.insideScvIdx();
        const auto& insideScv = problem.model().fvGeometries().subControlVolume(insideScvIdx);
        const auto insideK = problem.spatialParams().intrinsicPermeability(insideScv);
        Scalar ti = calculateOmega_(problem, scvFace, insideK, insideScv);

        if (!scvFace.boundary())
        {
            const auto outsideScvIdx = scvFace.outsideScvIdx();
            const auto& outsideScv = problem.model().fvGeometries().subControlVolume(outsideScvIdx);
            const auto outsideK = problem.spatialParams().intrinsicPermeability(outsideScv);
            Scalar tj = -1.0*calculateOmega_(problem, scvFace, outsideK, outsideScv);

            tij = scvFace.area()*(ti * tj)/(ti + tj);
        }
        else
        {
            tij = scvFace.area()*ti;
        }

        return TransmissivityVector(1, tij);
    }

private:

    static Scalar calculateOmega_(const Problem& problem, const SubControlVolumeFace& scvFace, const DimWorldMatrix &K, const SubControlVolume &scv)
    {
        GlobalPosition Knormal;
        K.mv(scvFace.unitOuterNormal(), Knormal);

        auto distanceVector = scvFace.center();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = Knormal * distanceVector;
        omega *= problem.model().curVolVars(scv).extrusionFactor();

        return omega;
    }

    static Scalar calculateOmega_(const Problem& problem, const SubControlVolumeFace& scvFace, Scalar K, const SubControlVolume &scv)
    {
        auto distanceVector = scvFace.center();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = K * (distanceVector * scvFace.unitOuterNormal());
        omega *= problem.model().curVolVars(scv).extrusionFactor();

        return omega;
    }
};

} // end namespace

#endif
