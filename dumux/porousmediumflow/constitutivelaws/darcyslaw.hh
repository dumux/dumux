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
 *        of the Darcy approximation. Specializations are provided for the different discretization methods.
 */
#ifndef DUMUX_POROUSMEDIUMFLOW_DARCYS_LAW_HH
#define DUMUX_POROUSMEDIUMFLOW_DARCYS_LAW_HH

#include <memory>

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>

#include <dumux/implicit/properties.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>


namespace Dumux
{

namespace Properties
{
// forward declaration of properties
NEW_PROP_TAG(ProblemEnableGravity);
}

/*!
 * \ingroup DarcysLaw
 * \brief Evaluates the normal component of the Darcy velocity
 * on a (sub)control volume face. Specializations are provided
 * for the different discretization methods.
 */
template <class TypeTag, typename DiscretizationMethod = void>
class DarcysLaw
{};

// Specialization for the CC-Tpfa method
template <class TypeTag>
class DarcysLaw<TypeTag, typename std::enable_if<GET_PROP_VALUE(TypeTag, DiscretizationMethod) == GET_PROP(TypeTag, DiscretizationMethods)::CCTpfa>::type >
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace) SubControlVolumeFace;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::IndexSet::IndexType IndexType;
    typedef std::vector<IndexType> Stencil;

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimWorldMatrix;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:

    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const SubControlVolumeFace& scvFace,
                       const IndexType phaseIdx)
    {
        const auto& tij = problem.model().fluxVarsCache(scvFace).tij();

        // Get the inside and outside volume variables
        const auto& insideScv = problem.model().fvGeometries().subControlVolume(scvFace.insideScvIdx());
        const auto& insideVolVars = problem.model().curVolVars(insideScv);
        const auto& outsideVolVars = problem.model().curVolVars(scvFace.outsideScvIdx());

        auto hInside = insideVolVars.pressure(phaseIdx);
        auto hOutside = outsideVolVars.pressure(phaseIdx);

        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
        {
            // do averaging for the density
            const auto rhoInside = insideVolVars.density(phaseIdx);
            const auto rhoOutide = outsideVolVars.density(phaseIdx);
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

        return tij*(hInside - hOutside);
    }

    static Stencil stencil(const Problem& problem, const Element& element, const SubControlVolumeFace& scvFace)
    {
        Stencil stencil;
        if (!scvFace.boundary())
        {
            stencil.push_back(scvFace.insideScvIdx());
            stencil.push_back(scvFace.outsideScvIdx());
        }
        else
            stencil.push_back(scvFace.insideScvIdx());

        return stencil;
    }

    // The flux variables cache has to be bound to an element prior to flux calculations
    // During the binding, the transmissibilities will be computed and stored using the method below.
    static Scalar calculateTransmissibilities(const Problem& problem, const Element& element, const SubControlVolumeFace& scvFace)
    {
        Scalar tij;

        const auto insideScvIdx = scvFace.insideScvIdx();
        const auto& insideScv = problem.model().fvGeometries().subControlVolume(insideScvIdx);
        const auto insideK = problem.spatialParams().intrinsicPermeability(insideScv);
        Scalar ti = calculateOmega_(problem, scvFace, insideK, element, insideScv);

        if (!scvFace.boundary())
        {
            const auto outsideScvIdx = scvFace.outsideScvIdx();
            const auto& outsideScv = problem.model().fvGeometries().subControlVolume(outsideScvIdx);
            const auto outsideElement = problem.model().fvGeometries().element(outsideScvIdx);
            const auto outsideK = problem.spatialParams().intrinsicPermeability(outsideScv);
            Scalar tj = -1.0*calculateOmega_(problem, scvFace, outsideK, outsideElement, outsideScv);

            tij = scvFace.area()*(ti * tj)/(ti + tj);
        }
        else
        {
            tij = scvFace.area()*ti;
        }

        return tij;
    }

private:

    static Scalar calculateOmega_(const Problem& problem, const SubControlVolumeFace& scvFace, const DimWorldMatrix &K, const Element& element, const SubControlVolume &scv)
    {
        GlobalPosition Knormal;
        K.mv(scvFace.unitOuterNormal(), Knormal);

        auto distanceVector = scvFace.center();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = Knormal * distanceVector;
        omega *= problem.boxExtrusionFactor(element, scv);

        return omega;
    }

    static Scalar calculateOmega_(const Problem& problem, const SubControlVolumeFace& scvFace, const Scalar K, const Element& element, const SubControlVolume &scv)
    {
        auto distanceVector = scvFace.center();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = K * (distanceVector * scvFace.unitOuterNormal());
        omega *= problem.boxExtrusionFactor(element, scv);

        return omega;
    }
};

// Specialization for the Box Method
template <class TypeTag>
class DarcysLaw<TypeTag, typename std::enable_if<GET_PROP_VALUE(TypeTag, DiscretizationMethod) == GET_PROP(TypeTag, DiscretizationMethods)::Box>::type >
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace) SubControlVolumeFace;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariablesCache) FluxVariablesCache;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::IndexSet::IndexType IndexType;
    typedef typename GridView::ctype CoordScalar;
    typedef std::vector<IndexType> Stencil;

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimWorldMatrix;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

    typedef Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1> FeCache;
    typedef typename FeCache::FiniteElementType::Traits::LocalBasisType FeLocalBasis;
    typedef typename FluxVariablesCache::FaceData FaceData;

public:

    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const SubControlVolumeFace& scvFace,
                       const IndexType phaseIdx)
    {
        // get the precalculated local jacobian and shape values at the integration point
        const auto& faceData = problem.model().fluxVarsCache()[scvFace].faceData();

        const auto& fvGeometry = problem.model().fvGeometries(element);
        const auto& insideScv = problem.model().fvGeometries().subControlVolume(scvFace.insideScvIdx());
        const auto extrusionFactor = problem.model().curVolVars(insideScv).extrusionFactor();
        const auto K = problem.spatialParams().intrinsicPermeability(insideScv);

        // evaluate gradP - rho*g at integration point
        DimVector gradP(0.0);
        for (const auto& scv : fvGeometry.scvs())
        {
            // the global shape function gradient
            DimVector gradI;
            faceData.jacInvT.mv(faceData.localJacobian[scv.indexInElement()][0], gradI);

            gradI *= problem.model().curVolVars(scv).pressure(phaseIdx);
            gradP += gradI;
        }
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
        {
            // gravitational acceleration
            DimVector g(problem.gravityAtPos(scvFace.center()));

            // interpolate the density at the IP
            const auto& insideVolVars = problem.model().curVolVars(insideScv);
            const auto& outsideScv = problem.model().fvGeometries().subControlVolume(scvFace.outsideScvIdx());
            const auto& outsideVolVars = problem.model().curVolVars(outsideScv);
            Scalar rho = 0.5*(insideVolVars.density(phaseIdx) + outsideVolVars.density(phaseIdx));

            // turn gravity into a force
            g *= rho;

            // subtract from pressure gradient
            gradP -= g;
        }

        // apply the permeability and return the flux
        auto KGradP = applyPermeability(K, gradP);
        return -1.0*(KGradP*scvFace.unitOuterNormal())*scvFace.area()*extrusionFactor;
    }

    static DimVector applyPermeability(const DimWorldMatrix& K, const DimVector& gradI)
    {
        DimVector result(0.0);
        K.mv(gradI, result);

        return result;
    }

    static DimVector applyPermeability(const Scalar k, const DimVector& gradI)
    {
        DimVector result(gradI);
        result *= k;
        return result;
    }

    // This is for compatibility with the cc methods. The flux stencil info is obsolete for the box method.
    static Stencil stencil(const Problem& problem, const Element& element, const SubControlVolumeFace& scvFace)
    {
        return Stencil(0);
    }

    static FaceData calculateFaceData(const Problem& problem, const Element& element, const typename Element::Geometry& geometry, const FeLocalBasis& localBasis, const SubControlVolumeFace& scvFace)
    {
        FaceData faceData;

        // evaluate shape functions and gradients at the integration point
        const auto ipLocal = geometry.local(scvFace.center());
        faceData.jacInvT = geometry.jacobianInverseTransposed(ipLocal);
        localBasis.evaluateJacobian(ipLocal, faceData.localJacobian);
        localBasis.evaluateFunction(ipLocal, faceData.shapeValues);

        return std::move(faceData);
    }
};

} // end namespace

#endif
