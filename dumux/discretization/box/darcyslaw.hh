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
#ifndef DUMUX_DISCRETIZATION_BOX_DARCYS_LAW_HH
#define DUMUX_DISCRETIZATION_BOX_DARCYS_LAW_HH

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
 * \brief Specialization of Darcy's Law for the box method.
 */
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
        for (const auto& scv : scvs(fvGeometry))
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
