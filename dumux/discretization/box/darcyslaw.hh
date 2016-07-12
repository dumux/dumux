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
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using CoordScalar = typename GridView::ctype;
    using Stencil = std::vector<IndexType>;

    enum { dim = GridView::dimension};
    enum { dimWorld = GridView::dimensionworld};

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using FaceData = typename FluxVariablesCache::FaceData;

public:

    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const SubControlVolumeFace& scvf,
                       const IndexType phaseIdx)
    {
        // get the precalculated local jacobian and shape values at the integration point
        const auto& faceData = problem.model().fluxVarsCache()[scvf].faceData();

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto extrusionFactor = problem.model().curVolVars(insideScv).extrusionFactor();
        const auto K = problem.spatialParams().intrinsicPermeability(insideScv);

        // evaluate gradP - rho*g at integration point
        DimVector gradP(0.0);
        for (auto&& scv : scvs(fvGeometry))
        {
            // the global shape function gradient
            DimVector gradI;
            faceData.jacInvT.mv(faceData.shapeJacobian[scv.index()][0], gradI);

            gradI *= problem.model().curVolVars(scv).pressure(phaseIdx);
            gradP += gradI;
        }
        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
        {
            // gravitational acceleration
            DimVector g(problem.gravityAtPos(scvf.center()));

            // interpolate the density at the IP
            const auto& insideVolVars = problem.model().curVolVars(insideScv);
            const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
            const auto& outsideVolVars = problem.model().curVolVars(outsideScv);
            Scalar rho = 0.5*(insideVolVars.density(phaseIdx) + outsideVolVars.density(phaseIdx));

            // turn gravity into a force
            g *= rho;

            // subtract from pressure gradient
            gradP -= g;
        }

        // apply the permeability and return the flux
        auto KGradP = applyPermeability(K, gradP);
        return -1.0*(KGradP*scvf.unitOuterNormal())*scvf.area()*extrusionFactor;
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
    static Stencil stencil(const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const SubControlVolumeFace& scvf)
    {
        return Stencil(0);
    }

    static FaceData calculateFaceData(const Problem& problem,
                                      const Element& element,
                                      const FVElementGeometry& fvGeometry,
                                      const SubControlVolumeFace& scvf)
    {
        FaceData faceData;
        const auto geometry = element.geometry();
        const auto& localBasis = fvGeometry.feLocalBasis();

        // evaluate shape functions and gradients at the integration point
        const auto ipLocal = geometry.local(scvf.center());
        faceData.jacInvT = geometry.jacobianInverseTransposed(ipLocal);
        localBasis.evaluateJacobian(ipLocal, faceData.shapeJacobian);
        localBasis.evaluateFunction(ipLocal, faceData.shapeValue);

        return faceData;
    }
};

} // end namespace

#endif
