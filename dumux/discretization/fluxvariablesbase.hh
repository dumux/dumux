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
 * \brief Base class for the flux variables
 */
#ifndef DUMUX_DISCRETIZATION_FLUXVARIABLESBASE_HH
#define DUMUX_DISCRETIZATION_FLUXVARIABLESBASE_HH

#include <dumux/implicit/properties.hh>

namespace Dumux
{

namespace Properties
{
// forward declaration
NEW_PROP_TAG(ImplicitMassUpwindWeight);
}

/*!
 * \ingroup Discretization
 * \brief Base class for the flux variables
 *        Actual flux variables inherit from this class
 */
template<class TypeTag, class Implementation>
class FluxVariablesBase
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Stencil = std::vector<IndexType>;
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);

public:

    void init(const Problem& problem,
              const Element& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars,
              const SubControlVolumeFace &scvFace,
              const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        problemPtr_ = &problem;
        elementPtr_ = &element;
        scvFacePtr_ = &scvFace;
        fvGeometryPtr_ = &fvGeometry;
        elemVolVarsPtr_ = &elemVolVars;
        elemFluxVarsCachePtr_ = &elemFluxVarsCache;
        // retrieve the upwind weight for the mass conservation equations. Use the value
        // specified via the property system as default, and overwrite
        // it by the run-time parameter from the Dune::ParameterTree
        upwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
    }

    const Problem& problem() const
    { return *problemPtr_; }

    const Element& element() const
    { return *elementPtr_; }

    const SubControlVolumeFace& scvFace() const
    { return *scvFacePtr_; }

    const FVElementGeometry& fvGeometry() const
    { return *fvGeometryPtr_; }

    const ElementVolumeVariables& elemVolVars() const
    { return *elemVolVarsPtr_; }

    const ElementFluxVariablesCache& elemFluxVarsCache() const
    { return *elemFluxVarsCachePtr_; }

    Stencil computeStencil(const Problem& problem, const Element& element, const SubControlVolumeFace& scvFace)
    { DUNE_THROW(Dune::InvalidStateException, "computeStencil() routine is not provided by the implementation."); }

    // For cell-centered surface and network grids (dim < dimWorld) we have to do a special upwind scheme
    template<typename FunctionType, class T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(T, ImplicitIsBox) &&
                             GET_PROP_TYPE(T, Grid)::dimension < GET_PROP_TYPE(T, Grid)::dimensionworld, Scalar>::type
    upwindScheme(Scalar flux, int phaseIdx, const FunctionType& upwindTerm)
    {
        const auto& insideScv = this->fvGeometry().scv(this->scvFace().insideScvIdx());
        const auto& insideVolVars = this->elemVolVars()[insideScv];

        // check if this is a branching point
        if (this->scvFace().numOutsideScvs() > 1)
        {
            // more complicated upwind scheme
            // we compute a flux-weighted average of all inflowing branches
            Scalar branchingPointUpwindTerm = 0.0;
            Scalar sumUpwindFluxes = 0.0;

            // the inside flux
            if (!std::signbit(flux))
                branchingPointUpwindTerm += upwindTerm(insideVolVars)*flux;
            else
                sumUpwindFluxes += flux;

            for (unsigned int i = 0; i < this->scvFace().numOutsideScvs(); ++i)
            {
                 // compute the outside flux
                const auto outsideScvIdx = this->scvFace().outsideScvIdx(i);
                const auto outsideElement = this->fvGeometry().globalFvGeometry().element(outsideScvIdx);
                const auto& flippedScvf = this->fvGeometry().flipScvf(this->scvFace().index(), i);

                const auto outsideFlux = AdvectionType::flux(this->problem(),
                                                             outsideElement,
                                                             this->fvGeometry(),
                                                             this->elemVolVars(),
                                                             flippedScvf,
                                                             phaseIdx,
                                                             this->elemFluxVarsCache());

                if (!std::signbit(outsideFlux))
                    branchingPointUpwindTerm += upwindTerm(this->elemVolVars()[outsideScvIdx])*outsideFlux;
                else
                    sumUpwindFluxes += outsideFlux;
            }

            // the flux might be zero
            if (sumUpwindFluxes != 0.0)
                branchingPointUpwindTerm /= -sumUpwindFluxes;
            else
                branchingPointUpwindTerm = 0.0;

            // upwind scheme (always do fully upwind at branching points)
            // a weighting here would lead to an error since the derivation is based on a fully upwind scheme
            // TODO How to implement a weight of e.g. 0.5
            if (std::signbit(flux))
                return flux*branchingPointUpwindTerm;
            else
                return flux*upwindTerm(insideVolVars);
        }
        // non-branching points and boundaries
        else
        {
            // upwind scheme
            const auto& outsideVolVars = this->elemVolVars()[this->scvFace().outsideScvIdx()];
            if (std::signbit(flux))
                return flux*(upwindWeight_*upwindTerm(outsideVolVars)
                             + (1.0 - upwindWeight_)*upwindTerm(insideVolVars));
            else
                return flux*(upwindWeight_*upwindTerm(insideVolVars)
                             + (1.0 - upwindWeight_)*upwindTerm(outsideVolVars));
        }
    }

    // For grids with dim == dimWorld or the box-method we use a simple upwinding scheme
    template<typename FunctionType, class T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, ImplicitIsBox) ||
                            GET_PROP_TYPE(T, Grid)::dimension == GET_PROP_TYPE(T, Grid)::dimensionworld, Scalar>::type
    upwindScheme(Scalar flux, int phaseIdx, const FunctionType& upwindTerm)
    {
        const auto& insideScv = this->fvGeometry().scv(this->scvFace().insideScvIdx());
        const auto& insideVolVars = this->elemVolVars()[insideScv];

        // upwind scheme
        const auto& outsideVolVars = this->elemVolVars()[this->scvFace().outsideScvIdx()];
        if (std::signbit(flux))
            return flux*(upwindWeight_*upwindTerm(outsideVolVars)
                         + (1.0 - upwindWeight_)*upwindTerm(insideVolVars));
        else
            return flux*(upwindWeight_*upwindTerm(insideVolVars)
                         + (1.0 - upwindWeight_)*upwindTerm(outsideVolVars));
    }

private:

    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }

    const Problem* problemPtr_;              //! Pointer to the problem
    const Element* elementPtr_;              //! Pointer to the element at hand
    const FVElementGeometry* fvGeometryPtr_;
    const SubControlVolumeFace* scvFacePtr_; //! Pointer to the sub control volume face for which the flux variables are created
    const ElementVolumeVariables* elemVolVarsPtr_;
    const ElementFluxVariablesCache* elemFluxVarsCachePtr_;
    Scalar upwindWeight_;
};

} // end namespace

#endif
