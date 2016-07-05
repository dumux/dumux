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
#ifndef DUMUX_POROUSMEDIUM_IMPLICIT_FLUXVARIABLESCACHE_HH
#define DUMUX_POROUSMEDIUM_IMPLICIT_FLUXVARIABLESCACHE_HH

#include <dumux/implicit/properties.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief The flux variables cache classes for porous media.
 *        Store flux stencils and data required for flux calculation
 */
template<class TypeTag, typename DiscretizationMethod = void>
class PorousMediumFluxVariablesCache {};

// specialization for the Box Method
template<class TypeTag>
class PorousMediumFluxVariablesCache<TypeTag, typename std::enable_if<GET_PROP_VALUE(TypeTag, DiscretizationMethod) == GET_PROP(TypeTag, DiscretizationMethods)::Box>::type >
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;
    using TransmissibilityVector = std::vector<IndexType>;

    typedef typename GridView::ctype CoordScalar;
    static const int dim = GridView::dimension;
    typedef Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1> FeCache;
    typedef typename FeCache::FiniteElementType::Traits::LocalBasisType FeLocalBasis;
    typedef typename FeLocalBasis::Traits::JacobianType ShapeJacobian;
    typedef typename Dune::FieldVector<Scalar, 1> ShapeValue;
    typedef typename Element::Geometry::JacobianInverseTransposed JacobianInverseTransposed;

public:
    // stores required data for the flux calculation
    struct FaceData
    {
        std::vector<ShapeJacobian> localJacobian;
        std::vector<ShapeValue> shapeValues;
        JacobianInverseTransposed jacInvT;
    };

    void update(const Problem& problem,
                const Element& element,
                const typename Element::Geometry& geometry,
                const FeLocalBasis& localBasis,
                const SubControlVolumeFace &scvFace)
    {
        faceData_ = AdvectionType::calculateFaceData(problem, element, geometry, localBasis, scvFace);

        // The stencil info is obsolete for the box method.
        // It is here for compatibility with cc methods
        stencil_ = Stencil(0);
    }

    const FaceData& faceData() const
    { return faceData_; }

    const Stencil& stencil() const
    {
        return stencil_;
    }

private:
    FaceData faceData_;
    Stencil stencil_;
};

// specialization for the cell centered tpfa method
template<class TypeTag>
class PorousMediumFluxVariablesCache<TypeTag, typename std::enable_if<GET_PROP_VALUE(TypeTag, DiscretizationMethod) == GET_PROP(TypeTag, DiscretizationMethods)::CCTpfa>::type >
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;

public:
    void update(const Problem& problem,
                const Element& element,
                const SubControlVolumeFace &scvFace)
    {
        FluxVariables fluxVars;
        stencil_ = fluxVars.computeStencil(problem, element, scvFace);
        tij_ = AdvectionType::calculateTransmissibilities(problem, element, scvFace);
    }

    const Stencil& stencil() const
    { return stencil_; }

    const Scalar& tij() const
    { return tij_; }

private:
    Stencil stencil_;
    Scalar tij_;
};

} // end namespace

#endif
