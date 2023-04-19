// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Geometry
 * \author Timo Koch
 * \brief A quadrature based on refinement
 */

#ifndef DUMUX_GEOMETRY_REFINEMENT_QUADRATURERULE_HH
#define DUMUX_GEOMETRY_REFINEMENT_QUADRATURERULE_HH

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/virtualrefinement.hh>
#include <dumux/common/math.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief A "quadrature" based on virtual refinement
 */
template<typename ct, int mydim>
class RefinementQuadratureRule final : public Dune::QuadratureRule<ct, mydim>
{
public:
    //! The space dimension
    static constexpr int dim = mydim;

    //! The highest quadrature order available
    static constexpr int highest_order = 1;

    ~RefinementQuadratureRule() final {}

    RefinementQuadratureRule(Dune::GeometryType type, int levels)
    : Dune::QuadratureRule<ct, dim>(type, highest_order)
    {
        auto& refinement = Dune::buildRefinement<dim, ct>(type, type);
        const auto tag = Dune::refinementLevels(levels);

        // get the vertices
        const auto numVertices = refinement.nVertices(tag);
        std::vector<Dune::FieldVector<ct, dim>> points; points.reserve(numVertices);
        auto vIt = refinement.vBegin(tag);
        const auto vItEnd = refinement.vEnd(tag);
        for (; vIt != vItEnd; ++vIt)
            points.emplace_back(vIt.coords());

        // go over all elements and get center and volume
        const auto numElements = refinement.nElements(tag);
        this->reserve(numElements);
        auto eIt = refinement.eBegin(tag);
        const auto eItEnd = refinement.eEnd(tag);
        for (; eIt != eItEnd; ++eIt)
        {
            // quadrature point in the centroid of the sub-element and weight is its volume
            const auto weight = computeVolume_(type, points, eIt.vertexIndices());
            this->emplace_back(eIt.coords(), weight);
        }
    }

private:
    template<class VertexIndices>
    ct computeVolume_(Dune::GeometryType t, const std::vector<Dune::FieldVector<ct, dim>>& points, const VertexIndices& indices) const
    {
        if (t != Dune::GeometryTypes::simplex(2))
            DUNE_THROW(Dune::NotImplemented, "Only implemented for 2d simplices");

        const auto ab = points[indices[1]] - points[indices[0]];
        const auto ac = points[indices[2]] - points[indices[0]];
        using std::abs;
        return abs(crossProduct(ab, ac))*0.5;
    }
};

} // end namespace Dumux

#endif
