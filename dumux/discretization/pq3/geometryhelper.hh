// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file includes code adapted from dune-functions.
// These parts are clearly marked with inline comments and were originally
// licensed under LGPL-3.0-or-later. They are adapted and relicensed here
// under the terms of the GNU GPL-3.0-or-later as permitted by LGPLv3 Sec. 3.
//
// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
//
/*!
 * \file
 * \ingroup PQ3Discretization
 * \brief Helper class constructing the dual grid finite volume geometries
 *        for the pq3 cvfe discretization method
 */
#ifndef DUMUX_DISCRETIZATION_PQ3_GEOMETRY_HELPER_HH
#define DUMUX_DISCRETIZATION_PQ3_GEOMETRY_HELPER_HH

#include <array>

#include <dune/common/exceptions.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/common/reservedvector.hh>

#include <dumux/common/math.hh>
#include <dumux/geometry/volume.hh>
#include <dumux/discretization/box/boxgeometryhelper.hh>
// Reuse PQ2 corner storage traits (same structure, different order)
#include <dumux/discretization/pq2/geometryhelper.hh>

namespace Dumux {

/*!
 * \ingroup PQ3Discretization
 * \brief A class to create sub control volume and sub control volume face geometries per element
 *        for the order-3 hybrid CVFE scheme.
 *
 * Control volumes are defined only for vertex DOFs (same box dual mesh as PQ2).
 * Edge, face, and element interior DOFs are non-CV DOFs.
 * The key differences from PQ2:
 *  - dofIndex adds localKey.index() to handle multiple DOFs per edge/face/element entity
 *  - localDofPosition returns correct equidistant Lagrange node positions for order 3
 */
template <class GridView, class ScvType, class ScvfType>
class HybridPQ3GeometryHelper
{
    using Scalar = typename GridView::ctype;
    using GlobalPosition = typename Dune::FieldVector<Scalar, GridView::dimensionworld>;
    using ScvCornerStorage = typename ScvType::Traits::CornerStorage;
    using ScvfCornerStorage = typename ScvfType::Traits::CornerStorage;
    using LocalIndexType = typename ScvType::Traits::LocalIndexType;

    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    static constexpr auto dim = GridView::dimension;
    static constexpr auto dimWorld = GridView::dimensionworld;

    using BoxHelper = Dumux::BoxGeometryHelper<GridView, dim, ScvType, ScvfType>;
public:

    HybridPQ3GeometryHelper(const typename Element::Geometry& geometry)
    : geo_(geometry)
    , boxHelper_(geometry)
    {}

    //! Create a vector with the scv corners
    ScvCornerStorage getScvCorners(unsigned int localScvIdx) const
    {
        return getScvCorners(geo_.type(), [&](const auto& local){ return geo_.global(local); }, localScvIdx);
    }

    //! Create a vector with the scv corners
    template<class Transformation>
    static ScvCornerStorage getScvCorners(Dune::GeometryType type, Transformation&& trans, unsigned int localScvIdx)
    {
        const auto& ref = Dune::referenceElement<Scalar, dim>(type);
        const auto numBoxScv = ref.size(dim);
        if (localScvIdx < numBoxScv)
            return BoxHelper::getScvCorners(type, trans, localScvIdx);

        DUNE_THROW(Dune::NotImplemented, "PQ3 scv corners call for hybrid dofs");
    }

    Dune::GeometryType getScvGeometryType(unsigned int localScvIdx) const
    {
        const auto numBoxScv = boxHelper_.numScv();
        if (localScvIdx < numBoxScv)
            return Dune::GeometryTypes::cube(dim);

        DUNE_THROW(Dune::NotImplemented, "PQ3 scv geometry call for hybrid dofs");
    }

    //! Create a vector with the corners of sub control volume faces
    ScvfCornerStorage getScvfCorners(unsigned int localScvfIdx) const
    {
        return getScvfCorners(geo_.type(), [&](const auto& local){ return geo_.global(local); }, localScvfIdx);
    }

    //! Create a vector with the corners of sub control volume faces
    template<class Transformation>
    static ScvfCornerStorage getScvfCorners(Dune::GeometryType type, Transformation&& trans, unsigned int localScvfIdx)
    {
        const auto& ref = Dune::referenceElement<Scalar, dim>(type);
        const auto numBoxScvf = ref.size(dim-1);
        if (localScvfIdx < numBoxScvf)
            return BoxHelper::getScvfCorners(type, trans, localScvfIdx);

        DUNE_THROW(Dune::NotImplemented, "PQ3 scvf corners call for hybrid dofs");
    }

    Dune::GeometryType getInteriorScvfGeometryType(unsigned int localScvfIdx) const
    {
        const auto numBoxScvf = boxHelper_.numInteriorScvf();
        if (localScvfIdx < numBoxScvf)
            return Dune::GeometryTypes::cube(dim-1);

        DUNE_THROW(Dune::NotImplemented, "PQ3 interior scvf geometry type call for hybrid dofs");
    }

    //! Create the sub control volume face geometries on the boundary
    ScvfCornerStorage getBoundaryScvfCorners(unsigned int localFacetIndex, unsigned int indexInFacet) const
    {
        return boxHelper_.getBoundaryScvfCorners(localFacetIndex, indexInFacet);
    }

    Dune::GeometryType getBoundaryScvfGeometryType(unsigned int localScvfIdx) const
    {
        return Dune::GeometryTypes::cube(dim-1);
    }

    template<int d = dimWorld, std::enable_if_t<(d==3), int> = 0>
    GlobalPosition normal(const ScvfCornerStorage& p, const std::array<LocalIndexType, 2>& scvPair)
    {
        auto normal = Dumux::crossProduct(p[1]-p[0], p[2]-p[0]);
        normal /= normal.two_norm();

        GlobalPosition v = geo_.corner(scvPair[1]) - geo_.corner(scvPair[0]);

        const auto s = v*normal;
        if (std::signbit(s))
            normal *= -1;

        return normal;
    }

    template<int d = dimWorld, std::enable_if_t<(d==2), int> = 0>
    GlobalPosition normal(const ScvfCornerStorage& p, const std::array<LocalIndexType, 2>& scvPair)
    {
        const auto t = p[1] - p[0];
        GlobalPosition normal({-t[1], t[0]});
        normal /= normal.two_norm();

        GlobalPosition v = geo_.corner(scvPair[1]) - geo_.corner(scvPair[0]);

        const auto s = v*normal;
        if (std::signbit(s))
            normal *= -1;

        return normal;
    }

    //! the wrapped element geometry
    const typename Element::Geometry& elementGeometry() const
    { return geo_; }

    //! number of interior sub control volume faces
    static auto numInteriorScvf(Dune::GeometryType type)
    {
        return BoxHelper::numInteriorScvf(type);
    }

    //! number of boundary sub control volume faces for face localFacetIndex
    static auto numBoundaryScvf(Dune::GeometryType type, unsigned int localFacetIndex)
    {
        return Dune::referenceElement<Scalar, dim>(type).size(localFacetIndex, 1, dim);
    }

    //! number of sub control volumes (one per vertex)
    std::size_t numScv() const
    {
        return boxHelper_.numScv();
    }

    //! get scv volume
    Scalar scvVolume(unsigned int localScvIdx, const ScvCornerStorage& p) const
    {
        const auto scvType = getScvGeometryType(localScvIdx);
        return Dumux::convexPolytopeVolume<dim>(
            scvType,
            [&](unsigned int i){ return p[i]; }
        );
    }

    //! Local dof index related to a localDof, with index localDofIdx, on an intersection with index iIdx
    template<class LocalKey>
    static auto localDofOnIntersection(Dune::GeometryType type, unsigned int iIdx, const LocalKey& localKey)
    {
        const auto& refElement = Dune::referenceElement<Scalar, dim>(type);
        const auto numEntitiesIntersection = refElement.size(iIdx, 1, localKey.codim());
        for (std::size_t idx = 0; idx < numEntitiesIntersection; idx++)
            if (localKey.subEntity() == refElement.subEntity(iIdx, 1, idx, localKey.codim()))
                return true;
        return false;
    }

    //! Global DOF index with orientation-consistent permutation for higher-order Lagrange DOFs.
    //!
    //! Two adjacent elements sharing an edge (2D/3D) or face (3D) may assign
    //! different local orderings to the shared sub-entity's interior DOFs.
    //! This function maps every (element, localKey) pair to a globally consistent
    //! DOF index by applying the same orientation logic as dune-functions LagrangeBasis.
    //!
    //! Edge DOFs (codim == dim-1): flip the intra-entity index when the edge's
    //! first reference-element vertex has a larger global index than its second.
    //! For k=3 this means: index → 1 - index.
    //!
    //! Quad face DOFs in 3D (codim == 1, quad face): apply one of 8 rotation/
    //! reflection permutations determined by the face's global vertex ordering.
    //! Triangle face DOFs (1 interior DOF for k=3): no permutation needed.
    //!
    //! The GlobalIdSet (from gridView.grid().globalIdSet()) must be passed to
    //! obtain globally unique, MPI-consistent vertex identifiers. Using rank-local
    //! DOF mapper indices here would cause different flip decisions on different ranks
    //! for shared border edges, leading to wrong parallel assembly.
    //!
    //! Adapted from dune-functions Experimental::FaceOrientations and
    //! LagrangeFaceDOFPermutation (dune/functions/functionspacebases/lagrangebasis.hh).
    //! See license remarks and copyright notice at the top of this file.
    template<class DofMapper, class LocalKey, class GlobalIdSet>
    static auto dofIndex(const DofMapper& dofMapper, const Element& element,
                         const LocalKey& localKey, const GlobalIdSet& gidSet)
    {
        const auto baseIndex = dofMapper.subIndex(element, localKey.subEntity(), localKey.codim());
        const auto& ref = Dune::referenceElement<Scalar, dim>(element.type());

        // Vertex (codim==dim) and cell interior (codim==0): no permutation
        if (localKey.codim() == 0 || localKey.codim() == dim)
            return baseIndex + localKey.index();

        // Edge DOFs (codim==dim-1 in both 2D and 3D):
        // flip intra-entity index when local edge direction disagrees with global orientation.
        // Use gidSet (globally unique across MPI ranks) — NOT dofMapper.subIndex which
        // gives rank-local indices that can differ on different processes for the same entity.
        if (localKey.codim() == dim - 1)
        {
            const int rv0 = ref.subEntity(localKey.subEntity(), dim-1, 0, dim);
            const int rv1 = ref.subEntity(localKey.subEntity(), dim-1, 1, dim);
            const bool flip = gidSet.subId(element, rv0, dim)
                            > gidSet.subId(element, rv1, dim);
            return baseIndex + (flip ? (1 - localKey.index()) : localKey.index());
        }

        // 3D face DOFs (codim==1): permute based on face orientation.
        if constexpr (dim == 3)
        {
            if (localKey.codim() == 1)
            {
                // P3 on simplex: 1 interior DOF per triangular face, no permutation needed.
                if (ref.type(localKey.subEntity(), 1).isTriangle())
                    return baseIndex + localKey.index();

                // Q3 on cube: 4 interior DOFs per quad face — apply rotation/reflection.
                if (ref.type(localKey.subEntity(), 1).isQuadrilateral())
                {
                    const auto j = quadFaceOrientIdx_(element, localKey.subEntity(), ref, gidSet);
                    return baseIndex + quadFaceDofPermutation_[j][localKey.index()];
                }
            }
        }

        return baseIndex + localKey.index();
    }

private:
    // Permutation table for Q3 quad face interior DOFs (k=3, 4 DOFs, 8 orientations).
    // quadFaceDofPermutation_[j][i] = globally-oriented index for local DOF i, orientation j.
    // Row j encodes: bits 0-1 of j = number of CCW rotations, bit 2 of j = reflection.
    // Adapted from dune-functions LagrangeFaceDOFPermutation::globallyOrientedQuadrilateralDOFTable(3).
    // See license remarks and copyright notice at the top of this file.
    static constexpr std::array<std::array<unsigned char, 4>, 8> quadFaceDofPermutation_ = {{
        {0, 1, 2, 3},  // j=0: identity
        {1, 3, 0, 2},  // j=1: 1 CCW rotation
        {3, 2, 1, 0},  // j=2: 2 CCW rotations (180°)
        {2, 0, 3, 1},  // j=3: 3 CCW rotations
        {0, 2, 1, 3},  // j=4: reflection
        {2, 3, 0, 1},  // j=5: 1 CCW rotation + reflection
        {3, 1, 2, 0},  // j=6: 2 CCW rotations + reflection
        {1, 0, 3, 2},  // j=7: 3 CCW rotations + reflection
    }};

    // Compute orientation index j in 0..7 for a quad face (codim=1 in 3D).
    // j encodes: bits 0-1 = CCW rotations needed to align min-vertex to position 0,
    //            bit 2 = additional reflection.
    // Uses gidSet for globally consistent vertex ordering (not rank-local mapper indices).
    // Adapted from dune-functions FaceOrientations::computeQuadrilateralOrientation.
    // See license remarks and copyright notice at the top of this file.
    template<class GlobalIdSet>
    static unsigned int quadFaceOrientIdx_(const Element& element,
                                           unsigned int faceSubEntity, const auto& ref,
                                           const GlobalIdSet& gidSet)
    {
        // Global vertex ids of the face's 4 corners in reference-element order.
        std::array<std::size_t, 4> vg;
        for (int k = 0; k < 4; ++k)
            vg[k] = gidSet.subId(element, ref.subEntity(faceSubEntity, 1, k, dim), dim);

        // Edge flip bits for the 4 edges of this face.
        // ref.subEntity(faceSubEntity, 1, e, dim-1) gives the element-level edge index.
        const auto edgeFlip = [&](int e) -> bool {
            const int edgeIdx = ref.subEntity(faceSubEntity, 1, e, dim-1);
            return gidSet.subId(element, ref.subEntity(edgeIdx, dim-1, 0, dim), dim)
                 > gidSet.subId(element, ref.subEntity(edgeIdx, dim-1, 1, dim), dim);
        };
        const std::size_t eo = (std::size_t(edgeFlip(0)) << 0) | (std::size_t(edgeFlip(1)) << 1)
                             | (std::size_t(edgeFlip(2)) << 2) | (std::size_t(edgeFlip(3)) << 3);

        // Determine the face vertex with the smallest global index (i_min).
        // Lookup table maps 4-bit edge orientation bitfield to i_min (2 bits each, 16 entries).
        // Special cases eo=5 and eo=10 need an extra vertex comparison.
        constexpr uint32_t eoToImin = 0b11'11'01'01'11'00'00'00'10'00'00'01'10'00'10'00;
        std::size_t i_min;
        if (eo == 5)       i_min = (vg[1] < vg[2]) ? 1 : 2;
        else if (eo == 10) i_min = (vg[0] < vg[3]) ? 0 : 3;
        else               i_min = (eoToImin >> (2 * eo)) & 3;

        // Encode orientation index j = (rotations in bits 0-1) | (reflection in bit 2).
        if (i_min == 0) return 0u | (unsigned((vg[2] < vg[1])) << 2);
        if (i_min == 1) return 3u | (unsigned((vg[0] < vg[3])) << 2);
        if (i_min == 2) return 1u | (unsigned((vg[3] < vg[0])) << 2);
      /*i_min == 3*/   return 2u | (unsigned((vg[1] < vg[2])) << 2);
    }

public:

    template<class Geometry, class LocalKey>
    static GlobalPosition dofPosition(const Geometry& geo, const LocalKey& localKey)
    {
        if (localKey.codim() == dim)
            return geo.corner(localKey.subEntity());
        else
            return geo.global(localDofPosition(geo.type(), localKey));
    }

    template<class LocalKey>
    GlobalPosition dofPosition(const LocalKey& localKey) const
    {
        return dofPosition(geo_, localKey);
    }

    //! Local coordinate of a DOF for order-3 Lagrange basis
    template<class LocalKey>
    static typename Element::Geometry::LocalCoordinate localDofPosition(Dune::GeometryType type, const LocalKey& localKey)
    {
        using LocalCoord = typename Element::Geometry::LocalCoordinate;
        const auto ref = Dune::referenceElement<Scalar, dim>(type);

        // Vertices: exact reference element corner positions
        if (localKey.codim() == dim)
            return ref.position(localKey.subEntity(), dim);

        // Edges (codim == dim-1): 2 equidistant interior nodes at (j+1)/3 from edge start
        if (localKey.codim() == dim - 1)
        {
            const int v0Idx = ref.subEntity(localKey.subEntity(), dim-1, 0, dim);
            const int v1Idx = ref.subEntity(localKey.subEntity(), dim-1, 1, dim);
            auto pos = LocalCoord(ref.position(v0Idx, dim));
            const auto dir = LocalCoord(ref.position(v1Idx, dim)) - pos;
            pos.axpy(Scalar(localKey.index() + 1) / Scalar(3), dir);
            return pos;
        }

        // Face interior DOFs (3D only, codim=1)
        if constexpr (dim == 3)
        {
            if (localKey.codim() == 1)
            {
                const auto faceGeom = ref.template geometry<1>(localKey.subEntity());
                if (type.isSimplex())
                {
                    // P3 tet: triangle faces have 1 interior node at the face centroid
                    // Face centroid in face-local coords for a reference triangle: (1/3, 1/3)
                    const Dune::FieldVector<Scalar, 2> c{Scalar(1)/Scalar(3), Scalar(1)/Scalar(3)};
                    return faceGeom.global(c);
                }
                else
                {
                    // Q3 hex: quad faces have 4 interior nodes at tensor product positions
                    const auto ix = localKey.index() % 2;
                    const auto iy = localKey.index() / 2;
                    const Dune::FieldVector<Scalar, 2> c{Scalar(ix+1)/Scalar(3), Scalar(iy+1)/Scalar(3)};
                    return faceGeom.global(c);
                }
            }
        }

        // Element interior DOFs (codim == 0)
        if (localKey.codim() == 0)
        {
            if (type.isSimplex())
            {
                // P3: 1 interior node for triangle (at centroid), 0 for tet
                return ref.position(0, 0);
            }
            else
            {
                // Q3 cube: tensor product of interior nodes at (i+1)/3 per direction
                if constexpr (dim == 1)
                {
                    return LocalCoord{Scalar(localKey.index() + 1) / Scalar(3)};
                }
                else if constexpr (dim == 2)
                {
                    const auto ix = localKey.index() % 2;
                    const auto iy = localKey.index() / 2;
                    return LocalCoord{Scalar(ix+1)/Scalar(3), Scalar(iy+1)/Scalar(3)};
                }
                else if constexpr (dim == 3)
                {
                    const auto ix = localKey.index() % 2;
                    const auto iy = (localKey.index() / 2) % 2;
                    const auto iz = localKey.index() / 4;
                    return LocalCoord{Scalar(ix+1)/Scalar(3), Scalar(iy+1)/Scalar(3), Scalar(iz+1)/Scalar(3)};
                }
            }
        }

        DUNE_THROW(Dune::NotImplemented, "PQ3 localDofPosition for codim=" << localKey.codim());
    }

    std::array<LocalIndexType, 2> getScvPairForScvf(unsigned int localScvfIndex) const
    {
        const auto numBoxFaces = boxHelper_.numInteriorScvf();
        if (localScvfIndex < numBoxFaces)
        {
            return {
                static_cast<LocalIndexType>(referenceElement(geo_).subEntity(localScvfIndex, dim-1, 0, dim)),
                static_cast<LocalIndexType>(referenceElement(geo_).subEntity(localScvfIndex, dim-1, 1, dim))
            };
        }

        DUNE_THROW(Dune::NotImplemented, "PQ3 scv pair call for hybrid dofs");
    }

    std::array<LocalIndexType, 2> getScvPairForBoundaryScvf(unsigned int localFacetIndex, unsigned int localIsScvfIndex) const
    {
        const auto numBoxScvf = referenceElement(geo_).size(localFacetIndex, 1, dim);
        if (localIsScvfIndex < numBoxScvf)
        {
            const LocalIndexType insideScvIdx
                = static_cast<LocalIndexType>(referenceElement(geo_).subEntity(localFacetIndex, 1, localIsScvfIndex, dim));
            return { insideScvIdx, insideScvIdx };
        }

        DUNE_THROW(Dune::NotImplemented, "PQ3 scv boundary pair call for hybrid dofs");
    }

    bool isOverlappingScvf(unsigned int localScvfIndex) const
    { return false; }

    bool isOverlappingBoundaryScvf(unsigned int localFacetIndex) const
    { return false; }

    bool isOverlappingScv(unsigned int localScvIndex) const
    { return false; }

private:
    const typename Element::Geometry& geo_;
    BoxHelper boxHelper_;
};

} // end namespace Dumux

#endif
