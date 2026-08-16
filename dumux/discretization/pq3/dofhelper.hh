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
 * \brief Lightweight DOF helper for order-3 Lagrange elements.
 *
 * Provides `dofIndex` and `dofPosition` as pure static methods, including
 * the orientation-consistent permutation logic needed for shared edge and
 * face DOFs (adapted from HybridPQ3GeometryHelper / dune-functions).
 * The uniform 4-argument `dofIndex` signature lets PQ2 and PQ3 helpers
 * be used interchangeably by callers.
 */
#ifndef DUMUX_DISCRETIZATION_PQ3_DOF_HELPER_HH
#define DUMUX_DISCRETIZATION_PQ3_DOF_HELPER_HH

#include <array>
#include <cstdint>

#include <dune/common/fvector.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/discretization/fem/fedofhelper.hh>

namespace Dumux {

/*!
 * \ingroup PQ3Discretization
 * \brief DOF index and position helper for order-3 Lagrange discretizations.
 *
 * Edge and face DOFs require orientation-consistent permutation so that
 * adjacent elements assign the same global index to shared sub-entity DOFs.
 * The algorithm is adapted from dune-functions' LagrangeBasis orientation
 * logic (see license remarks at the top of this file).
 *
 * \tparam GridView  The Dune grid view type.
 */
template<class GridView>
class PQ3LagrangeDofHelper : public Dumux::Experimental::FEDofHelper<GridView>
{
private:
    using ParentType = Dumux::Experimental::FEDofHelper<GridView>;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using Scalar = typename GridView::ctype;
    static constexpr int dim = GridView::dimension;
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;

public:
    /*!
     * \brief Orientation-consistent global DOF index.
     *
     * For edge DOFs the within-entity index is flipped when the edge's global
     * vertex ordering disagrees with the reference-element ordering.  For 3D
     * quad-face DOFs one of 8 rotation/reflection permutations is applied.
     * A GlobalIdSet must be passed so that the orientation decision is
     * consistent across MPI ranks.
     */
    template<class DofMapper, class Element, class LocalKey, class IdSet>
    static std::size_t dofIndex(const DofMapper& m, const Element& e,
                                const LocalKey& lk, const IdSet& idSet)
    {
        const auto base = m.subIndex(e, lk.subEntity(), lk.codim());
        const auto& ref = Dune::referenceElement<Scalar, dim>(e.type());

        if (lk.codim() == 0 || lk.codim() == dim)
            return base + lk.index();

        if (lk.codim() == dim - 1)
        {
            const int rv0 = ref.subEntity(lk.subEntity(), dim-1, 0, dim);
            const int rv1 = ref.subEntity(lk.subEntity(), dim-1, 1, dim);
            const bool flip = idSet.subId(e, rv0, dim) > idSet.subId(e, rv1, dim);
            return base + (flip ? 1 - (int)lk.index() : (int)lk.index());
        }

        if constexpr (dim == 3)
        {
            if (lk.codim() == 1)
            {
                if (ref.type(lk.subEntity(), 1).isTriangle())
                    return base + lk.index();
                if (ref.type(lk.subEntity(), 1).isQuadrilateral())
                {
                    const auto j = quadFaceOrientIdx_(e, lk.subEntity(), ref, idSet);
                    return base + quadFacePerm_[j][lk.index()];
                }
            }
        }

        return base + lk.index();
    }

    /*!
     * \brief Iterator range over all local dofs on an element.
     * \param elemDisc the element discretization (must be bound)
     */
    template<class ElemDisc>
    static auto localDofs(const ElemDisc& elemDisc)
    {
        const auto& gridDisc = elemDisc.gridDiscretization();

        return Dune::transformedRangeView(
            Dune::range(elemDisc.numLocalDofs()),
            [&](const auto i) {
                return CVFE::LocalDof{
                    static_cast<LocalIndexType>(i),
                    static_cast<GridIndexType>(dofIndex(
                        gridDisc.dofMapper(),
                        elemDisc.element(),
                        elemDisc.feLocalCoefficients().localKey(i),
                        gridDisc.gridView().grid().globalIdSet())),
                    static_cast<GridIndexType>(elemDisc.elementIndex())
                };
            }
        );
    }

    /*!
     * \brief Iterator range over all local dofs on a given boundary face.
     *         Uses a filter over all local dofs via localDofOnIntersection.
     * \param elemDisc the element discretization (must be bound)
     * \param boundaryFace the boundary face
     */
    template<class ElemDisc, class BoundaryFace>
    static auto localDofsOnBoundaryFace(const ElemDisc& elemDisc, const BoundaryFace& boundaryFace)
    {
        const auto& gridDisc = elemDisc.gridDiscretization();

        return std::views::iota(std::size_t(0), elemDisc.numLocalDofs())
            | std::views::filter([&](std::size_t i) {
                return ParentType::localDofOnIntersection(
                    elemDisc.element().type(),
                    boundaryFace.intersectionIndex(),
                    elemDisc.feLocalCoefficients().localKey(i));
            })
            | std::views::transform([&](std::size_t i) {
                return CVFE::LocalDof(
                    static_cast<LocalIndexType>(i),
                    static_cast<GridIndexType>(dofIndex(
                        gridDisc.dofMapper(),
                        elemDisc.element(),
                        elemDisc.feLocalCoefficients().localKey(i),
                        gridDisc.gridView().grid().globalIdSet())),
                    static_cast<GridIndexType>(elemDisc.elementIndex())
                );
            });
    }

    //! Physical position of a DOF in global coordinates.
    template<class Geometry, class LocalKey>
    static GlobalPosition dofPosition(const Geometry& geo, const LocalKey& lk)
    { return geo.global(localDofPos_(geo.type(), lk)); }

    //! Reference-element position of a DOF for order-3 Lagrange basis.
    template<class LocalKey>
    static typename GridView::template Codim<0>::Entity::Geometry::LocalCoordinate
    localDofPosition(Dune::GeometryType gt, const LocalKey& lk)
    { return localDofPos_(gt, lk); }

private:
    // Permutation table for Q3 quad-face interior DOFs (4 DOFs, 8 orientations).
    // Adapted from dune-functions LagrangeFaceDOFPermutation::globallyOrientedQuadrilateralDOFTable(3).
    static constexpr std::array<std::array<unsigned char, 4>, 8> quadFacePerm_ = {{
        {0,1,2,3}, {1,3,0,2}, {3,2,1,0}, {2,0,3,1},
        {0,2,1,3}, {2,3,0,1}, {3,1,2,0}, {1,0,3,2},
    }};

    // Compute orientation index j in [0,8) for a quad face (3D, codim=1).
    // Adapted from dune-functions FaceOrientations::computeQuadrilateralOrientation.
    template<class Element, class IdSet>
    static unsigned int quadFaceOrientIdx_(const Element& e, unsigned int face,
                                           const auto& ref, const IdSet& idSet)
    {
        std::array<typename IdSet::IdType, 4> vg;
        for (int i = 0; i < 4; ++i)
            vg[i] = idSet.subId(e, ref.subEntity(face, 1, i, dim), dim);

        const auto flip = [&](int ed) {
            const int ei = ref.subEntity(face, 1, ed, dim-1);
            return idSet.subId(e, ref.subEntity(ei, dim-1, 0, dim), dim)
                 > idSet.subId(e, ref.subEntity(ei, dim-1, 1, dim), dim);
        };
        const std::size_t eo = (std::size_t(flip(0))<<0) | (std::size_t(flip(1))<<1)
                             | (std::size_t(flip(2))<<2) | (std::size_t(flip(3))<<3);

        constexpr uint32_t eoToImin = 0b11'11'01'01'11'00'00'00'10'00'00'01'10'00'10'00;
        std::size_t i_min;
        if (eo == 5)       i_min = (vg[1] < vg[2]) ? 1 : 2;
        else if (eo == 10) i_min = (vg[0] < vg[3]) ? 0 : 3;
        else               i_min = (eoToImin >> (2*eo)) & 3;

        if (i_min == 0) return 0u | (unsigned(vg[2] < vg[1]) << 2);
        if (i_min == 1) return 3u | (unsigned(vg[0] < vg[3]) << 2);
        if (i_min == 2) return 1u | (unsigned(vg[3] < vg[0]) << 2);
        /*i_min==3*/    return 2u | (unsigned(vg[1] < vg[2]) << 2);
    }

    template<class LocalKey>
    static typename GridView::template Codim<0>::Entity::Geometry::LocalCoordinate
    localDofPos_(Dune::GeometryType gt, const LocalKey& lk)
    {
        using LocalCoord = typename GridView::template Codim<0>::Entity::Geometry::LocalCoordinate;
        const auto ref = Dune::referenceElement<Scalar, dim>(gt);

        if (lk.codim() == dim)
            return ref.position(lk.subEntity(), dim);

        if (lk.codim() == dim - 1)
        {
            const int v0 = ref.subEntity(lk.subEntity(), dim-1, 0, dim);
            const int v1 = ref.subEntity(lk.subEntity(), dim-1, 1, dim);
            auto pos = LocalCoord(ref.position(v0, dim));
            const auto dir = LocalCoord(ref.position(v1, dim)) - pos;
            pos.axpy(Scalar(lk.index() + 1) / Scalar(3), dir);
            return pos;
        }

        if constexpr (dim == 3)
        {
            if (lk.codim() == 1)
            {
                const auto faceGeom = ref.template geometry<1>(lk.subEntity());
                if (gt.isSimplex())
                {
                    const Dune::FieldVector<Scalar, 2> c{Scalar(1)/3, Scalar(1)/3};
                    return faceGeom.global(c);
                }
                else
                {
                    const auto ix = lk.index() % 2, iy = lk.index() / 2;
                    const Dune::FieldVector<Scalar, 2> c{Scalar(ix+1)/3, Scalar(iy+1)/3};
                    return faceGeom.global(c);
                }
            }
        }

        if (lk.codim() == 0)
        {
            if (gt.isSimplex())
                return ref.position(0, 0);

            if constexpr (dim == 2)
            {
                const auto ix = lk.index() % 2, iy = lk.index() / 2;
                return LocalCoord{Scalar(ix+1)/3, Scalar(iy+1)/3};
            }
            else if constexpr (dim == 3)
            {
                const auto ix = lk.index()%2, iy = (lk.index()/2)%2, iz = lk.index()/4;
                return LocalCoord{Scalar(ix+1)/3, Scalar(iy+1)/3, Scalar(iz+1)/3};
            }
        }

        DUNE_THROW(Dune::NotImplemented, "PQ3 local DOF position for codim=" << lk.codim());
    }
};

} // namespace Dumux

#endif // DUMUX_DISCRETIZATION_PQ3_DOF_HELPER_HH
