// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ1NonconformingDiscretization
 * \brief Helper class providing degree of freedom information for the PQ1 nonconforming scheme.
 */
#ifndef DUMUX_DISCRETIZATION_PQ1_NONCONFORMING_DOF_HELPER_HH
#define DUMUX_DISCRETIZATION_PQ1_NONCONFORMING_DOF_HELPER_HH

#include <ranges>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/fem/fedofhelper.hh>
#include <dumux/discretization/cvfe/localdof.hh>

namespace Dumux {

/*!
 * \ingroup PQ1NonconformingDiscretization
 * \brief Helper class providing degree of freedom information for the PQ1 nonconforming scheme.
 *        It uses face-centered dofs: exactly one dof per face/intersection.
 */
template<class GridView>
class PQ1NonconformingDofHelper : public Dumux::Experimental::FEDofHelper<GridView>
{
    using LocalIndexType = typename IndexTraits<GridView>::SmallLocalIndex;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
public:
    /*!
     * \brief Iterator range over all local dofs on a given boundary face.
     *       For the PQ1 nonconforming scheme, each boundary face has exactly one dof.
     */
    template<class ElemDisc, class BoundaryFace>
    static auto localDofsOnBoundaryFace(const ElemDisc& elemDisc, const BoundaryFace& boundaryFace)
    {
        const auto localDofIdx = boundaryFace.intersectionIndex();
        return std::views::single(CVFE::LocalDof(
            static_cast<LocalIndexType>(localDofIdx),
            static_cast<GridIndexType>(elemDisc.scv(localDofIdx).dofIndex()),
            static_cast<GridIndexType>(elemDisc.elementIndex())
        ));
    }
};

} // namespace Dumux

#endif // DUMUX_DISCRETIZATION_PQ1_NONCONFORMING_DOF_HELPER_HH
