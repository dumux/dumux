// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEDiscretization
 * \brief Class representing dofs on elements for control-volume finite element schemes
 */
#ifndef DUMUX_CVFE_LOCAL_DOF_HH
#define DUMUX_CVFE_LOCAL_DOF_HH

#include <ranges>

#include <dumux/common/indextraits.hh>

namespace Dumux::CVFE {

/*!
 * \ingroup CVFEDiscretization
 * \brief A local degree of freedom from an element perspective
 */
template<class LocalIndex, class GridIndex>
class LocalDof
{
public:
    using LocalIndexType = LocalIndex;
    using GridIndexType = GridIndex;

    LocalDof(LocalIndex index, GridIndex dofIndex, GridIndex eIdx)
    : index_(index), dofIndex_(dofIndex), eIdx_(eIdx) {}
    LocalIndex index() const { return index_; }
    GridIndex dofIndex() const { return dofIndex_; }
    GridIndex elementIndex() const { return eIdx_; }
private:
    LocalIndex index_;
    GridIndex dofIndex_;
    GridIndex eIdx_;
};


} // end namespace Dumux::CVFE

namespace Dumux {

//! range over local dofs
template<class FVElementGeometry>
inline auto localDofs(const FVElementGeometry& fvGeometry)
{
    using LocalIndexType = typename IndexTraits<typename FVElementGeometry::GridGeometry::GridView>::LocalIndex;
    using GridIndexType = typename IndexTraits<typename FVElementGeometry::GridGeometry::GridView>::GridIndex;

    return std::views::iota(0u, fvGeometry.numScv())
        | std::views::transform([&](const auto i) { return CVFE::LocalDof
            {
                static_cast<LocalIndexType>(i),
                static_cast<GridIndexType>(fvGeometry.scv(i).dofIndex()),
                static_cast<GridIndexType>(fvGeometry.scv(i).elementIndex())
            }; }
        );
}

//! range over control-volume local dofs
template<class FVElementGeometry>
inline auto cvLocalDofs(const FVElementGeometry& fvGeometry)
{
    // As default all dofs are cv dofs
    return localDofs(fvGeometry);
}

//! range over sub control volumes related to a local dof.
//! this is the default range where a one-to-one mapping between scvs and localDofs is assumed.
//! If multiple scvs are related to a localDof, this range needs to be overwritten
//! within the fvElementGeometry class
template<class FVElementGeometry, class LocalDof>
inline auto
scvs(const FVElementGeometry& fvGeometry, const LocalDof& localDof)
{
    assert(fvGeometry.numScv() > localDof.index());
    return std::views::single(1) | std::views::transform([&](const auto i) { return fvGeometry.scv(localDof.index()); });
}

} // end namespace Dumux
#endif
