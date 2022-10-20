// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup Discretization
 * \copydoc Dumux::FVGridGeometryStorage
 */
#ifndef DUMUX_DISCRETIZATION_FV_GRID_GEOMETRY_STORAGE_HH
#define DUMUX_DISCRETIZATION_FV_GRID_GEOMETRY_STORAGE_HH

#include <utility>

#include <dumux/experimental/new_assembly/dumux/common/size.hh>
#include <dumux/experimental/new_assembly/dumux/common/storage.hh>

namespace Dumux {

struct DefaultFVGridGeometryStorageSizes
{
    static constexpr auto maxNumScv = Size::dynamic;
    static constexpr auto maxNumScvf = Size::dynamic;
    static constexpr auto maxNumFaces = Size::dynamic;
};

/*!
 * \ingroup Discretization
 * \brief Class to store finite-volume grid geometries.
 *        Allows insertion of geometries while retrieving their local storage
 *        index. Access to the geometries afterwards occurs via their index.
 */
template<typename F,
         typename Scv,
         typename Scvf,
         typename Sizes = DefaultFVGridGeometryStorageSizes>
class FVGridGeometryStorage
{
public:
    using Face = F;
    using SubControlVolume = Scv;
    using SubControlVolumeFace = Scvf;

    struct Size
    {
        std::size_t numFaces;
        std::size_t numScvfs;
        std::size_t numScvs;
    };

    void reserve(const Size& size)
    {
        faces_.reserve(size.numFaces);
        scvfs_.reserve(size.numScvfs);
        scvs_.reserve(size.numScvs);
    }

    void clear()
    {
        faces_.clear();
        scvfs_.clear();
        scvs_.clear();
    }

    void shrinkToFit()
    {
        faces_.shrink_to_fit();
        scvfs_.shrink_to_fit();
        scvs_.shrink_to_fit();
    }

    std::size_t numFaces() const { return faces_.size(); }
    std::size_t numScvfs() const { return scvfs_.size(); }
    std::size_t numScvs() const { return scvs_.size(); }

    std::size_t pushFace(Face&& f) { return add_(faces_, std::move(f)); }
    std::size_t pushScvf(Scvf&& f) { return add_(scvfs_, std::move(f)); }
    std::size_t pushScv(Scv&& s) { return add_(scvs_, std::move(s)); }

    const Face& face(const std::size_t i) const { return faces_[i]; }
    const Scvf& scvf(const std::size_t i) const { return scvfs_[i]; }
    const Scv& scv(const std::size_t i) const { return scvs_[i]; }

    Face& face(const std::size_t i) { return faces_[i]; }
    Scvf& scvf(const std::size_t i) { return scvfs_[i]; }
    Scv& scv(const std::size_t i) { return scvs_[i]; }

private:
    template<typename Container, typename Entity>
    std::size_t add_(Container& c, Entity&& e)
    {
        const std::size_t nextIdx = c.size();
        c.push_back(std::move(e));
        return nextIdx;
    }

    DefaultStorage<Face, Sizes::maxNumFaces> faces_;
    DefaultStorage<Scvf, Sizes::maxNumScvf> scvfs_;
    DefaultStorage<Scv, Sizes::maxNumScv> scvs_;
};

} // end namespace Dumux

#endif
