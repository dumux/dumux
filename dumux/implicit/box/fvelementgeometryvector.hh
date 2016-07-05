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
 * \brief Base class for the finite volume geometry vector for box models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 */
#ifndef DUMUX_IMPLICIT_BOX_FV_GEOMETRY_VECTOR_HH
#define DUMUX_IMPLICIT_BOX_FV_GEOMETRY_VECTOR_HH

#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/implicit/box/properties.hh>
#include <dumux/implicit/subcontrolvolume.hh>
#include <dumux/implicit/subcontrolvolumeface.hh>
#include <dumux/implicit/fvelementgeometry.hh>
#include <dumux/common/elementmap.hh>
#include <dumux/common/math.hh>

namespace Dumux
{

//! A class to create sub control volume and sub control volume face geometries per element
template <class GridView>
class BoxGeometryHelper
{
private:
    using Scalar = typename GridView::ctype;
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using ScvGeometry = Dune::MultiLinearGeometry<Scalar, dim, dimWorld>;
    using ScvfGeometry = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld>;

    using ScvGeometryVector = std::vector<std::shared_ptr<ScvGeometry>>;
    using ScvfGeometryVector = std::vector<std::shared_ptr<ScvfGeometry>>;
    using GeometryPair = std::pair<ScvGeometryVector, ScvfGeometryVector>;

    using GlobalPosition = typename ScvGeometry::GlobalCoordinate;
    using CornerList = std::vector<GlobalPosition>;
    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    using ReferenceElements = typename Dune::ReferenceElements<Scalar, dim>;
    using FaceReferenceElements = typename Dune::ReferenceElements<Scalar, dim-1>;

public:
    //! get sub control volume geometries from element of dimension 1
    template <int d = dim>
    static typename std::enable_if<d == 1, GeometryPair>::type
    getScvAndScvfGeometries(const typename Element::Geometry& geometry)
    {
        // the sub control volumes
        ScvGeometryVector scvGeometries(2);
        scvGeometries[0] = std::make_shared<ScvGeometry>(Dune::GeometryType(1), std::vector<GlobalPosition>({geometry.corner(0), geometry.center()}));
        scvGeometries[1] = std::make_shared<ScvGeometry>(Dune::GeometryType(1), std::vector<GlobalPosition>({geometry.center(), geometry.corner(1)}));

        // the sub control volume faces
        ScvfGeometryVector scvfGeometries(1);
        scvfGeometries[0] = std::make_shared<ScvfGeometry>(Dune::GeometryType(0), std::vector<GlobalPosition>({geometry.center()}));

        return std::make_pair(scvGeometries, scvfGeometries);
    }

    //! get sub control volume geometries from element of dimension 2
    template <int d = dim>
    static typename std::enable_if<d == 2, GeometryPair>::type
    getScvAndScvfGeometries(const typename Element::Geometry& geometry)
    {
        // sub control volume geometries in 2D are always quadrilaterals
        Dune::GeometryType scvGeometryType; scvGeometryType.makeQuadrilateral();
        // sub control volume face geometries in 2D are always lines
        Dune::GeometryType scvfGeometryType; scvfGeometryType.makeLine();

        // extract the corners of the sub control volumes
        const auto& referenceElement = ReferenceElements::general(geometry.type());

        // vertices
        CornerList v;
        for (int i = 0; i < geometry.corners(); ++i)
            v.emplace_back(geometry.corner(i));

        // face midpoints
        CornerList f;
        for (int i = 0; i < referenceElement.size(1); ++i)
            f.emplace_back(geometry.global(referenceElement.position(i, 1)));

        auto c = geometry.center();

        // proceed according to number of corners
        switch (geometry.corners())
        {
        case 3: // triangle
        {
            // the sub control volumes
            ScvGeometryVector scvGeometries(3);
            scvGeometries[0] = std::make_shared<ScvGeometry>(scvGeometryType, std::vector<GlobalPosition>({v[0], f[0], f[1], c}));
            scvGeometries[1] = std::make_shared<ScvGeometry>(scvGeometryType, std::vector<GlobalPosition>({f[0], v[1], c, f[2]}));
            scvGeometries[2] = std::make_shared<ScvGeometry>(scvGeometryType, std::vector<GlobalPosition>({f[1], c, v[2], f[2]}));

            // the sub control volume faces
            ScvfGeometryVector scvfGeometries(3);
            scvfGeometries[0] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({f[0], c}));
            scvfGeometries[1] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({f[1], c}));
            scvfGeometries[2] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({f[2], c}));

            return std::make_pair(scvGeometries, scvfGeometries);
        }
        case 4: // quadrilateral
        {
            // the sub control volumes
            ScvGeometryVector scvGeometries(4);
            scvGeometries[0] = std::make_shared<ScvGeometry>(scvGeometryType, std::vector<GlobalPosition>({v[0], f[2], f[0], c}));
            scvGeometries[1] = std::make_shared<ScvGeometry>(scvGeometryType, std::vector<GlobalPosition>({f[2], v[1], c, f[1]}));
            scvGeometries[2] = std::make_shared<ScvGeometry>(scvGeometryType, std::vector<GlobalPosition>({f[0], c, v[2], f[3]}));
            scvGeometries[3] = std::make_shared<ScvGeometry>(scvGeometryType, std::vector<GlobalPosition>({c, f[1], f[3], v[3]}));

             // the sub control volume faces
            ScvfGeometryVector scvfGeometries(4);
            scvfGeometries[0] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({f[0], c}));
            scvfGeometries[1] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({f[1], c}));
            scvfGeometries[2] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({f[2], c}));
            scvfGeometries[3] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({f[3], c}));

            return std::make_pair(scvGeometries, scvfGeometries);
        }
        default:
            DUNE_THROW(Dune::NotImplemented, "Box scv geometries for " << geometry.type());
        }
    }

    //! get sub control volume geometries from element of dimension 3
    template <int d = dim>
    static typename std::enable_if<d == 3, GeometryPair>::type
    getScvAndScvfGeometries(const typename Element::Geometry& geometry)
    {
        // sub control volume geometries in 3D are always hexahedrons
        Dune::GeometryType scvGeometryType; scvGeometryType.makeHexahedron();
        // sub control volume face geometries in 3D are always quadrilaterals
        Dune::GeometryType scvfGeometryType; scvfGeometryType.makeQuadrilateral();

        // extract the corners of the sub control volumes
        const auto& referenceElement = ReferenceElements::general(geometry.type());

        // vertices
        CornerList v;
        for (int i = 0; i < geometry.corners(); ++i)
            v.emplace_back(geometry.corner(i));

        // edge midpoints
        CornerList e;
        for (int i = 0; i < referenceElement.size(dim-1); ++i)
            e.emplace_back(geometry.global(referenceElement.position(i, dim-1)));

        // face midpoints
        CornerList f;
        for (int i = 0; i < referenceElement.size(1); ++i)
            f.emplace_back(geometry.global(referenceElement.position(i, 1)));

        auto c = geometry.center();

        // procees according to number of corners
        // \todo prisms (corners == 6) and pyramids (corners == 5)
        switch (geometry.corners())
        {
        case 4: // tetrahedron
        {
            // sub control volumes
            ScvGeometryVector scvGeometries(4);
            scvGeometries[0] = std::make_shared<ScvGeometry>(scvGeometryType, std::vector<GlobalPosition>({v[0], e[0], e[1], f[0], e[3], f[1], f[2], c}));
            scvGeometries[1] = std::make_shared<ScvGeometry>(scvGeometryType, std::vector<GlobalPosition>({v[1], e[2], e[0], f[0], f[3], e[4], c, f[1]}));
            scvGeometries[2] = std::make_shared<ScvGeometry>(scvGeometryType, std::vector<GlobalPosition>({v[2], e[1], e[2], f[0], e[5], f[2], f[3], c}));
            scvGeometries[3] = std::make_shared<ScvGeometry>(scvGeometryType, std::vector<GlobalPosition>({v[3], e[3], e[5], f[2], e[4], f[1], f[3], c}));

            // sub control volume faces
            ScvfGeometryVector scvfGeometries(6);
            scvfGeometries[0] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({e[0], f[0], f[1], c}));
            scvfGeometries[1] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({f[0], e[1], c, f[2]}));
            scvfGeometries[2] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({e[2], f[0], f[3], c}));
            scvfGeometries[3] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({f[2], e[3], c, f[1]}));
            scvfGeometries[4] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({f[3], c, e[4], f[1]}));
            scvfGeometries[5] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({e[5], f[2], f[3], c}));

            return std::make_pair(scvGeometries, scvfGeometries);
        }
        case 8: // hexahedron
        {
            // sub control volumes
            ScvGeometryVector scvGeometries(8);
            scvGeometries[0] = std::make_shared<ScvGeometry>(scvGeometryType, std::vector<GlobalPosition>({v[0], e[6], e[4], f[4], e[0], f[2], f[0], c}));
            scvGeometries[1] = std::make_shared<ScvGeometry>(scvGeometryType, std::vector<GlobalPosition>({e[6], v[1], f[4], e[5], f[2], e[1], c, f[1]}));
            scvGeometries[2] = std::make_shared<ScvGeometry>(scvGeometryType, std::vector<GlobalPosition>({e[4], f[4], v[2], e[7], f[0], c, e[2], f[3]}));
            scvGeometries[3] = std::make_shared<ScvGeometry>(scvGeometryType, std::vector<GlobalPosition>({f[4], e[5], e[7], v[3], c, f[1], f[3], e[3]}));
            scvGeometries[4] = std::make_shared<ScvGeometry>(scvGeometryType, std::vector<GlobalPosition>({e[0], f[2], f[0], c, v[4], e[10], e[8], f[5]}));
            scvGeometries[5] = std::make_shared<ScvGeometry>(scvGeometryType, std::vector<GlobalPosition>({f[2], e[1], c, f[1], e[10], v[5], f[5], e[9]}));
            scvGeometries[6] = std::make_shared<ScvGeometry>(scvGeometryType, std::vector<GlobalPosition>({f[0], c, e[2], f[3], e[8], f[5], v[6], e[11]}));
            scvGeometries[7] = std::make_shared<ScvGeometry>(scvGeometryType, std::vector<GlobalPosition>({c, f[1], f[3], e[3], f[5], e[9], e[11], v[7]}));

            // sub control volume faces
            ScvfGeometryVector scvfGeometries(12);
            scvfGeometries[0] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({f[0], e[0], c, f[2]}));
            scvfGeometries[1] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({f[1], c, e[1], f[2]}));
            scvfGeometries[2] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({f[3], e[2], c, f[0]}));
            scvfGeometries[3] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({e[3], f[3], f[1], c}));
            scvfGeometries[4] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({f[4], e[4], c, f[0]}));
            scvfGeometries[5] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({e[5], f[4], f[1], c}));
            scvfGeometries[6] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({e[6], f[4], f[2], c}));
            scvfGeometries[7] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({f[4], e[7], c, f[3]}));
            scvfGeometries[8] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({c, f[0], f[5], e[8]}));
            scvfGeometries[9] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({e[9], f[1], f[5], c}));
            scvfGeometries[10] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({e[10], f[2], f[5], c}));
            scvfGeometries[11] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({e[11], f[5], f[3], c}));

            return std::make_pair(scvGeometries, scvfGeometries);
        }
        default:
            DUNE_THROW(Dune::NotImplemented, "Box scv geometries for " << geometry.type());
        }
    }

    //! get sub control volume geometries from element of dimension 1
    template <int d = dim>
    static typename std::enable_if<d == 1, ScvfGeometryVector>::type
    getBoundaryScvfGeometries(const typename Intersection::Geometry& geometry)
    {
        return {std::make_shared<ScvfGeometry>(ScvfGeometry(Dune::GeometryType(0), std::vector<GlobalPosition>({geometry.center()})))};
    }

    //! get sub control volume geometries from element of dimension 2
    template <int d = dim>
    static typename std::enable_if<d == 2, ScvfGeometryVector>::type
    getBoundaryScvfGeometries(const typename Intersection::Geometry& geometry)
    {
        return {std::make_shared<ScvfGeometry>(ScvfGeometry(Dune::GeometryType(1), std::vector<GlobalPosition>({geometry.corner(0), geometry.center()}))),
                std::make_shared<ScvfGeometry>(ScvfGeometry(Dune::GeometryType(1), std::vector<GlobalPosition>({geometry.center(), geometry.corner(1)})))};
    }

    //! get sub control volume geometries from element of dimension 3
    template <int d = dim>
    static typename std::enable_if<d == 3, ScvfGeometryVector>::type
    getBoundaryScvfGeometries(const typename Intersection::Geometry& geometry)
    {
        // sub control volume face geometries in 3D are always quadrilaterals
        Dune::GeometryType scvfGeometryType; scvfGeometryType.makeQuadrilateral();

        // extract the corners of the sub control volumes
        const auto& referenceElement = FaceReferenceElements::general(geometry.type());

        // vertices
        CornerList v;
        for (int i = 0; i < geometry.corners(); ++i)
            v.emplace_back(geometry.corner(i));

        // edge midpoints
        CornerList e;
        for (int i = 0; i < referenceElement.size(1); ++i)
            e.emplace_back(geometry.global(referenceElement.position(i, 1)));

        // face midpoint
        auto c = geometry.center();

        // procees according to number of corners
        switch (geometry.corners())
        {
        case 3: // triangle
        {
            ScvfGeometryVector scvfGeometries(3);
            scvfGeometries[0] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({v[0], e[0], e[1], c}));
            scvfGeometries[1] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({v[1], e[2], e[0], c}));
            scvfGeometries[2] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({v[2], e[1], e[2], c}));

            return scvfGeometries;
        }
        case 4: // quadrilateral
        {
            ScvfGeometryVector scvfGeometries(4);
            scvfGeometries[0] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({v[0], e[2], e[0], c}));
            scvfGeometries[1] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({v[1], e[1], e[2], c}));
            scvfGeometries[2] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({v[2], e[0], e[3], c}));
            scvfGeometries[3] = std::make_shared<ScvfGeometry>(scvfGeometryType, std::vector<GlobalPosition>({v[3], e[3], e[1], c}));

            return scvfGeometries;
        }
        default:
            DUNE_THROW(Dune::NotImplemented, "Box scvf boundary geometries for " << geometry.type());
        }
    }

    //! get scvf normal vector for dim == 3, dimworld == 3
    template <int d = dim, int w = dimWorld>
    static typename std::enable_if<d == 3 && w == 3, GlobalPosition>::type
    normal(const typename Element::Geometry& geometry, const ScvfGeometry& scvfGeometry)
    {
        const auto t1 = scvfGeometry.corner(1) - scvfGeometry.corner(0);
        const auto t2 = scvfGeometry.corner(2) - scvfGeometry.corner(0);
        GlobalPosition normal = Dumux::crossProduct(t1, t2);
        normal /= normal.two_norm();
        return normal;
    }

    //! get scvf normal vector for dim == 2, dimworld == 3
    template <int d = dim, int w = dimWorld>
    static typename std::enable_if<d == 2 && w == 3, GlobalPosition>::type
    normal(const typename Element::Geometry& geometry, const ScvfGeometry& scvfGeometry)
    {
        const auto v1 = geometry.corner(1) - geometry.corner(0);
        const auto v2 = geometry.corner(2) - geometry.corner(0);
        const auto v3 = Dumux::crossProduct(v1, v2);
        const auto t = scvfGeometry.corner(1) - scvfGeometry.corner(0);
        GlobalPosition normal = Dumux::crossProduct(v3, t);
        normal /= normal.two_norm();
        return normal;
    }

    //! get scvf normal vector for dim == 2, dimworld == 2
    template <int d = dim, int w = dimWorld>
    static typename std::enable_if<d == 2 && w == 2, GlobalPosition>::type
    normal(const typename Element::Geometry& geometry, const ScvfGeometry& scvfGeometry)
    {
        const auto t = scvfGeometry.corner(1) - scvfGeometry.corner(0);
        GlobalPosition normal({-t[1], t[0]});
        normal /= normal.two_norm();
        return normal;
    }

    //! get scvf normal vector for dim == 2, dimworld == 2
    template <int d = dim>
    static typename std::enable_if<d == 1, GlobalPosition>::type
    normal(const typename Element::Geometry& geometry, const ScvfGeometry& scvfGeometry)
    {
        GlobalPosition normal = geometry.corner(1) - geometry.corner(0);
        normal /= normal.two_norm();
        return normal;
    }
};

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the finite volume geometry vector for box models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 */
template<class TypeTag, bool EnableFVElementGeometryCache>
class BoxFVElementGeometryVector
{};

// specialization in case the FVElementGeometries are stored
template<class TypeTag>
class BoxFVElementGeometryVector<TypeTag, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using Element = typename GridView::template Codim<0>::Entity;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CoordScalar = typename GridView::ctype;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using FeCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    using FeLocalBasis = typename FeCache::FiniteElementType::Traits::LocalBasisType;
    using ReferenceElements = typename Dune::ReferenceElements<CoordScalar, dim>;

public:
    //! Constructor
    BoxFVElementGeometryVector(const GridView gridView)
    : gridView_(gridView), elementMap_(gridView) {}

    /* \brief Get the finite volume geometry of an element
     * \note The finite volume geometry offers iterators over the sub control volumes
     *       and the sub control volume faces of an element.
     */
    const FVElementGeometry& fvGeometry(IndexType eIdx) const
    {
        return fvGeometries_[eIdx];
    }

    //! Get a sub control volume with a global scv index
    const SubControlVolume& subControlVolume(IndexType scvIdx) const
    {
        return *scvs_[scvIdx];
    }

    //! Get a sub control volume face with a global scvf index
    const SubControlVolumeFace& subControlVolumeFace(IndexType scvfIdx) const
    {
        return *scvfs_[scvfIdx];
    }

    //! Get a local finite element basis
    const FeLocalBasis& feLocalBasis(const Dune::GeometryType& type) const
    {
        return feCache_.get(type).localBasis();
    }

    //! The total number of sub control volumes
    std::size_t numScv() const
    {
        return scvs_.size();
    }

    //! The total number of sun control volume faces
    std::size_t numScvf() const
    {
        return scvfs_.size();
    }

    // Get an element from a sub control volume contained in it
    Element element(const SubControlVolume& scv) const
    { return elementMap_.element(scv.elementIndex()); }

    // Get an element from a global element index
    Element element(IndexType eIdx) const
    { return elementMap_.element(eIdx); }

    //! update all fvElementGeometries (do this again after grid adaption)
    void update(const Problem& problem)
    {
        scvs_.clear();
        scvfs_.clear();
        fvGeometries_.clear();
        elementMap_.clear();

        // Build the SCV and SCV faces
        IndexType scvIdx = 0;
        IndexType scvfIdx = 0;
        elementMap_.resize(gridView_.size(0));

        for (const auto& element : elements(gridView_))
        {
            // the element-wise index sets for finite volume geometry
            std::vector<IndexType> scvIndexSet;
            std::vector<IndexType> scvfIndexSet;

            // fill the element map with seeds
            auto eIdx = problem.elementMapper().index(element);
            elementMap_[eIdx] = element.seed();

            // get the element geometry
            auto elementGeometry = element.geometry();
            const auto& referenceElement = ReferenceElements::general(elementGeometry.type());

            // get the sub control volume geometries of this element
            auto boxGeometries = BoxGeometryHelper<GridView>::getScvAndScvfGeometries(elementGeometry);
            // define some aliases
            auto& scvGeometries = boxGeometries.first;
            auto& scvfGeometries = boxGeometries.second;

            // construct the sub control volumes
            scvIndexSet.reserve(scvGeometries.size());
            scvs_.reserve(scvs_.size() + scvGeometries.size());
            IndexType scvLocalIdx = 0;
            for (auto&& scvGeometry : scvGeometries)
            {
                scvIndexSet.push_back(scvIdx);
                auto dofIdxGlobal = problem.vertexMapper().subIndex(element, dim, scvLocalIdx);
                scvs_.emplace_back(std::make_shared<SubControlVolume>(*scvGeometry,
                                                                      scvIdx++,
                                                                      eIdx,
                                                                      scvLocalIdx++,
                                                                      dofIdxGlobal));
            }

            // construct the sub control volume faces
            scvfIndexSet.reserve(scvfGeometries.size());
            scvfs_.reserve(scvfs_.size() + scvfGeometries.size());
            IndexType scvfLocalIdx = 0;
            for (auto&& scvfGeometry : scvfGeometries)
            {
                // add this scvf to the element's scvf list
                scvfIndexSet.push_back(scvfIdx);

                // find the global scv indices this scvf is belonging to
                std::vector<IndexType> scvIndices(2);
                scvIndices[0] = scvIndexSet[referenceElement.subEntity(scvfLocalIdx, dim-1, 0, dim)];
                scvIndices[1] = scvIndexSet[referenceElement.subEntity(scvfLocalIdx, dim-1, 1, dim)];

                // compute the scvf normal unit outer normal
                auto normal = BoxGeometryHelper<GridView>::normal(elementGeometry, *scvfGeometry);

                const auto v = elementGeometry.corner(scvIndices[1]) - elementGeometry.corner(scvIndices[0]);
                const auto s = v*normal;
                if (std::signbit(s))
                    normal *= -1;

                scvfs_.emplace_back(std::make_shared<SubControlVolumeFace>(*scvfGeometry,
                                                                           scvfGeometry->center(),
                                                                           normal,
                                                                           scvfIdx++,
                                                                           scvIndices,
                                                                           false));
            }

            // construct the sub control volume faces on the domain boundary
            for (const auto& intersection : intersections(gridView_, element))
            {
                scvfLocalIdx = 0;
                if (intersection.boundary())
                {
                    auto isGeometry = intersection.geometry();
                    auto boundaryScvfGeometries = BoxGeometryHelper<GridView>::getBoundaryScvfGeometries(isGeometry);

                    for (auto&& scvfGeometry : scvfGeometries)
                    {
                        // add this scvf to the element's scvf list
                        scvfIndexSet.push_back(scvfIdx);

                        // find the scvs this scvf is belonging to
                        std::vector<IndexType> scvIndices = {static_cast<IndexType>(scvIndexSet[referenceElement.subEntity(intersection.indexInInside(), 1, scvfLocalIdx++, dim)])};

                        // get the unit outer normal through the intersection
                        auto normal = intersection.unitOuterNormal(isGeometry.local(scvfGeometry->center()));

                        scvfs_.emplace_back(std::make_shared<SubControlVolumeFace>(*scvfGeometry,
                                                                                   scvfGeometry->center(),
                                                                                   normal,
                                                                                   scvfIdx++,
                                                                                   scvIndices,
                                                                                   true));
                    }
                }
            }

            // Compute the finite volume element geometries
            fvGeometries_.push_back(FVElementGeometry(*this, scvIndexSet, scvfIndexSet));
        }
    }

private:
    GridView gridView_;
    Dumux::ElementMap<GridView> elementMap_;
    std::vector<std::shared_ptr<SubControlVolume>> scvs_;
    std::vector<std::shared_ptr<SubControlVolumeFace>> scvfs_;
    std::vector<FVElementGeometry> fvGeometries_;
    const FeCache feCache_;
};

} // end namespace

#endif
