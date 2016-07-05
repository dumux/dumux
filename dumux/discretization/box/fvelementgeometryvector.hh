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
#ifndef DUMUX_DISCRETIZATION_BOX_FV_GEOMETRY_VECTOR_HH
#define DUMUX_DISCRETIZATION_BOX_FV_GEOMETRY_VECTOR_HH

#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/implicit/box/properties.hh>
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
    const FVElementGeometry& fvGeometry(const Element& element) const
    {
        return fvGeometries_[problem_().elementMapper().index(element)];
    }

    //! Get a sub control volume with a global scv index
    const SubControlVolume& subControlVolume(IndexType scvIdx) const
    {
        return *scvs_[globalScvIndices_[eIdxBound_][scvIdx]];
    }

    //! Get a sub control volume face with a global scvf index
    const SubControlVolumeFace& subControlVolumeFace(IndexType scvfIdx) const
    {
        return *scvfs_[globalScvfIndices_[eIdxBound_][scvfIdx]];
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

    //! The total number of boundary sub control volume faces
    //! For compatibility reasons with cc methods
    std::size_t numBoundaryScvf() const
    {
        return 0;
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
        problemPtr_ = &problem;

        scvs_.clear();
        scvfs_.clear();
        fvGeometries_.clear();
        elementMap_.clear();

        IndexType numElements = gridView_.size(0);
        IndexType numScvs = 0;
        IndexType numScvfs = 0;
        for (const auto& element : elements(gridView_))
        {
            numScvs += element.subEntities(dim);
            numScvfs += element.subEntities(dim-1);

            // add number of boundary faces
            for (const auto& is : intersections(gridView_, element))
            {
                if (is.boundary())
                    numScvfs += is.geometry().corners();
            }
        }

        scvs_.reserve(numScvs);
        scvfs_.reserve(numScvfs);
        globalScvIndices_.resize(numElements);
        globalScvfIndices_.resize(numElements);
        elementMap_.resize(numElements);

        // Build the SCV and SCV faces
        IndexType scvIdx = 0;
        IndexType scvfIdx = 0;
        for (const auto& element : elements(gridView_))
        {
            // the element-wise index sets for finite volume geometry
            std::vector<IndexType> scvLocalIndexSet;
            std::vector<IndexType> scvfLocalIndexSet;

            // fill the element map with seeds
            auto eIdx = problem.elementMapper().index(element);
            elementMap_[eIdx] = element.seed();

            // get the element geometry
            auto elementGeometry = element.geometry();
            const auto& referenceElement = ReferenceElements::general(elementGeometry.type());

            // get the sub control volume geometries of this element
            auto boxGeometries = std::move(BoxGeometryHelper<GridView>::getScvAndScvfGeometries(elementGeometry));
            // define some aliases
            auto& scvGeometries = boxGeometries.first;
            auto& scvfGeometries = boxGeometries.second;

            // construct the sub control volumes
            globalScvIndices_[eIdx].resize(scvGeometries.size());
            scvLocalIndexSet.reserve(scvGeometries.size());
            IndexType scvLocalIdx = 0;
            for (auto&& scvGeometry : scvGeometries)
            {
                scvLocalIndexSet.push_back(scvLocalIdx);
                auto dofIdxGlobal = problem.vertexMapper().subIndex(element, scvLocalIdx, dim);
                scvs_.emplace_back(std::make_shared<SubControlVolume>(*scvGeometry,
                                                                      scvIdx,
                                                                      eIdx,
                                                                      scvLocalIdx,
                                                                      dofIdxGlobal));

                globalScvIndices_[eIdx][scvLocalIdx] = scvIdx;

                // increment local counters
                scvIdx++; scvLocalIdx++;
            }

            // construct the sub control volume faces
            globalScvfIndices_[eIdx].resize(scvfGeometries.size());
            scvfLocalIndexSet.reserve(scvfGeometries.size());
            IndexType scvfLocalIdx = 0;
            for (auto&& scvfGeometry : scvfGeometries)
            {
                // add this scvf to the element's scvf list
                scvfLocalIndexSet.push_back(scvfLocalIdx);

                // find the global and local scv indices this scvf is belonging to
                std::vector<IndexType> localScvIndices(2);
                localScvIndices[0] = referenceElement.subEntity(scvfLocalIdx, dim-1, 0, dim);
                localScvIndices[1] = referenceElement.subEntity(scvfLocalIdx, dim-1, 1, dim);

                // compute the scvf normal unit outer normal
                auto normal = std::move(BoxGeometryHelper<GridView>::normal(elementGeometry, *scvfGeometry));
                const auto v = elementGeometry.corner(localScvIndices[1]) - elementGeometry.corner(localScvIndices[0]);
                const auto s = v*normal;
                if (std::signbit(s))
                    normal *= -1;

                scvfs_.emplace_back(std::make_shared<SubControlVolumeFace>(*scvfGeometry,
                                                                           scvfGeometry->center(),
                                                                           normal,
                                                                           scvfIdx,
                                                                           scvfLocalIdx,
                                                                           localScvIndices,
                                                                           false));

                globalScvfIndices_[eIdx][scvfLocalIdx] = scvfIdx;

                // increment local counters
                scvfIdx++; scvfLocalIdx++;
            }

            // construct the sub control volume faces on the domain boundary
            for (const auto& intersection : intersections(gridView_, element))
            {
                if (intersection.boundary())
                {
                    auto isGeometry = intersection.geometry();
                    auto boundaryScvfGeometries = std::move(BoxGeometryHelper<GridView>::getBoundaryScvfGeometries(isGeometry));

                    IndexType isScvfLocalIdx = 0;
                    for (auto&& scvfGeometry : boundaryScvfGeometries)
                    {
                        // add this scvf to the element's scvf list
                        scvfLocalIndexSet.push_back(scvfLocalIdx);

                        // find the scvs this scvf is belonging to
                        std::vector<IndexType> localScvIndices = {static_cast<IndexType>(referenceElement.subEntity(intersection.indexInInside(), 1, isScvfLocalIdx, dim))};

                        // get the unit outer normal through the intersection
                        auto normal = intersection.unitOuterNormal(isGeometry.local(scvfGeometry->center()));

                        scvfs_.emplace_back(std::make_shared<SubControlVolumeFace>(*scvfGeometry,
                                                                                   scvfGeometry->center(),
                                                                                   normal,
                                                                                   scvfIdx,
                                                                                   scvfLocalIdx++,
                                                                                   localScvIndices,
                                                                                   true));

                        globalScvfIndices_[eIdx].push_back(scvfIdx++);
                    }
                }
            }

            // Compute the finite volume element geometries
            fvGeometries_.push_back(FVElementGeometry(*this, scvLocalIndexSet, scvfLocalIndexSet));
        }
    }

    // Binding of an element, has to be called before using the fvgeometries
    // Prepares all the volume variables within the element
    // For compatibility reasons with the FVGeometry cache being disabled
    void bindElement(const Element& element)
    {
        eIdxBound_ = problem_().elementMapper().index(element);
    }
    // this function is for compatibility reasons with cc methods
    void bind(const Element& element)
    {
        bindElement(element);
    }

private:
    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;

    IndexType eIdxBound_;
    GridView gridView_;
    Dumux::ElementMap<GridView> elementMap_;
    std::vector<std::shared_ptr<SubControlVolume>> scvs_;
    std::vector<std::shared_ptr<SubControlVolumeFace>> scvfs_;
    std::vector<FVElementGeometry> fvGeometries_;

    std::vector<std::vector<IndexType>> globalScvIndices_;
    std::vector<std::vector<IndexType>> globalScvfIndices_;
    const FeCache feCache_;
};

// specialization in case the FVElementGeometries are not stored
template<class TypeTag>
class BoxFVElementGeometryVector<TypeTag, false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CoordScalar = typename GridView::ctype;

    using FeCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    using FeLocalBasis = typename FeCache::FiniteElementType::Traits::LocalBasisType;
    using ReferenceElements = typename Dune::ReferenceElements<CoordScalar, dim>;

public:
    //! Constructor
    BoxFVElementGeometryVector(const GridView gridView)
    : gridView_(gridView), elementMap_(gridView), eIdxBound_(-1) {}

    /* \brief Get the finite volume geometry of an element
     * \note The finite volume geometry offers iterators over the sub control volumes
     *       and the sub control volume faces of an element.
     */
    const FVElementGeometry& fvGeometry(const Element& element) const
    {
        auto eIdx = problem_().elementMapper().index(element);
        if (eIdx != eIdxBound_)
            DUNE_THROW(Dune::InvalidStateException, "Trying to get FVElementGeometry for element with index " << eIdx
                                                    << ", but bound is element " << eIdxBound_ << ". Please call bindElement() before using the fvGeometry.");
        return fvGeometry_[0];
    }

    //! Get a sub control volume with a global scv index
    const SubControlVolume& subControlVolume(IndexType scvIdx) const
    {
        return stencilScvs_[scvIdx];
    }

    //! Get a sub control volume face with a global scvf index
    const SubControlVolumeFace& subControlVolumeFace(IndexType scvfIdx) const
    {
        return stencilScvfs_[scvfIdx];
    }

    //! Get a local finite element basis
    const FeLocalBasis& feLocalBasis(const Dune::GeometryType& type) const
    {
        return feCache_.get(type).localBasis();
    }

    //! The total number of sub control volumes
    std::size_t numScv() const
    {
        return numScvs_;
    }

    //! The total number of sun control volume faces
    std::size_t numScvf() const
    {
        return numScvf_;
    }

    //! The total number of boundary sub control volume faces
    //! For compatibility reasons with cc methods
    std::size_t numBoundaryScvf() const
    {
        return 0;
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
        problemPtr_ = &problem;
        elementMap_.clear();
        release_();
        eIdxBound_ = -1;

        auto numElems = gridView_.size(0);
        elementMap_.resize(numElems);
        scvIndices_.resize(numElems);
        scvFaceIndices_.resize(numElems);

        // save data on the grid's scvs and scvfs
        IndexType scvIdx = 0;
        IndexType scvfIdx = 0;
        for (const auto& element : elements(gridView_))
        {
            // the element-wise index sets for finite volume geometry
            std::vector<IndexType> scvIndexSet;
            std::vector<IndexType> scvfIndexSet;

            // fill the element map with seeds
            auto eIdx = problem.elementMapper().index(element);
            elementMap_[eIdx] = element.seed();

            // store scv indices of the element
            for (int scvLocalIdx = 0; scvLocalIdx < element.subEntities(dim); ++scvLocalIdx)
                scvIndexSet.push_back(scvIdx++);

            // store scv face indices of the element
            for (int scvfLocalIdx = 0; scvfLocalIdx < element.subEntities(dim-1); ++scvfLocalIdx)
                scvfIndexSet.push_back(scvfIdx++);

            // store the sub control volume face indices on the domain boundary
            for (const auto& intersection : intersections(gridView_, element))
            {
                if (intersection.boundary())
                {
                    for (int localNodeIdx = 0; localNodeIdx < intersection.geometry().corners(); ++localNodeIdx)
                    {
                        // add this scvf to the element's scvf list
                        scvfIndexSet.push_back(scvfIdx++);
                    }
                }
            }

            // store the sets of indices in the data containers
            scvIndices_[eIdx] = scvIndexSet;
            scvFaceIndices_[eIdx] = scvfIndexSet;
        }

        numScvs_ = scvIdx;
        numScvf_ = scvfIdx;
    }

    // build the geometries inside an element
    void bindElement(const Element& element)
    {
        eIdxBound_ = problem_().elementMapper().index(element);
        release_();
        makeElementGeometries_(element);
    }

    // for cc methods, this method binds all the elements in the stencil.
    // Here, we simply forward to the bindElement() routine.
    void bind(const Element& element)
    {
        bindElement(element);
    }

private:
    void release_()
    {
        stencilScvs_.clear();
        stencilScvfs_.clear();
        fvGeometry_.clear();
    }

    void makeElementGeometries_(const Element& element)
    {
        const auto eIdx = problem_().elementMapper().index(element);

        // the element-wise index sets for finite volume geometry
        auto scvIndexSet = scvIndices_[eIdx];
        auto scvfIndexSet = scvFaceIndices_[eIdx];
        std::vector<IndexType> localScvIndexSet;
        std::vector<IndexType> localScvfIndexSet;
        localScvIndexSet.reserve(scvIndexSet.size());
        localScvfIndexSet.reserve(scvfIndexSet.size());

        // get the element geometry
        auto elementGeometry = element.geometry();
        const auto& referenceElement = ReferenceElements::general(elementGeometry.type());

        // get the sub control volume geometries of this element
        auto boxGeometries = std::move(BoxGeometryHelper<GridView>::getScvAndScvfGeometries(elementGeometry));
        // define some aliases
        auto& scvGeometries = boxGeometries.first;
        auto& scvfGeometries = boxGeometries.second;

        // construct the sub control volumes
        IndexType scvLocalIdx = 0;
        for (auto&& scvGeometry : scvGeometries)
        {
            localScvIndexSet.push_back(scvLocalIdx);
            auto dofIdxGlobal = problem_().vertexMapper().subIndex(element, scvLocalIdx, dim);

            // add scv to the local container
            stencilScvs_.emplace_back(SubControlVolume(*scvGeometry,
                                                       scvIndexSet[scvLocalIdx],
                                                       eIdx,
                                                       scvLocalIdx,
                                                       dofIdxGlobal));
            // increment local counter
            scvLocalIdx++;
        }

        // construct the sub control volume faces
        IndexType scvfLocalIdx = 0;
        for (auto&& scvfGeometry : scvfGeometries)
        {
            localScvfIndexSet.push_back(scvfLocalIdx);

            // find the local scv indices this scvf is connsected to
            std::vector<IndexType> localScvIndices(2);
            localScvIndices[0] = referenceElement.subEntity(scvfLocalIdx, dim-1, 0, dim);
            localScvIndices[1] = referenceElement.subEntity(scvfLocalIdx, dim-1, 1, dim);

            // compute the scvf normal unit outer normal
            auto normal = std::move(BoxGeometryHelper<GridView>::normal(elementGeometry, *scvfGeometry));
            const auto v = elementGeometry.corner(localScvIndices[1]) - elementGeometry.corner(localScvIndices[0]);
            const auto s = v*normal;
            if (std::signbit(s))
                normal *= -1;

            stencilScvfs_.emplace_back(SubControlVolumeFace(*scvfGeometry,
                                                            scvfGeometry->center(),
                                                            normal,
                                                            scvfIndexSet[scvfLocalIdx],
                                                            scvfLocalIdx,
                                                            localScvIndices,
                                                            false));

            // increment local counter
            scvfLocalIdx++;
        }

        // construct the sub control volume faces on the domain boundary
        for (const auto& intersection : intersections(gridView_, element))
        {
            if (intersection.boundary())
            {
                auto isGeometry = intersection.geometry();
                auto boundaryScvfGeometries = std::move(BoxGeometryHelper<GridView>::getBoundaryScvfGeometries(isGeometry));

                IndexType boundaryFaceLocalIdx = 0;
                for (auto&& scvfGeometry : boundaryScvfGeometries)
                {
                    localScvfIndexSet.push_back(scvfLocalIdx);

                    // find the scv this scvf is connected to
                    std::vector<IndexType> localScvIndices = {static_cast<IndexType>(referenceElement.subEntity(intersection.indexInInside(), 1, boundaryFaceLocalIdx++, dim))};

                    // get the unit outer normal through the intersection
                    auto normal = intersection.unitOuterNormal(isGeometry.local(scvfGeometry->center()));

                    stencilScvfs_.emplace_back(SubControlVolumeFace(*scvfGeometry,
                                                                    scvfGeometry->center(),
                                                                    normal,
                                                                    scvfIndexSet[scvfLocalIdx],
                                                                    scvfLocalIdx,
                                                                    localScvIndices,
                                                                    true));

                    // increment local counter
                    scvfLocalIdx++;
                }
            }
        }

        // Compute the finite volume element geometries
        fvGeometry_.push_back(FVElementGeometry(*this, localScvIndexSet, localScvfIndexSet));
    }

    const Problem& problem_() const
    { return *problemPtr_; }

    GridView gridView_;
    const Problem* problemPtr_;
    const FeCache feCache_;

    // Information on the global number of geometries
    IndexType numScvs_;
    IndexType numScvf_;

    // vectors that store the global data
    Dumux::ElementMap<GridView> elementMap_;
    std::vector<std::vector<IndexType>> scvIndices_;
    std::vector<std::vector<IndexType>> scvFaceIndices_;

    // vectors to store the geometries temporarily after binding an element
    IndexType eIdxBound_;
    std::vector<SubControlVolume> stencilScvs_;
    std::vector<SubControlVolumeFace> stencilScvfs_;
    std::vector<FVElementGeometry> fvGeometry_;
};

} // end namespace

#endif
