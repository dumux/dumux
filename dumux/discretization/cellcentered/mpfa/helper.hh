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
 * \brief Helper class to get the required information on an interaction volume.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_HELPER_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_HELPER_HH

#include "methods.hh"
#include "facetypes.hh"

namespace Dumux
{

// Mpfa method-specific implementation of the helper class (again dimension-dependent)
template<class TypeTag, MpfaMethods Method, int dim, int dimWorld>
class MpfaMethodHelper;

// dimension-specific implementation of the helper class (common for all methods)
template<class TypeTag, int dim, int dimWorld>
class MpfaDimensionHelper;

// Specialization for dim == 2, dimWorld == 2
template<class TypeTag>
class MpfaDimensionHelper<TypeTag, /*dim*/2, /*dimWorld*/2>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using LocalIndexType = typename InteractionVolume::Traits::LocalIndexType;

    // We know that dim = 2 and dimworld = 2, but
    // the dim = 2, dimWorld = 3 specialization forwards to some methods of this class
    // Thus, we need to get the GlobalPosition vector etc right
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using ScvfVector = std::array<const SubControlVolumeFace*, dim>;
    using LocalBasis = std::array<GlobalPosition, dim>;

    using PointVector = std::vector<GlobalPosition>;
    using FaceReferenceElements = typename Dune::ReferenceElements<Scalar, dim-1>;

public:
    // Gets the two scv faces in the outer element, that share the vertex
    // orders them to form a right hand system, local indices can be deduced from on the rotational direction
    static ScvfVector getCommonAndNextScvFace(const SubControlVolumeFace& outsideScvf,
                                              const FVElementGeometry& fvGeometry,
                                              const bool clockWise)
    {
        LocalIndexType commonFaceIdx = clockWise ? 0 : 1;
        LocalIndexType nextFaceIdx = clockWise ? 1 : 0;
        auto vIdxGlobal = outsideScvf.vertexIndex();

        unsigned int count = 0;
        ScvfVector scvfVector({nullptr});
        for (auto&& scvf : scvfs(fvGeometry))
        {
            if (scvf.vertexIndex() == vIdxGlobal)
            {
                if (scvf.outsideScvIdx() == outsideScvf.insideScvIdx())
                {
                    scvfVector[commonFaceIdx] = &scvf;
                    count++;
                }
                else
                {
                    scvfVector[nextFaceIdx] = &scvf;
                    count++;
                }
            }
        }

        // make sure we found two faces
        assert(count == 2 && "did not find two scv faces sharing the vertex in the outside element");
        return scvfVector;
    }

    // calculates the inner normal vectors
    static LocalBasis calculateInnerNormals(const LocalBasis& localBasis)
    {
        static const Dune::FieldMatrix<Scalar, dim, dim> R = {{0.0, 1.0}, {-1.0, 0.0}};
        // make sure the basis forms a right hand system
        assert(isRightHandSystem(localBasis) > 0 && "Local basis does not form a right hand system");

        LocalBasis innerNormals;
        R.mv(localBasis[1], innerNormals[0]);
        R.mv(localBasis[0], innerNormals[1]);
        innerNormals[1] *= -1;

        return innerNormals;
    }

    // the determinant of the local basis is equal to the cross product for dim = dimWorld = 2
    static Scalar calculateDetX(const LocalBasis& localBasis)
    {
        // make sure the basis forms a right hand system
        assert(isRightHandSystem(localBasis) > 0 && "Local basis does not form a right hand system");
        return crossProduct<Scalar>(localBasis[0], localBasis[1]);
    }

    // returns the global number of scvfs in the grid
    static std::size_t getGlobalNumScvf(const GridView& gridView)
    {
        Dune::GeometryType triangle, quadrilateral;
        triangle.makeTriangle();
        quadrilateral.makeQuadrilateral();

        return gridView.size(triangle)*6 + gridView.size(quadrilateral)*8;
    }

    //! Check whether or not the local basis forms a right hand system
    static bool isRightHandSystem(const LocalBasis& localBasis)
    { return !std::signbit(crossProduct<Scalar>(localBasis[0], localBasis[1])); }

    //! get sub control volume face corners on a given face geometry for the given local index
    static PointVector getScvfCorners(const PointVector& isCorners, unsigned int cornerIdx)
    {
        // compute intersection center (in 2d intersections have two corners)
        GlobalPosition center = isCorners[0] + isCorners[1];
        center /= 2;

        if (cornerIdx == 0)
            return PointVector({center, isCorners[0]});
        else if (cornerIdx == 1)
            return PointVector({center, isCorners[1]});
        else
            DUNE_THROW(Dune::InvalidStateException, "local index exceeds the number of corners of 2d intersections");
    }

    //! calculate integration point on an scvf
    static GlobalPosition getScvfIntegrationPoint(const PointVector& scvfCorners, Scalar q)
    {
        // ordering -> first corner: facet center, last corner: vertex
        if (q == 0.0)
            return scvfCorners[0];

        auto d = scvfCorners[1] - scvfCorners[0];
        d *= q;
        return scvfCorners[0] + d;
    }

    //! calculate the area of a scvf
    static Scalar getScvfArea(const PointVector& scvfCorners)
    { return (scvfCorners[1]-scvfCorners[0]).two_norm(); }

    static std::size_t getNumLocalScvfs(const Dune::GeometryType gt)
    {
        if (gt == Dune::GeometryType(Dune::GeometryType::simplex, 2))
            return 6;
        else if (gt == Dune::GeometryType(Dune::GeometryType::cube, 2))
            return 8;
        else
            DUNE_THROW(Dune::InvalidStateException, "unknown 2d geometry type " << gt);
    }
};

// Specialization for dim == 2, dimWorld == 3
template<class TypeTag>
class MpfaDimensionHelper<TypeTag, /*dim*/2, /*dimWorld*/3>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using LocalIndexType = typename InteractionVolume::Traits::LocalIndexType;

    // We know that dim = 2 and dimworld = 3
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using ScvfVector = std::array<const SubControlVolumeFace*, dim>;
    using LocalBasis = std::array<GlobalPosition, dim>;
    using PointVector = std::vector<GlobalPosition>;

public:
    // calculates the inner normal vectors
    static LocalBasis calculateInnerNormals(const LocalBasis& localBasis)
    {
        // make sure the basis forms a right hand system
        assert(isRightHandSystem(localBasis) > 0 && "Local basis does not form a right hand system");

        // compute vector normal to the basis plane
        auto normal = crossProduct<Scalar>(localBasis[0], localBasis[1]);
        normal /= normal.two_norm();

        // compute inner normals using the normal vector
        LocalBasis innerNormals;
        innerNormals[0] = crossProduct<Scalar>(localBasis[1], normal);
        innerNormals[1] = crossProduct<Scalar>(normal, localBasis[0]);

        return innerNormals;
    }

    // This is actually not the determinant of the basis for dim = 2 < dimWorld = 3
    // It simply is the area of the parallelogram spanned by the basis vectors
    static Scalar calculateDetX(const LocalBasis& localBasis)
    { return std::abs( crossProduct<Scalar>(localBasis[0], localBasis[1]).two_norm() ); }

    // returns the global number of scvfs in the grid
    static std::size_t getGlobalNumScvf(const GridView& gridView)
    {
        Dune::GeometryType triangle, quadrilateral;
        triangle.makeTriangle();
        quadrilateral.makeQuadrilateral();

        return gridView.size(triangle)*6 + gridView.size(quadrilateral)*8;
    }

    //! For 2d in 3d there is no unique basis forming a right hand system
    static bool isRightHandSystem(const LocalBasis& localBasis)
    { return true; }

    //! get sub control volume face corners on a given face geometry for the given local index
    static PointVector getScvfCorners(const PointVector& isCorners, unsigned int cornerIdx)
    { return MpfaDimensionHelper<TypeTag, dim, dim>::getScvfCorners(isCorners, cornerIdx); }

    //! calculate integration point on an scvf
    static GlobalPosition getScvfIntegrationPoint(const PointVector& scvfCorners, Scalar q)
    { return MpfaDimensionHelper<TypeTag, dim, dim>::getScvfIntegrationPoint(scvfCorners, q); }

    //! calculate the area of a scvf
    static Scalar getScvfArea(const PointVector& scvfCorners)
    { return MpfaDimensionHelper<TypeTag, dim, dim>::getScvfArea(scvfCorners); }

    static std::size_t getNumLocalScvfs(const Dune::GeometryType gt)
    { return MpfaDimensionHelper<TypeTag, dim, dim>::getNumLocalScvfs(gt); }

    static ScvfVector getCommonAndNextScvFace(const SubControlVolumeFace& outsideScvf,
                                              const FVElementGeometry& fvGeometry,
                                              const bool clockWise)
    { return MpfaDimensionHelper<TypeTag, dim, dim>::getCommonAndNextScvFace(outsideScvf, fvGeometry, clockWise); }
};

// Specialization for dim == 3 (dimWorld has to be 3)
template<class TypeTag>
class MpfaDimensionHelper<TypeTag, /*dim*/3, /*dimWorld*/3>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using LocalIndexType = typename InteractionVolume::Traits::LocalIndexType;

    // We know that dim = 3 and dimworld = 3
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using ScvfVector = std::array<const SubControlVolumeFace*, dim>;
    using LocalBasis = std::array<GlobalPosition, dim>;

    using PointVector = std::vector<GlobalPosition>;
    using FaceReferenceElements = typename Dune::ReferenceElements<Scalar, dim-1>;

public:

    // calculates the inner normal vectors
    static LocalBasis calculateInnerNormals(const LocalBasis& localBasis)
    {
        LocalBasis innerNormals;

        innerNormals[0] = crossProduct<Scalar>(localBasis[1], localBasis[2]);
        innerNormals[1] = crossProduct<Scalar>(localBasis[2], localBasis[0]);
        innerNormals[2] = crossProduct<Scalar>(localBasis[0], localBasis[1]);

        return innerNormals;
    }

    // calculates the determinant of the local basis
    static Scalar calculateDetX(const LocalBasis& localBasis)
    {
        assert(isRightHandSystem(localBasis) && "Local basis does not form a right hand system");
        return tripleProduct<Scalar>(localBasis[0], localBasis[1], localBasis[2]);
    }

    // returns the global number of scvfs in the grid
    static std::size_t getGlobalNumScvf(const GridView& gridView)
    {
        Dune::GeometryType simplex, pyramid, prism, cube;
        simplex.makeTetrahedron();
        pyramid.makePyramid();
        prism.makePrism();
        cube.makeHexahedron();

        return gridView.size(simplex)*12 + gridView.size(pyramid)*16 + gridView.size(prism)*18 + gridView.size(cube)*24;
    }

    //! Check whether or not the local basis forms a right hand system
    static bool isRightHandSystem(const LocalBasis& localBasis)
    { return !std::signbit(tripleProduct<Scalar>(localBasis[0], localBasis[1], localBasis[2])); }

    //! get sub control volume face corners on a given face geometry for the given local index
    static PointVector getScvfCorners(const PointVector& isCorners, unsigned int cornerIdx)
    {
        // maximum number of necessary points is 9 (for quadrilateral)
        GlobalPosition p[9];
        auto numCorners = isCorners.size();

        // fill point vector with center and corners
        p[0] = 0.0;
        for (int i = 0; i < numCorners; ++i)
        {
            p[0] += isCorners[i];
            p[i+1] = isCorners[i];
        }
        p[0] /= numCorners;

        // proceed according to number of corners
        switch (numCorners)
        {
        case 3: // triangle
        {
            // add edge midpoints, there are 3 for triangles
            p[numCorners+1] = isCorners[1] + isCorners[0];
            p[numCorners+1] /= 2;
            p[numCorners+2] = isCorners[2] + isCorners[0];
            p[numCorners+2] /= 2;
            p[numCorners+3] = isCorners[2] + isCorners[1];
            p[numCorners+3] /= 2;

            //! Only build the maps the first time we encounter a triangle
            static const std::uint8_t vo = 1; //! vertex offset in point vector p
            static const std::uint8_t eo = 4; //! edge offset in point vector p
            static const std::uint8_t map[3][4] =
            {
                {0, eo+1, eo+0, vo+0},
                {0, eo+0, eo+2, vo+1},
                {0, eo+2, eo+1, vo+2}
            };

            return PointVector( {p[map[cornerIdx][0]],
                                 p[map[cornerIdx][1]],
                                 p[map[cornerIdx][2]],
                                 p[map[cornerIdx][3]]} );
        }
        case 4: // quadrilateral
        {
            // add edge midpoints, there are 4 for quadrilaterals
            p[numCorners+1] = isCorners[2] + isCorners[0];
            p[numCorners+1] /= 2;
            p[numCorners+2] = isCorners[3] + isCorners[1];
            p[numCorners+2] /= 2;
            p[numCorners+3] = isCorners[1] + isCorners[0];
            p[numCorners+3] /= 2;
            p[numCorners+4] = isCorners[3] + isCorners[2];
            p[numCorners+4] /= 2;

            //! Only build the maps the first time we encounter a quadrilateral
            static const std::uint8_t vo = 1; //! vertex offset in point vector p
            static const std::uint8_t eo = 5; //! face offset in point vector p
            static const std::uint8_t map[4][4] =
            {
                {0, eo+0, eo+2, vo+0},
                {0, eo+2, eo+1, vo+1},
                {0, eo+3, eo+0, vo+2},
                {0, eo+1, eo+3, vo+3}
            };

            return PointVector( {p[map[cornerIdx][0]],
                                 p[map[cornerIdx][1]],
                                 p[map[cornerIdx][2]],
                                 p[map[cornerIdx][3]]} );
        }
        default:
            DUNE_THROW(Dune::NotImplemented, "Mpfa scvf corners for dim=" << dim
                                                            << " dimWorld=" << dimWorld
                                                            << " corners=" << numCorners);
        }
    }

    //! calculate integration point on a scvf
    static GlobalPosition getScvfIntegrationPoint(const PointVector& scvfCorners, Scalar q)
    {
        // scvfs in 3d are always quadrilaterals
        // ordering -> first corner: facet center, last corner: vertex
        if (q == 0.0)
            return scvfCorners[0];

        auto d = scvfCorners[3] - scvfCorners[0];
        d *= q;
        return scvfCorners[0] + d;
    }

    //! calculate scvf area
    static Scalar getScvfArea(const PointVector& scvfCorners)
    {
        // after Wolfram alpha quadrilateral area
        return 0.5*Dumux::crossProduct(scvfCorners[3]-scvfCorners[0], scvfCorners[2]-scvfCorners[1]).two_norm();
    }

    //! determine the number of scvfs of an element
    static std::size_t getNumLocalScvfs(const Dune::GeometryType gt)
    {
        if (gt == Dune::GeometryType(Dune::GeometryType::simplex, 3))
            return 12;
        else if (gt == Dune::GeometryType(Dune::GeometryType::pyramid, 3))
            return 16;
        else if (gt == Dune::GeometryType(Dune::GeometryType::prism, 3))
            return 18;
        else if (gt == Dune::GeometryType(Dune::GeometryType::cube, 3))
            return 24;
        else
            DUNE_THROW(Dune::InvalidStateException, "unknown 3d geometry type " << gt);
    }
};

/*!
 * \ingroup Mpfa
 * \brief Helper class to get the required information on an interaction volume.
 *        Specializations depending on the method and dimension are provided.
 */
template<class TypeTag, MpfaMethods Method, int dim, int dimWorld>
class CCMpfaHelperImplementation : public MpfaDimensionHelper<TypeTag, dim, dimWorld>,
                                   public MpfaMethodHelper<TypeTag, Method, dim, dimWorld>
{
    using Implementation = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalIndexType = typename GridView::IndexSet::IndexType;
    using LocalIndexType = typename InteractionVolume::Traits::LocalIndexType;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using ScvfVector = std::array<const SubControlVolumeFace*, dim>;
    using LocalBasis = std::array<GlobalPosition, dim>;

    using CoordScalar = typename GridView::ctype;
    using ReferenceElements = typename Dune::ReferenceElements<CoordScalar, dim>;

    static constexpr bool enableInteriorBoundaries = GET_PROP_VALUE(TypeTag, EnableInteriorBoundaries);

public:
    // returns shared pointers to the scv faces that share a vertex in the order of a right hand system
    static ScvfVector getScvFacesAtVertex(const GlobalIndexType vIdxGlobal,
                                          const Element& element,
                                          const FVElementGeometry& fvGeometry)
    {
        ScvfVector scvfVector({nullptr});
        LocalBasis basisVectors;

        // The element center
        auto elementCenter = element.geometry().center();

        LocalIndexType count = 0;
        for (auto&& scvf : scvfs(fvGeometry))
        {
            if (scvf.vertexIndex() == vIdxGlobal)
            {
                scvfVector[count] = &scvf;
                basisVectors[count] = scvf.ipGlobal();
                basisVectors[count] -= elementCenter;
                count++;
            }
        }

        // We should always find dim faces sharing the vertex
        assert(count == dim);

        // sort the scv faces to form a right hand system
        if (!Implementation::isRightHandSystem(basisVectors))
            std::swap(scvfVector[0], scvfVector[1]);

        return scvfVector;
    }

    // Finds the local index in an ScvfVector that corresponds to the face that shares a facet with the given outsideScvf
    static LocalIndexType getCommonFaceLocalIndex(const SubControlVolumeFace& outsideScvf,
                                                  const ScvfVector& insideScvFaces)
    {
        auto outsideScvIdx = outsideScvf.insideScvIdx();
        for (int i = 0; i < insideScvFaces.size(); i++)
            if (contains(insideScvFaces[i]->outsideScvIndices(), outsideScvIdx))
                return i;

        DUNE_THROW(Dune::InvalidStateException, "Could not find corresponding scvf in the provided vector of scvfs.");
    }

    // Finds the local index in an ScvfVector that corresponds to a given scvf
    static LocalIndexType getLocalFaceIndex(const SubControlVolumeFace& scvf,
                                            const ScvfVector& scvFaces)
    {
        for (int i = 0; i < scvFaces.size(); i++)
            if (scvFaces[i]->index() == scvf.index())
                return i;

        DUNE_THROW(Dune::InvalidStateException, "Could not find corresponding scvf in the provided vector of scvfs.");
    }

    // Returns the MpfaFaceType of an scv face
    static MpfaFaceTypes getMpfaFaceType(const Problem& problem,
                                         const Element& element,
                                         const SubControlVolumeFace& scvf)
    {
        // We do simplified checks if interior boundaries are disabled
        if (!enableInteriorBoundaries)
        {
            if (!scvf.boundary())
                return MpfaFaceTypes::interior;

            const auto bcTypes = problem.boundaryTypes(element, scvf);
            if (bcTypes.hasOnlyNeumann())
                return MpfaFaceTypes::neumann;
            if (bcTypes.hasOnlyDirichlet())
                return MpfaFaceTypes::dirichlet;

            // throw exception
            return throwBoundaryExceptionMessage(bcTypes);
        }
        else
        {
            const auto bcTypes = problem.boundaryTypes(element, scvf);

            // if we are on an interior boundary return interior types
            if (problem.model().globalFvGeometry().isOnInteriorBoundary(scvf))
            {
                if (bcTypes.hasOnlyNeumann())
                    return MpfaFaceTypes::interiorNeumann;
                if (bcTypes.hasOnlyDirichlet())
                    return MpfaFaceTypes::interiorDirichlet;

                // throw exception
                return throwBoundaryExceptionMessage(bcTypes);
            }

            if (scvf.boundary())
            {
                if (bcTypes.hasOnlyNeumann())
                    return MpfaFaceTypes::neumann;
                if (bcTypes.hasOnlyDirichlet())
                    return MpfaFaceTypes::dirichlet;

                // throw exception
                return throwBoundaryExceptionMessage(bcTypes);
            }

            // This is an interior scvf
            return MpfaFaceTypes::interior;
        }
    }

    // returns a vector, which maps a bool (true if ghost vertex) to each vertex index
    static std::vector<bool> findGhostVertices(const Problem& problem, const GridView& gridView)
    {
        std::vector<bool> ghostVertices(gridView.size(dim), false);

        // if not run in parallel, skip the rest
        if (Dune::MPIHelper::getCollectiveCommunication().size() == 1)
            return ghostVertices;

        // mpfa methods can not yet handle ghost cells
        if (gridView.ghostSize(0) > 0)
            DUNE_THROW(Dune::InvalidStateException, "Mpfa methods in parallel do not work with ghost cells. Use overlap cells instead.");

        // mpfa methods have to have overlapping cells
        if (gridView.overlapSize(0) == 0)
            DUNE_THROW(Dune::InvalidStateException, "Grid no overlaping cells. This is required by mpfa methods in parallel.");

        for (const auto& element : elements(gridView))
        {
            for (const auto& is : intersections(gridView, element))
            {
                if (!is.neighbor() && !is.boundary())
                {
                    const auto& refElement = ReferenceElements::general(element.geometry().type());
                    for (unsigned int isVertex = 0; isVertex < is.geometry().corners(); ++isVertex)
                    {
                        const auto vIdxLocal = refElement.subEntity(is.indexInInside(), 1, isVertex, dim);
                        const auto vIdxGlobal = problem.vertexMapper().subIndex(element, vIdxLocal, dim);

                        ghostVertices[vIdxGlobal] = true;
                    }
                }
            }
        }

        return ghostVertices;
    }

    //! returns whether or not a value exists in a vector
    template<typename V1, typename V2>
    static bool contains(const std::vector<V1>& vector, const V2 value)
    { return std::find(vector.begin(), vector.end(), value) != vector.end(); }

    // calculates the product of a transposed vector n, a Matrix M and another vector v (n^T M v)
    static Scalar nT_M_v(const GlobalPosition& n,
                         const DimWorldMatrix& M,
                         const GlobalPosition& v)
    {
        GlobalPosition tmp;
        M.mv(v, tmp);
        return n*tmp;
    }

    // calculates the product of a transposed vector n, a Scalar M and another vector v (n^T M v)
    static Scalar nT_M_v(const GlobalPosition& n,
                         const Scalar M,
                         const GlobalPosition& v)
    {
        return M*(n*v);
    }

private:
    static MpfaFaceTypes throwBoundaryExceptionMessage(const BoundaryTypes& bcTypes)
    {
        // throw for outflow or mixed boundary conditions
        if (bcTypes.hasOutflow())
            DUNE_THROW(Dune::NotImplemented, "outflow BC for mpfa schemes");
        if (bcTypes.hasDirichlet() && bcTypes.hasNeumann())
            DUNE_THROW(Dune::InvalidStateException, "Mixed BC are not allowed for cellcentered schemes");

        DUNE_THROW(Dune::InvalidStateException, "unknown boundary condition type");
    }
};

/*!
 * \ingroup Mpfa
 * \brief Helper class for the mpfa methods for the construction of the interaction regions etc.
 *        It inherits from dimension-, dimensionworld- and method-specific implementations.
 */
template<class TypeTag>
using CCMpfaHelper = CCMpfaHelperImplementation<TypeTag,
                                                GET_PROP_VALUE(TypeTag, MpfaMethod),
                                                GET_PROP_TYPE(TypeTag, GridView)::dimension,
                                                GET_PROP_TYPE(TypeTag, GridView)::dimensionworld>;

} // end namespace

// The implemented helper classes need to be included here
#include <dumux/discretization/cellcentered/mpfa/lmethod/helper.hh>
#include <dumux/discretization/cellcentered/mpfa/omethod/helper.hh>
#include <dumux/discretization/cellcentered/mpfa/omethodfps/helper.hh>

#endif
