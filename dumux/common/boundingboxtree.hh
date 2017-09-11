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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief An axis-aligned bounding box volume hierarchy for dune grids
 *
 * Dumux implementation of an AABB tree
 * adapted from implementation in FEniCS by Anders Logg
 */
#ifndef DUMUX_BOUNDINGBOXTREE_HH
#define DUMUX_BOUNDINGBOXTREE_HH

#include <algorithm>
#include <memory>
#include <type_traits>
#include <dune/geometry/referenceelements.hh>
#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>
#include <dumux/common/math.hh>
#include <dumux/common/geometrycollision.hh>

namespace Dumux {

/*!
 * \brief A helper class to optimize dimension-dependent methods
 */
template <int dimworld>
class BoundingBoxTreeHelper {};

//! Helper methods for world dimension three
template <>
class BoundingBoxTreeHelper<3>
{
    // An epsilon for floating point operations
    static constexpr double eps_ = 1.0e-7;
    typedef Dune::FieldVector<double, 3> GlobalPosition;

public:
    //! Check whether a point is inside a given three-dimensional geometry
    template <class Geometry>
    static typename std::enable_if<Geometry::mydimension == 3, bool>::type
    pointInGeometry(const Geometry& geometry, const GlobalPosition& point)
    {
        // get the geometry type
        auto type = geometry.type();

        // if it's a tetrahedron we can check directly
        if (type.isTetrahedron())
        {
            return pointInTetrahedron(geometry.corner(0),
                                      geometry.corner(1),
                                      geometry.corner(2),
                                      geometry.corner(3),
                                      point);
        }
        // split hexahedrons into five tetrahedrons
        else if (type.isHexahedron())
        {
            if (pointInTetrahedron(geometry.corner(0),
                                   geometry.corner(1),
                                   geometry.corner(3),
                                   geometry.corner(5),
                                   point))
                return true;
            if (pointInTetrahedron(geometry.corner(0),
                                   geometry.corner(5),
                                   geometry.corner(6),
                                   geometry.corner(4),
                                   point))
                return true;
            if (pointInTetrahedron(geometry.corner(5),
                                   geometry.corner(3),
                                   geometry.corner(6),
                                   geometry.corner(7),
                                   point))
                return true;
            if (pointInTetrahedron(geometry.corner(0),
                                   geometry.corner(3),
                                   geometry.corner(2),
                                   geometry.corner(6),
                                   point))
                return true;
            if (pointInTetrahedron(geometry.corner(5),
                                   geometry.corner(3),
                                   geometry.corner(0),
                                   geometry.corner(6),
                                   point))
                return true;
            return false;
        }
        // split pyramids into two tetrahedrons
        else if (type.isPyramid())
        {
            if (pointInTetrahedron(geometry.corner(0),
                                   geometry.corner(1),
                                   geometry.corner(2),
                                   geometry.corner(4),
                                   point))
                return true;
            if (pointInTetrahedron(geometry.corner(1),
                                   geometry.corner(3),
                                   geometry.corner(2),
                                   geometry.corner(4),
                                   point))
                return true;
            return false;
        }
        // split prisms into three tetrahedrons
        else if (type.isPrism())
        {
            if (pointInTetrahedron(geometry.corner(0),
                                   geometry.corner(1),
                                   geometry.corner(2),
                                   geometry.corner(4),
                                   point))
                return true;
            if (pointInTetrahedron(geometry.corner(3),
                                   geometry.corner(0),
                                   geometry.corner(2),
                                   geometry.corner(4),
                                   point))
                return true;
            if (pointInTetrahedron(geometry.corner(2),
                                   geometry.corner(5),
                                   geometry.corner(3),
                                   geometry.corner(4),
                                   point))
                return true;
            return false;
        }
        else
            DUNE_THROW(Dune::NotImplemented,
                        "Intersection for point and geometry "
                        << type << " in three-dimensional world.");
    }

    //! Check whether a point is inside a given two-dimensional geometry
    template <class Geometry>
    static typename std::enable_if<Geometry::mydimension == 2, bool>::type
    pointInGeometry(const Geometry& geometry, const GlobalPosition& point)
    {
        // get the geometry type
        auto type = geometry.type();

        // if it's a triangle we can check directly
        if (type.isTriangle())
        {
            return pointInTriangle(geometry.corner(0),
                                   geometry.corner(1),
                                   geometry.corner(2),
                                   point);
        }
        // split quadrilaterals into two triangles
        else if (type.isQuadrilateral())
        {
            if (pointInTriangle(geometry.corner(0),
                                geometry.corner(1),
                                geometry.corner(3),
                                point))
                return true;
            if (pointInTriangle(geometry.corner(0),
                                geometry.corner(3),
                                geometry.corner(2),
                                point))
                return true;
            return false;
        }
        else
            DUNE_THROW(Dune::NotImplemented,
                        "Intersection for point and geometry "
                        << type << " in three-dimensional world.");

    }

    //! Check whether a point is inside a given one-dimensional geometry
    template <class Geometry>
    static typename std::enable_if<Geometry::mydimension == 1, bool>::type
    pointInGeometry(const Geometry& geometry, const GlobalPosition& point)
    {
        return pointInInterval(geometry.corner(0),
                               geometry.corner(1),
                               point);
    }

    //! Find out whether a point is inside a tetrahedron (p0, p1, p2, p3)
    static bool pointInTetrahedron(const GlobalPosition& p0, const GlobalPosition& p1,
                                   const GlobalPosition& p2, const GlobalPosition& p3,
                                   const GlobalPosition& point)
    {
        // Algorithm from http://www.blackpawn.com/texts/pointinpoly/
        // See also "Real-Time Collision Detection" by Christer Ericson.
        // See also implementation of those algorithms by Anders Logg in FEniCS

        // put the tetrahedron points in an array
        const GlobalPosition *p[4] = {&p0, &p1, &p2, &p3};

        // iterate over all faces
        for (unsigned int i = 0; i < 4; ++i)
        {
            // compute all the vectors from vertex (local index 0) to the other points
            const GlobalPosition v1 = *p[(i + 1)%4] - *p[i];
            const GlobalPosition v2 = *p[(i + 2)%4] - *p[i];
            const GlobalPosition v3 = *p[(i + 3)%4] - *p[i];
            const GlobalPosition v = point - *p[i];
            // compute the normal to the facet (cross product)
            GlobalPosition n1 = crossProduct(v1, v2);
            n1 /= n1.two_norm();
            // find out on which side of the plane v and v3 are
            const double t1 = n1.dot(v);
            const double t2 = n1.dot(v3);
            // If the point is not exactly on the plane the
            // points have to be on the same side
            const double eps = eps_ * v1.two_norm();
            if ((t1 > eps || t1 < -eps) && std::signbit(t1) != std::signbit(t2))
                return false;
        }
        return true;
    }

    //! Find out whether a point is inside a triangle (p0, p1, p2)
    static bool pointInTriangle(const GlobalPosition& p0, const GlobalPosition& p1,
                                const GlobalPosition& p2, const GlobalPosition& point)
    {
        // Algorithm from http://www.blackpawn.com/texts/pointinpoly/
        // See also "Real-Time Collision Detection" by Christer Ericson.
        // See also implementation of those algorithms by Anders Logg in FEniCS

        // compute the vectors of the edges and the vector point-p0
        const GlobalPosition v1 = p0 - p2;
        const GlobalPosition v2 = p1 - p0;
        const GlobalPosition v3 = p2 - p1;
        const GlobalPosition v = point - p0;

        // compute the normal of the triangle
        const GlobalPosition n = crossProduct(v1, v2);

        // first check if we are in the plane of the triangle
        // if not we can return early
        const double t = v.dot(n);
        using std::abs;
        if (abs(t) > v1.two_norm()*eps_) // take |v1| as scale
            return false;

        // compute the normal to the triangle made of point and first edge
        // the dot product of this normal and the triangle normal has to
        // be positive because we defined the edges in the right orientation
        const GlobalPosition n1 = crossProduct(v, v1);
        const double t1 = n.dot(n1);
        if (t1 < 0) return false;

        const GlobalPosition n2 = crossProduct(v, v2);
        const double t2 = n.dot(n2);
        if (t2 < 0) return false;

        const GlobalPosition n3 = crossProduct(v, v3);
        const double t3 = n.dot(n3);
        if (t3 < 0) return false;

        return true;
    }

    //! Find out whether a point is inside an interval (p0, p1)
    static bool pointInInterval(const GlobalPosition& p0, const GlobalPosition& p1,
                                const GlobalPosition& point)
    {
        // compute the vectors between p0 and the other points
        const GlobalPosition v1 = p1 - p0;
        const GlobalPosition v2 = point - p0;

        // check if point and p0 are the same
        const double v1norm = v1.two_norm();
        const double v2norm = v2.two_norm();
        if (v2norm < v1norm*eps_)
            return true;

        // if not check if p0 and p1 are the same
        // then we know that point is not in the interval
        if (v1norm < eps_)
            return false;

        // if the cross product is zero the points are on a line
        const GlobalPosition n = crossProduct(v1, v2);

        // early return if the vector length is larger than zero
        if (n.two_norm() > v1norm*eps_)
            return false;

        // we know the points are aligned
        // if the dot product is positive and the length in range
        // the point is in the interval
        return (v1.dot(v2) > 0.0 && v2norm < v1norm*(1 + eps_));
    }

    /*!
     * \brief Check whether a point is in a bounding box
     * \param point The point
     * \param b Pointer to bounding box coordinates
     */
    static bool pointInBoundingBox(const Dune::FieldVector<double, 3>& point,
                                   const double* b)
    {
        const double eps0 = eps_*(b[3] - b[0]);
        const double eps1 = eps_*(b[4] - b[1]);
        const double eps2 = eps_*(b[5] - b[2]);
        return (b[0] - eps0 <= point[0] && point[0] <= b[3] + eps0 &&
                b[1] - eps1 <= point[1] && point[1] <= b[4] + eps1 &&
                b[2] - eps2 <= point[2] && point[2] <= b[5] + eps2);
    }

    /*!
     * \brief Check whether bounding box a collides with bounding box b
     * \param a, b Pointer to bounding box coordinates
     */
    static bool boundingBoxInBoundingBox(const double* a,
                                         const double* b)
    {
      const double eps0 = eps_*(b[3] - b[0]);
      const double eps1 = eps_*(b[4] - b[1]);
      const double eps2 = eps_*(b[5] - b[2]);
      return (b[0] - eps0 <= a[3] && a[0] <= b[3] + eps0 &&
              b[1] - eps1 <= a[4] && a[1] <= b[4] + eps1 &&
              b[2] - eps2 <= a[5] && a[2] <= b[5] + eps2);
    }

    //! Compute the bounding box of a vector of bounding boxes
    static void computeBBoxOfBBoxes(double* bBox,
                                    std::size_t& axis,
                                    const std::vector<double>& leafBoxes,
                                    const std::vector<unsigned int>::iterator& begin,
                                    const std::vector<unsigned int>::iterator& end)
    {
        // Copy the iterator and get coordinates of first box
        auto it = begin;
        const double* bFirst = leafBoxes.data() + 6*(*it);
        // Maybe write out the loop for optimization
        for (std::size_t coordIdx = 0; coordIdx < 6; ++coordIdx)
            bBox[coordIdx] = bFirst[coordIdx];

        // Compute min and max with the remaining boxes
        for (; it != end; ++it)
        {
            const double* b = leafBoxes.data() + 6*(*it);
            if (b[0] < bBox[0]) bBox[0] = b[0];
            if (b[1] < bBox[1]) bBox[1] = b[1];
            if (b[2] < bBox[2]) bBox[2] = b[2];
            if (b[3] > bBox[3]) bBox[3] = b[3];
            if (b[4] > bBox[4]) bBox[4] = b[4];
            if (b[5] > bBox[5]) bBox[5] = b[5];
        }

        // Compute the longest axis
        const double x = bBox[3] - bBox[0];
        const double y = bBox[4] - bBox[1];
        const double z = bBox[5] - bBox[2];

        if (x > y && x > z)
            axis = 0;
        else if (y > z)
            axis = 1;
        else
            axis = 2;
    }

    //! Sort the bounding boxes along the longest axis
    static void sortBoundingBoxes(std::size_t axis,
                                  const std::vector<double>& leafBoxes,
                                  const std::vector<unsigned int>::iterator& begin,
                                  const std::vector<unsigned int>::iterator& middle,
                                  const std::vector<unsigned int>::iterator& end)
    {
        switch (axis)
        {
            case 0:
                std::nth_element(begin, middle, end, lessXBox(leafBoxes));
            case 1:
                std::nth_element(begin, middle, end, lessYBox(leafBoxes));
            default:
                std::nth_element(begin, middle, end, lessZBox(leafBoxes));
        }
    }

    /*!
     * \brief Comparison function for sorting bounding boxes on the x-axis
     * \note This could be replaced by lambdas
     */
    struct lessXBox
    {
        const std::vector<double>& bBoxes;
        lessXBox(const std::vector<double>& bBoxes_): bBoxes(bBoxes_) {}
        inline bool operator()(unsigned int i, unsigned int j)
        {
            const double* bi = bBoxes.data() + 6*i;
            const double* bj = bBoxes.data() + 6*j;
            return bi[0] + bi[3] < bj[0] + bj[3];
        }
    };

    /*!
     * \brief Comparison function for sorting bounding boxes on the y-axis
     * \note This could be replaced by lambdas
     */
    struct lessYBox
    {
        const std::vector<double>& bBoxes;
        lessYBox(const std::vector<double>& bBoxes_): bBoxes(bBoxes_) {}
        inline bool operator()(unsigned int i, unsigned int j)
        {
            const double* bi = bBoxes.data() + 6*i;
            const double* bj = bBoxes.data() + 6*j;
            return bi[1] + bi[4] < bj[1] + bj[4];
        }
    };

    /*!
     * \brief Comparison function for sorting bounding boxes on the z-axis
     * \note This could be replaced by lambdas
     */
    struct lessZBox
    {
        const std::vector<double>& bBoxes;
        lessZBox(const std::vector<double>& bBoxes_): bBoxes(bBoxes_) {}
        inline bool operator()(unsigned int i, unsigned int j)
        {
            const double* bi = bBoxes.data() + 6*i;
            const double* bj = bBoxes.data() + 6*j;
            return bi[2] + bi[5] < bj[2] + bj[5];
        }
    };
};

//! Helper methods for world dimension two
template <>
class BoundingBoxTreeHelper<2>
{
    // An epsilon for floating point operations
    static constexpr double eps_ = 1.0e-7;
    typedef Dune::FieldVector<double, 2> GlobalPosition;

public:
    // Check whether a point is inside a given geometry
    template <class Geometry>
    static typename std::enable_if<Geometry::mydimension == 2, bool>::type
    pointInGeometry(const Geometry& geometry, const GlobalPosition& point)
    {
        // get the geometry type
        auto type = geometry.type();

        // if it's a triangle we can check directly
        if (type.isTriangle())
        {
            return pointInTriangle(geometry.corner(0),
                                   geometry.corner(1),
                                   geometry.corner(2),
                                   point);
        }
        // split quadrilaterals into two triangles
        else if (type.isQuadrilateral())
        {
            if (pointInTriangle(geometry.corner(0),
                                geometry.corner(1),
                                geometry.corner(3),
                                point))
                return true;
            if (pointInTriangle(geometry.corner(0),
                                geometry.corner(3),
                                geometry.corner(2),
                                point))
                return true;
            return false;
        }
        else
            DUNE_THROW(Dune::NotImplemented,
                        "Intersection for point and geometry "
                        << type << " in two-dimensional world.");
    }

    template <class Geometry>
    static typename std::enable_if<Geometry::mydimension == 1, bool>::type
    pointInGeometry(const Geometry& geometry, const GlobalPosition& point)
    {
        return pointInInterval(geometry.corner(0),
                               geometry.corner(1),
                               point);
    }

    //! find out if a point is inside a triangle
    static bool pointInTriangle(const GlobalPosition& p0, const GlobalPosition& p1,
                                const GlobalPosition& p2, const GlobalPosition& point)
    {
        // Use barycentric coordinates
        const double A = 0.5*(-p1[1]*p2[0] + p0[1]*(p2[0] - p1[0])
                                           + p0[0]*(p1[1] - p2[1]) + p1[0]*p2[1]);
        const double sign = std::copysign(1.0, A);
        const double s = sign*(p0[1]*p2[0] + point[0]*(p2[1]-p0[1])
                              -p0[0]*p2[1] + point[1]*(p0[0]-p2[0]));
        const double t = sign*(p0[0]*p1[1] + point[0]*(p0[1]-p1[1])
                              -p0[1]*p1[0] + point[1]*(p1[0]-p0[0]));
        const double eps = (p0 - p1).two_norm()*eps_;

        return (s > -eps
                && t > -eps
                && (s + t) < 2*A*sign + eps);
    }

    //! find out if a point is inside an interval
    static bool pointInInterval(const GlobalPosition& p0, const GlobalPosition& p1,
                                const GlobalPosition& point)
    {
        // compute the vectors between p0 and the other points
        const GlobalPosition v1 = p1 - p0;
        const GlobalPosition v2 = point - p0;

        // check if point and p0 are the same
        const double v1norm = v1.two_norm();
        const double v2norm = v2.two_norm();
        if (v2norm < v1norm*eps_)
            return true;

        // if not check if p0 and p1 are the same
        // then we know that point is not in the interval
        if (v1norm < eps_)
            return false;

        // if the cross product is zero the points are on a line
        const double n = crossProduct(v1, v2);

        // early return if the cross product is larger than zero
        using std::abs;
        if (abs(n) > v1norm*eps_)
            return false;

        // we know the points are aligned
        // if the dot product is positive and the length in range
        // the point is in the interval
        return (v1.dot(v2) > 0.0 && v2norm < v1norm*(1 + eps_));
    }

    /*!
     * \brief Check whether a point is in a bounding box
     * \param point The point
     * \param b Pointer to bounding box coordinates
     */
    static bool pointInBoundingBox(const Dune::FieldVector<double, 2>& point,
                                   const double* b)
    {
        const double eps0 = eps_*(b[2] - b[0]);
        const double eps1 = eps_*(b[3] - b[1]);
        return (b[0] - eps0 <= point[0] && point[0] <= b[2] + eps0 &&
                b[1] - eps1 <= point[1] && point[1] <= b[3] + eps1);
    }

    /*!
     * \brief Check whether bounding box a collides with bounding box b
     * \param a, b Pointer to bounding box coordinates
     */
    static bool boundingBoxInBoundingBox(const double* a,
                                         const double* b)
    {
      const double eps0 = eps_*(b[2] - b[0]);
      const double eps1 = eps_*(b[3] - b[1]);
      return (b[0] - eps0 <= a[2] && a[0] <= b[2] + eps0 &&
              b[1] - eps1 <= a[3] && a[1] <= b[3] + eps1);
    }

    //! Compute the bounding box of a vector of bounding boxes
    static void computeBBoxOfBBoxes(double* bBox,
                                    std::size_t& axis,
                                    const std::vector<double>& leafBoxes,
                                    const std::vector<unsigned int>::iterator& begin,
                                    const std::vector<unsigned int>::iterator& end)
    {
        // Copy the iterator and get coordinates of first box
        auto it = begin;
        const double* bFirst = leafBoxes.data() + 4*(*it);
        // Maybe write out the loop for optimization
        for (std::size_t coordIdx = 0; coordIdx < 4; ++coordIdx)
            bBox[coordIdx] = bFirst[coordIdx];

        // Compute min and max with the remaining boxes
        for (; it != end; ++it)
        {
            const double* b = leafBoxes.data() + 4*(*it);
            if (b[0] < bBox[0]) bBox[0] = b[0];
            if (b[1] < bBox[1]) bBox[1] = b[1];
            if (b[2] > bBox[2]) bBox[2] = b[2];
            if (b[3] > bBox[3]) bBox[3] = b[3];
        }

        // Compute the longest axis
        const double x = bBox[2] - bBox[0];
        const double y = bBox[3] - bBox[1];

        if (x > y)
            axis = 0;
        else
            axis = 1;
    }

    //! Sort the bounding boxes along the longest axis
    static void sortBoundingBoxes(std::size_t axis,
                                  const std::vector<double>& leafBoxes,
                                  const std::vector<unsigned int>::iterator& begin,
                                  const std::vector<unsigned int>::iterator& middle,
                                  const std::vector<unsigned int>::iterator& end)
    {
        if (axis == 0)
            std::nth_element(begin, middle, end, lessXBox(leafBoxes));
        else
            std::nth_element(begin, middle, end, lessYBox(leafBoxes));
    }

    /*!
     * \brief Comparison function for sorting bounding boxes on the x-axis
     * \note This could be replaced by lambdas
     */
    struct lessXBox
    {
        const std::vector<double>& bBoxes;
        lessXBox(const std::vector<double>& bBoxes_): bBoxes(bBoxes_) {}
        inline bool operator()(unsigned int i, unsigned int j)
        {
            const double* bi = bBoxes.data() + 4*i;
            const double* bj = bBoxes.data() + 4*j;
            return bi[0] + bi[2] < bj[0] + bj[2];
        }
    };

    /*!
     * \brief Comparison function for sorting bounding boxes on the y-axis
     * \note This could be replaced by lambdas
     */
    struct lessYBox
    {
        const std::vector<double>& bBoxes;
        lessYBox(const std::vector<double>& bBoxes_): bBoxes(bBoxes_) {}
        inline bool operator()(unsigned int i, unsigned int j)
        {
            const double* bi = bBoxes.data() + 4*i;
            const double* bj = bBoxes.data() + 4*j;
            return bi[1] + bi[3] < bj[1] + bj[3];
        }
    };
};

//! Helper methods for world dimension one
template <>
class BoundingBoxTreeHelper<1>
{
    // An epsilon for floating point operations
    static constexpr double eps_ = 1.0e-7;
    typedef Dune::FieldVector<double, 1> GlobalPosition;

public:
    // Check whether a point is inside a given geometry
    template <class Geometry>
    static bool pointInGeometry(const Geometry& geometry,
                                const GlobalPosition& point)
    {
        return pointInInterval(geometry.corner(0),
                               geometry.corner(1),
                               point);
    }

    //! find out if a point is inside an interval
    static bool pointInInterval(const GlobalPosition& p0, const GlobalPosition& p1,
                                const GlobalPosition& point)
    {
        const double v1 = point[0] - p0[0];
        const double v2 = p1[0] - p0[0];
        // the coordinates are the same
        if (v1 < v2*eps_)
            return true;
        // the point doesn't coincide with p0
        // so if p0 and p1 are equal it's not inside
        if (v2 < eps_)
            return false;

        // the point is inside if the length is
        // small than the interval length and the
        // sign of v1 & v2 are the same
        using std::abs;
        using std::signbit;
        return (signbit(v1) == signbit(v2)
                && abs(v1) < abs(v2)*(1 + eps_));
    }

    /*!
     * \brief Check whether a point is in a bounding box
     * \param point The point
     * \param b Pointer to bounding box coordinates
     */
    static bool pointInBoundingBox(const GlobalPosition& point,
                                   const double* b)
    {
        const double eps0 = eps_*(b[1] - b[0]);
        return b[0] - eps0 <= point[0] && point[0] <= b[1] + eps0;
    }

    /*!
     * \brief Check whether bounding box a collides with bounding box b
     * \param a, b Pointer to bounding box coordinates
     */
    static bool boundingBoxInBoundingBox(const double* a,
                                         const double* b)
    {
      const double eps0 = eps_*(b[1] - b[0]);
      return b[0] - eps0 <= a[1] && a[0] <= b[1] + eps0;
    }

    //! Compute the bounding box of a vector of bounding boxes
    static void computeBBoxOfBBoxes(double* bBox,
                                    std::size_t& axis,
                                    const std::vector<double>& leafBoxes,
                                    const std::vector<unsigned int>::iterator& begin,
                                    const std::vector<unsigned int>::iterator& end)
    {
        // Copy the iterator and get coordinates of first box
        auto it = begin;
        const double* bFirst = leafBoxes.data() + 2*(*it);
        bBox[0] = bFirst[0];
        bBox[1] = bFirst[1];

        // Compute min and max with the remaining boxes
        for (; it != end; ++it)
        {
            const double* b = leafBoxes.data() + 2*(*it);
            if (b[0] < bBox[0]) bBox[0] = b[0];
            if (b[1] > bBox[1]) bBox[1] = b[1];
        }

        // Compute the longest axis
        axis = 0;
    }

    //! Sort the bounding boxes along the longest axis
    static void sortBoundingBoxes(std::size_t axis,
                                  const std::vector<double>& leafBoxes,
                                  const std::vector<unsigned int>::iterator& begin,
                                  const std::vector<unsigned int>::iterator& middle,
                                  const std::vector<unsigned int>::iterator& end)
    { std::nth_element(begin, middle, end, lessXBox(leafBoxes)); }

    /*!
     * \brief Comparison function for sorting bounding boxes on the x-axis
     * \note This could be replaced by lambdas
     */
    struct lessXBox
    {
        const std::vector<double>& bBoxes;
        lessXBox(const std::vector<double>& bBoxes_): bBoxes(bBoxes_) {}
        inline bool operator()(unsigned int i, unsigned int j)
        {
            const double* bi = bBoxes.data() + 2*i;
            const double* bj = bBoxes.data() + 2*j;
            return bi[0] + bi[1] < bj[0] + bj[1];
        }
    };
};

/*!
 * \brief An intersection object resulting from the collision of two bounding box tree primitives
 *
 * After is has been found that two leaf bounding boxes intersect a primitive test has to be
 * performed to see if the actual entities inside the bounding box intersect too. The result
 * if such an intersection is found as an object of this class containing the indices of the
 * intersecting entities and the corners of the intersection object.
 */
template<class GridView1, class GridView2>
class BoundingBoxTreeIntersection
{
    const static int dimworld = GridView1::dimensionworld;
    using Scalar = typename GridView1::ctype;
    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

public:
    BoundingBoxTreeIntersection(unsigned int a, unsigned int b, std::vector<GlobalPosition>&& c)
    : a_(a), b_(b), corners_(c) {}

    //! Get the index of the intersecting entity belonging to this grid
    inline unsigned int first() const
    { return a_; }

    //! Get the index of the intersecting entity belonging to the other grid
    inline unsigned int second() const
    { return b_; }

    //! Get the corners of the intersection geometry
    inline std::vector<GlobalPosition> corners() const
    { return corners_; }

    /*!
     * \brief Check if the corners of this intersection match with the given corners
     * \note This is useful to check if the intersection geometry of two intersections coincide.
     */
    bool cornersMatch(const std::vector<GlobalPosition>& corners) const
    {
        if (corners.size() != corners_.size())
            return false;

        const auto eps = 1.5e-7*(corners_[1] - corners_[0]).two_norm();
        for (int i = 0; i < corners_.size(); ++i)
            if ((corners[0] - corners[1]).two_norm() > eps)
                return false;

        return true;
    }

private:
    unsigned int a_, b_; //!< Indices of the intersection elements
    std::vector<GlobalPosition> corners_; //!< the corner points of the intersection geometry
};

/*!
 * \brief A class mapping from an element index to elements using element seeds
 */
template <class GridView>
class IndexToElementMap
  : public std::vector<typename GridView::Traits::Grid::template Codim<0>::EntitySeed>
{
    typedef typename GridView::Traits::Grid Grid;
    typedef typename GridView::template Codim<0>::Entity Element;
public:
    IndexToElementMap(const GridView& gridView)
      : grid_(gridView.grid()) {}

    //! get an element from an index t
    template<class T>
    Element entity(T&& t)
    { return grid_.entity((*this)[std::forward<T>(t)]); }

private:
    const Grid& grid_;
};

/*!
 * \brief An axis-aligned bounding box volume tree implementation
 *
 * The class constructs a hierarchical structure of bounding box volumes around
 * grid entities. This class can be used to efficiently compute intersections
 * between a grid and other geometrical object. It only implements the intersection
 * of two of such bounding box trees, so that two independent grids can be intersected.
 */
template <class GridView>
class BoundingBoxTree
{
    //! be friends with all other kinds bounding box trees so that
    //! they can call each others private methods
    template <class OtherGridView> friend class BoundingBoxTree;

    static const int dim = GridView::dimension;
    static const int dimworld = GridView::dimensionworld;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = typename GridView::ctype;
    using GlobalPosition = typename Dune::FieldVector<Scalar, dimworld>;

public:
    //! Default Constructor
    BoundingBoxTree() {}

    //! Constructor with gridView
    BoundingBoxTree(const GridView& leafGridView)
    { build(leafGridView); }

    //! Build up bounding box tree for a grid with leafGridView
    void build(const GridView& leafGridView)
    {
        // clear data if any
        clear_();

        // Create bounding boxes for all elements
        const unsigned int numLeaves = leafGridView.size(0);

        // start the timer
        Dune::Timer timer;

        // create a vector for leaf boxes (min and max for all dims)
        std::vector<double> leafBoxes(2*dimworld*numLeaves);

        // Create index to element map
        indexToElementMap_ = std::make_shared<IndexToElementMap<GridView> >(leafGridView);
        indexToElementMap_->resize(numLeaves);

        for (auto&& element : elements(leafGridView))
        {
            unsigned int eIdx = leafGridView.indexSet().index(element);
            computeEntityBoundingBox_(leafBoxes.data() + 2*dimworld*eIdx, element);
            (*indexToElementMap_)[eIdx] = element.seed();
        }

        // Create the leaf partition, the set of available indices (to be sorted)
        std::vector<unsigned int> leafPartition(numLeaves);
        for (unsigned int i = 0; i < numLeaves; ++i)
            leafPartition[i] = i;

        // Recursively build the bounding box tree
        build_(leafBoxes, leafPartition.begin(), leafPartition.end());

        // We are done, log output
        std::cout << "Computed bounding box tree with " << numBoundingBoxes_()
                  << " nodes for " << numLeaves << " grid entites in "
                  << timer.stop() << " seconds." << std::endl;
    }

    //! Compute all intersections between entities and a point
    std::vector<unsigned int> computeEntityCollisions(const Dune::FieldVector<double, dimworld>& point) const
    {
        // Call the recursive find function to find candidates
        std::vector<unsigned int> entities;
        computeCollisions_(point, numBoundingBoxes_() - 1, entities);
        return entities;
    }

    //! Compute all intersections between entities and another bounding box tree
    template<class OtherGridView>
    std::vector<BoundingBoxTreeIntersection<GridView, OtherGridView>>
    computeEntityCollisions(const BoundingBoxTree<OtherGridView>& otherTree) const
    {
        // check if the world dimensions match
        static_assert(dimworld == OtherGridView::dimensionworld, "Can only collide bounding box trees of same world dimension");

        // Create data structure for return type
        std::vector<BoundingBoxTreeIntersection<GridView, OtherGridView>> intersections;

        // Call the recursive find function to find candidates
        computeCollisions_(otherTree,
                           this->numBoundingBoxes_() - 1,
                           otherTree.numBoundingBoxes_() -1,
                           intersections);

        return intersections;
    }

    //! Get an element from a global element index
    Element entity(unsigned int eIdx) const
    { return indexToElementMap_->entity(eIdx); }

private:
    /*!
     * \brief Bounding box data structure
     * Leaf nodes are indicated by setting child_0 to
     * the node itself and child_1 is the index of the entity in the bounding box.
     */
    struct BoundingBox
    {
        unsigned int child_0;
        unsigned int child_1;
    };

    //! Vector of bounding boxes
    std::vector<BoundingBox> boundingBoxes_;

    //! Vector of bounding box coordinates
    std::vector<double> boundingBoxCoordinates_;

    //! Shared pointer to the index to element map
    std::shared_ptr<IndexToElementMap<GridView> > indexToElementMap_;

    //! Clear all data
    void clear_()
    {
        boundingBoxes_.clear();
        boundingBoxCoordinates_.clear();
        if(indexToElementMap_) indexToElementMap_->clear();
    }

    //! Build bounding box tree for all entities recursively
    unsigned int build_(const std::vector<double>& leafBoxes,
                        const std::vector<unsigned int>::iterator& begin,
                        const std::vector<unsigned int>::iterator& end)
    {
        assert(begin < end);

        // Create empty bounding box
        BoundingBox bBox;

        // If we reached the leaf
        if (end - begin == 1)
        {
            // Get the bounding box coordinates for the leaf
            const unsigned int eIdx = *begin;
            const double* b = leafBoxes.data() + 2*dimworld*eIdx;

            // Store the data in the bounding box
            bBox.child_0 = numBoundingBoxes_();
            bBox.child_1 = eIdx;
            return addBoundingBox_(bBox, b);
        }

        // Compute the bounding box of all bounding boxes in the range [begin, end]
        double b[dimworld*2];
        std::size_t axis;
        BoundingBoxTreeHelper<dimworld>::computeBBoxOfBBoxes(b, axis, leafBoxes, begin, end);

        // Sort bounding boxes along the longest axis
        std::vector<unsigned int>::iterator middle = begin + (end - begin)/2;
        BoundingBoxTreeHelper<dimworld>::sortBoundingBoxes(axis, leafBoxes, begin, middle, end);

        // Split the bounding boxes into two branches and call build recursively
        bBox.child_0 = build_(leafBoxes, begin, middle);
        bBox.child_1 = build_(leafBoxes, middle, end);

        // Store the bounding box data. The root will be added at the end.
        return addBoundingBox_(bBox, b);
    }

    //! Compute collisions with point recursively
    void computeCollisions_(const Dune::FieldVector<double, dimworld>& point,
                            unsigned int node,
                            std::vector<unsigned int>& entities) const
    {
        // Get the bounding box for the current node
        const BoundingBox& bBox = getBoundingBox_(node);

        // if the point is not in the bounding box we can stop
        if (!BoundingBoxTreeHelper<dimworld>::pointInBoundingBox(point, getBoundingBoxCoordinates_(node)))
            return;

        // We know now it's inside. If the box is a leaf add it.
        else if (isLeaf_(bBox, node))
        {
            // but add it only if the point is also inside the entity
            const unsigned int eIdx = bBox.child_1;
            auto geometry = (indexToElementMap_->entity(eIdx)).geometry();
            if (BoundingBoxTreeHelper<dimworld>::pointInGeometry(geometry, point))
                entities.push_back(eIdx);

            // const ReferenceElement &refElement = ReferenceElements::general(geometry.type());
            // if (refElement.checkInside(geometry.local(point)))
            //     entities.push_back(eIdx);
        }

        // No leaf. Check both children.
        else
        {
            computeCollisions_(point, bBox.child_0, entities);
            computeCollisions_(point, bBox.child_1, entities);
        }
    }

    //! Compute collisions with other bounding box tree recursively
    template <class OtherGridView>
    void computeCollisions_(const BoundingBoxTree<OtherGridView>& treeB,
                            unsigned int nodeA,
                            unsigned int nodeB,
                            std::vector<BoundingBoxTreeIntersection<GridView, OtherGridView>>& intersections) const
    {
        // get alias
        const auto& treeA = *this;

        // Get the bounding box for the current node
        const auto& bBoxA = treeA.getBoundingBox_(nodeA);
        const auto& bBoxB = treeB.getBoundingBox_(nodeB);

        // if the two bounding boxes don't collide we can stop searching
        if (!BoundingBoxTreeHelper<dimworld>::
             boundingBoxInBoundingBox(treeA.getBoundingBoxCoordinates_(nodeA),
                                      treeB.getBoundingBoxCoordinates_(nodeB)))
            return;

        // Check if we have a leaf in treeA or treeB
        const bool isLeafA = treeA.isLeaf_(bBoxA, nodeA);
        const bool isLeafB = treeB.isLeaf_(bBoxB, nodeB);

        // If both boxes are leaves and collide add them
        if (isLeafA && isLeafB)
        {
            const unsigned int eIdxA = bBoxA.child_1;
            const unsigned int eIdxB = bBoxB.child_1;

            auto geometryA = treeA.entity(eIdxA).geometry();
            auto geometryB = treeB.entity(eIdxB).geometry();

            using CollisionType = GeometryCollision<decltype(geometryA), decltype(geometryB)>;
            std::vector<GlobalPosition> intersection;
            if (CollisionType::collide(geometryA, geometryB, intersection))
                intersections.emplace_back(eIdxA, eIdxB, std::move(intersection));
        }

        // if we reached the leaf in treeA, just continue in treeB
        else if (isLeafA)
        {
            computeCollisions_(treeB, nodeA, bBoxB.child_0, intersections);
            computeCollisions_(treeB, nodeA, bBoxB.child_1, intersections);
        }

        // if we reached the leaf in treeB, just continue in treeA
        else if (isLeafB)
        {
            computeCollisions_(treeB, bBoxA.child_0, nodeB, intersections);
            computeCollisions_(treeB, bBoxA.child_1, nodeB, intersections);
        }

        // we know now that both trees didn't reach the leaf yet so
        // we continue with the larger tree first (bigger node number)
        else if (nodeA > nodeB)
        {
            computeCollisions_(treeB, bBoxA.child_0, nodeB, intersections);
            computeCollisions_(treeB, bBoxA.child_1, nodeB, intersections);
        }
        else
        {
            computeCollisions_(treeB, nodeA, bBoxB.child_0, intersections);
            computeCollisions_(treeB, nodeA, bBoxB.child_1, intersections);
        }
    }

    //! Add a new bounding box to the tree
    inline unsigned int addBoundingBox_(const BoundingBox& bBox,
                                        const double* b)
    {
        // Add the bounding box
        boundingBoxes_.push_back(bBox);

        // Add the bounding box's coordinates
        for (std::size_t i = 0; i < 2*dimworld; ++i)
            boundingBoxCoordinates_.push_back(b[i]);

        // return the index of the new node
        return boundingBoxes_.size() - 1;
    }

    //! Get an existing bounding box for a given node
    inline const BoundingBox& getBoundingBox_(unsigned int node) const
    { return boundingBoxes_[node]; }

    //! Get an existing bounding box for a given node
    const double* getBoundingBoxCoordinates_(unsigned int node) const
    { return boundingBoxCoordinates_.data() + 2*dimworld*node; }

    //! Get the number of bounding boxes currently in the tree
    inline std::size_t numBoundingBoxes_() const
    { return boundingBoxes_.size(); }

    //! Check whether a bounding box is a leaf node
    //! Leaf nodes have itself as child_0
    inline bool isLeaf_(const BoundingBox& bBox, unsigned int node) const
    { return bBox.child_0 == node; }

    //! Compute the bounding box of a grid entity
    template <class Entity>
    void computeEntityBoundingBox_(double* b, const Entity& entity) const
    {
        // get the bounding box coordinates
        double* xMin = b;
        double* xMax = b + dimworld;

        // get mesh entity data
        auto geometry = entity.geometry();

        // Get coordinates of first vertex
        auto corner = geometry.corner(0);
        for (std::size_t dimIdx = 0; dimIdx < dimworld; ++dimIdx)
            xMin[dimIdx] = xMax[dimIdx] = corner[dimIdx];

        // Compute the min and max over the remaining vertices
        for (std::size_t vLocalIdx = 1; vLocalIdx < entity.subEntities(dim); ++vLocalIdx)
        {
            corner = geometry.corner(vLocalIdx);
            for (std::size_t dimIdx = 0; dimIdx < dimworld; ++dimIdx)
            {
                using std::max;
                using std::min;
                xMin[dimIdx] = min(xMin[dimIdx], corner[dimIdx]);
                xMax[dimIdx] = max(xMax[dimIdx], corner[dimIdx]);
            }
        }
    }
};

} // end namespace Dumux

#endif
