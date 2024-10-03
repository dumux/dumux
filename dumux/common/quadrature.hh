// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Quadrature
 * \brief Quadrature rules over different geometries.
 */

#ifndef DUMUX_COMMON_QUADRATURE_HH
#define DUMUX_COMMON_QUADRATURE_HH

#include <dune/geometry/quadraturerules.hh>

#include <dumux/geometry/intersectingentities.hh>
#include <dumux/geometry/volume.hh>
#include <dumux/geometry/center.hh>

namespace Dumux::Quadrature {

/*!
 * \ingroup Quadrature
 * \brief Quadrature point assigned to a given grid entity with given index
 */
template<typename ct, int dim, typename IndexType>
class GridQuadraturePoint : public Dune::QuadraturePoint<ct, dim>
{
    using ParentType = Dune::QuadraturePoint<ct, dim>;

public:
    GridQuadraturePoint (const typename ParentType::Vector& x, ct w, IndexType idx)
    : ParentType(x,w), index_(idx)
    {}

    const IndexType index() const
    { return index_;}

private:
    IndexType index_;
};

/*!
 * \ingroup Quadrature
 * \brief Helper class to integrate over intersecting entities
 * \note The algorithm creates quadrature points based on the sub-triangulated intersections
 */
template<class Scalar, int dim, typename IndexType>
class IntersectingEntitiesRule : public std::vector< GridQuadraturePoint<Scalar, dim, IndexType> >
{
    using QuadraturePoint = GridQuadraturePoint<Scalar, dim, IndexType>;

public:
    /*!
     * \brief Constructor
     * \note IntersectionInfo is assumed to be based on sub-triangulated intersections
     */
    template<class CoordTypeA, class CoordTypeB>
    IntersectingEntitiesRule(const std::vector<IntersectionInfo<dim, CoordTypeA, CoordTypeB>>& intersectingEntities)
    {
        numPoints_ = intersectingEntities.size();
        initialize_(intersectingEntities);
    }

    int size() const
    { return numPoints_;}

private:
    template<class IntersectionConterainer>
    void initialize_(const IntersectionConterainer& intersectingEntities)
    {
        for(const auto& intersection : intersectingEntities)
        {
            const auto& corners = intersection.corners();
            if(corners.size() != dim + 1)
                DUNE_THROW(Dune::NotImplemented, "Quadrature rule on intersecting entities assumes sub-triangulated intersections.");

            this->push_back(QuadraturePoint(Dumux::center(corners),
                                            Dumux::convexPolytopeVolume<dim>(Dune::GeometryTypes::simplex(dim),
                                                                             [corners](unsigned int i){ return corners[i];}),
                                            intersection.second()));
        }

    }

    int numPoints_;
};

} // end namespace Dumux::Quadrature

#endif
