// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FEMDiscretization
 * \brief The grid geometry class for models using finite element schemes.
 *        This is basically a wrapper around a function space basis.
 */
#ifndef DUMUX_DISCRETIZATION_FE_GRID_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_FE_GRID_GEOMETRY_HH

#include <unordered_map>

#include <dumux/common/indextraits.hh>
#include <dumux/common/defaultmappertraits.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/checkoverlapsize.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/fem/feelementgeometry.hh>

namespace Dumux {

/*!
 * \ingroup FEMDiscretization
 * \brief Default Traits class for the fem grid geometry.
 * \tparam The finite element function space basis
 * \tparam MapperTraits Traits class containing data types for mappers
 */
template<class FEBasis, class MapperTraits = DefaultMapperTraits<typename FEBasis::GridView>>
struct DefaultFEGridGeometryTraits : public MapperTraits
{
    template<class GridGeometry>
    using LocalView = FEElementGeometry<GridGeometry>;
};

/*!
 * \ingroup FEMDiscretization
 * \brief The grid geometry class for models using finite element schemes.
 *        This is basically a wrapper around a function space basis.
 * \tparam FEB The finite element function space basis
 * \tparam MapperTraits Traits class containing data types for mappers
 */
template<class FEB, class Traits = DefaultFEGridGeometryTraits<FEB>>
class FEGridGeometry
: public BaseGridGeometry< typename FEB::GridView, Traits >
{
    using ThisType = FEGridGeometry<FEB, Traits>;
    using ParentType = BaseGridGeometry<typename FEB::GridView, Traits>;

    using GridIndexType = typename IndexTraits<typename FEB::GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<typename FEB::GridView>::LocalIndex;

public:
    //! export the discretization method this geometry belongs to
    using DiscretizationMethod = DiscretizationMethods::FEM;
    static constexpr DiscretizationMethod discMethod{};

    //! export the grid view type
    using GridView = typename FEB::GridView;
    //! export the type of extrusion
    using Extrusion = Extrusion_t<Traits>;
    //! export the type of finite element basis
    using FEBasis = FEB;
    //! export local view
    using LocalView = typename Traits::template LocalView<ThisType>;

    //! Constructor
    FEGridGeometry(std::shared_ptr<FEBasis> feBasis)
    : ParentType(feBasis->gridView())
    , feBasis_(feBasis)
    {
        // Check if the overlap size is what we expect
        if (!CheckOverlapSize<DiscretizationMethod>::isValid(*feBasis))
            DUNE_THROW(Dune::InvalidStateException, "The finite element discretization method only works with zero overlap for parallel computations. "
                                                     << " Set the parameter \"Grid.Overlap\" in the input file.");
    }

    //! The total number of degrees of freedom
    auto numDofs() const
    { return feBasis_->size(); }

    //! The total number of degrees of freedom
    const FEBasis& feBasis() const
    { return *feBasis_; }

    //! If a vertex / d.o.f. is on a periodic boundary
    bool dofOnPeriodicBoundary(GridIndexType dofIdx) const
    { DUNE_THROW(Dune::NotImplemented, "Periodic BC support for FEM schemes"); }

    //! The index of the vertex / d.o.f. on the other side of the periodic boundary
    GridIndexType periodicallyMappedDof(GridIndexType dofIdx) const
    { DUNE_THROW(Dune::NotImplemented, "Periodic BC support for FEM schemes"); }

    //! Returns the map between dofs across periodic boundaries
    const std::unordered_map<GridIndexType, GridIndexType>& periodicVertexMap() const
    { DUNE_THROW(Dune::NotImplemented, "Periodic BC support for FEM schemes"); }

private:
    std::shared_ptr<FEBasis> feBasis_;
};

} // end namespace Dumux

#endif
