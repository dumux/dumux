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
 * \ingroup FEMDiscretization
 * \brief The grid geometry class for models using finite element schemes.
 *        This is basically a wrapper around a function space basis.
 */
#ifndef DUMUX_DISCRETIZATION_FE_GRID_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_FE_GRID_GEOMETRY_HH

#include <type_traits>

#include <dumux/common/indextraits.hh>
#include <dumux/common/defaultmappertraits.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/checkoverlapsize.hh>
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
 * \tparam ASB The finite element function space basis of the ansatz space
 * \tparam TSB The finite element function space basis of the trial space
 * \tparam MapperTraits Traits class containing data types for mappers
 */
template<class ASB, class TSB = ASB, class Traits = DefaultFEGridGeometryTraits<ASB>>
class FEGridGeometry
: public BaseGridGeometry< typename ASB::GridView, Traits >
{
    using ThisType = FEGridGeometry<ASB, TSB, Traits>;
    using ParentType = BaseGridGeometry<typename ASB::GridView, Traits>;

    using GridIndexType = typename IndexTraits<typename ASB::GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<typename ASB::GridView>::LocalIndex;

    static constexpr bool standardGalerkin = std::is_same<ASB, TSB>::value;
    static_assert(std::is_same<typename ASB::GridView, typename TSB::GridView>::value,
                  "Trial and Ansatz space must operate on the same grid!");

public:
    //! export discretization method
    static constexpr DiscretizationMethod discMethod = DiscretizationMethod::fem;

    //! export the grid view type
    using GridView = typename ASB::GridView;
    //! export the type of finite element basis of the ansatz space
    using FEBasis = ASB;
    //! export the type of finite element basis of the ansatz space
    using AnsatzSpaceBasis = FEBasis;
    //! export the type of finite element basis of the trial space
    using TrialSpaceBasis = TSB;
    //! export local view
    using LocalView = typename Traits::template LocalView<ThisType>;

    //! Constructor for standard Galerkin type settings
    template<bool isSG = standardGalerkin, std::enable_if_t<isSG, int> = 0>
    FEGridGeometry(std::shared_ptr<AnsatzSpaceBasis> ansatzSpaceBasis)
    : ParentType(ansatzSpaceBasis->gridView())
    , ansatzSpaceBasis_(ansatzSpaceBasis)
    , trialSpaceBasis_(ansatzSpaceBasis)
    {
        // Check if the overlap size is what we expect
        if (!CheckOverlapSize<DiscretizationMethod::fem>::isValid(ansatzSpaceBasis->gridView()))
            DUNE_THROW(Dune::InvalidStateException, "The finite element discretization method only works with zero overlap for parallel computations. "
                                                     << " Set the parameter \"Grid.Overlap\" in the input file.");
    }

    //! Constructor for standard Galerkin type settings
    template<bool isSG = standardGalerkin, std::enable_if_t<isSG, int> = 0>
    FEGridGeometry(std::shared_ptr<AnsatzSpaceBasis> ansatzSpaceBasis,
                   std::shared_ptr<TrialSpaceBasis> trialSpaceBasis)
    : ParentType(ansatzSpaceBasis->gridView())
    , ansatzSpaceBasis_(ansatzSpaceBasis)
    , trialSpaceBasis_(trialSpaceBasis)
    {
        if (ansatzSpaceBasis->gridView().size(0) != trialSpaceBasis->gridView().size(0))
            DUNE_THROW(Dune::InvalidStateException, "Trial and Ansatz space must operate on the same grids!");

        // Check if the overlap size is what we expect
        if (!CheckOverlapSize<DiscretizationMethod::fem>::isValid(ansatzSpaceBasis->gridView()))
            DUNE_THROW(Dune::InvalidStateException, "The finite element discretization method only works with zero overlap for parallel computations. "
                                                     << " Set the parameter \"Grid.Overlap\" in the input file.");
    }

    //! Constructor for petrov Galerkin type settings
    template<bool isSG = standardGalerkin, std::enable_if_t<!isSG, int> = 0>
    FEGridGeometry(std::shared_ptr<AnsatzSpaceBasis> ansatzSpaceBasis,
                   std::shared_ptr<TrialSpaceBasis> trialSpaceBasis)
    : ParentType(ansatzSpaceBasis->gridView())
    , ansatzSpaceBasis_(ansatzSpaceBasis)
    , trialSpaceBasis_(trialSpaceBasis)
    {
        if (ansatzSpaceBasis->gridView().size(0) != trialSpaceBasis->gridView().size(0))
            DUNE_THROW(Dune::InvalidStateException, "Trial and Ansatz space must operate on the same grids!");

        // Check if the overlap size is what we expect
        if (!CheckOverlapSize<DiscretizationMethod::fem>::isValid(ansatzSpaceBasis->gridView()))
            DUNE_THROW(Dune::InvalidStateException, "The finite element discretization method only works with zero overlap for parallel computations. "
                                                     << " Set the parameter \"Grid.Overlap\" in the input file.");
    }

    //! Returns true if the undelying approach is of standard Galerkin type.
    static constexpr bool isStandardGalerkin()
    { return standardGalerkin; }

    //! update the mappers etc
    void update()
    {
        ParentType::update();
    }

    //! The total number of degrees of freedom
    auto numDofs() const
    { return ansatzSpaceBasis_->size(); }

    //! Return the ansatz space basis
    const AnsatzSpaceBasis& feBasis() const
    { return *ansatzSpaceBasis_; }

    //! Return the ansatz space basis
    const AnsatzSpaceBasis& ansatzSpaceBasis() const
    { return *ansatzSpaceBasis_; }

    //! Return the ansatz space basis
    const TrialSpaceBasis& trialSpaceBasis() const
    { return *trialSpaceBasis_; }

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
    std::shared_ptr<AnsatzSpaceBasis> ansatzSpaceBasis_;
    std::shared_ptr<TrialSpaceBasis> trialSpaceBasis_;
};

} // end namespace Dumux

#endif
