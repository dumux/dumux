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
 * \ingroup TracerTests
 * \brief The properties for the incompressible test
 */
#ifndef DUMUX_INCOMPRESSIBLE_ONEP_TEST_PROBLEM_HH
#define DUMUX_INCOMPRESSIBLE_ONEP_TEST_PROBLEM_HH

#include <dumux/discretization/cellcentered/mpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/omethod/staticinteractionvolume.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include "1ptestspatialparams.hh"

namespace Dumux
{
/*!
 * \ingroup TracerTests
 * \brief The properties for the incompressible test
 */
// forward declarations
template<class TypeTag>
class OnePTestProblem;

namespace Properties
{
NEW_TYPE_TAG(IncompressibleTestProblem, INHERITS_FROM(CCMpfaModel, OneP));

// Set the grid type
SET_TYPE_PROP(IncompressibleTestProblem, Grid, Dune::UGGrid<3>);

// Set the problem type
SET_TYPE_PROP(IncompressibleTestProblem, Problem, OnePTestProblem<TypeTag>);
SET_TYPE_PROP(IncompressibleTestProblem, SpatialParams, OnePTestSpatialParams<TypeTag>);
SET_TYPE_PROP(IncompressibleTestProblem, LocalResidual, OnePIncompressibleLocalResidual<TypeTag>);

// the fluid system
SET_PROP(IncompressibleTestProblem, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::LiquidPhase<Scalar, Dumux::Components::Constant<0, Scalar> >;
};

// Enable caching
SET_BOOL_PROP(IncompressibleTestProblem, EnableGridVolumeVariablesCache, false);
SET_BOOL_PROP(IncompressibleTestProblem, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(IncompressibleTestProblem, EnableFVGridGeometryCache, true);

SET_BOOL_PROP(IncompressibleTestProblem, SolutionDependentAdvection, false);

// use the static interaction volume type with sizes known at compile time
SET_PROP(IncompressibleTestProblem, PrimaryInteractionVolume)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using NodalIndexSet = typename GET_PROP_TYPE(TypeTag, DualGridNodalIndexSet);

    // use the default traits
    using Traits = CCMpfaODefaultStaticInteractionVolumeTraits< NodalIndexSet, Scalar, 8, 12 >;
public:
    using type = CCMpfaOStaticInteractionVolume< Traits >;
};

// sizes of the ivs are known, use Dune::ReservedVector in index sets
SET_PROP(IncompressibleTestProblem, DualGridNodalIndexSet)
{
    using GV = typename GET_PROP_TYPE(TypeTag, GridView);
    static constexpr int dim = GV::dimension;
    static constexpr int dimWorld = GV::dimensionworld;
private:
    struct Traits
    {
        using GridView = GV;
        using GridIndexType = typename GV::IndexSet::IndexType;
        using LocalIndexType = std::uint8_t;

        //! per default, we use dynamic data containers (iv size unknown)
        template< class T > using NodalScvDataStorage = Dune::ReservedVector< T, 8 >;
        template< class T > using NodalScvfDataStorage = Dune::ReservedVector< T, 24 >;

        //! store data on neighbors of scvfs in static containers if possible
        template< class T >
        using ScvfNeighborDataStorage = typename std::conditional_t< (dim<dimWorld),
                                                                     std::vector< T >,
                                                                     Dune::ReservedVector< T, 2 > >;
    };
public:
    using type = CCMpfaDualGridNodalIndexSet< Traits >;
};

} // end namespace Properties

template<class TypeTag>
class OnePTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);

    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;

public:
    OnePTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry) : ParentType(fvGridGeometry)
    {
        headAtInlet_ = getParam<Scalar>("Problem.HeadAtInlet");
    }

    //! Specifies which kind of boundary condition should be used
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        if (isOnInflowBoundary(globalPos) || isOnOutflowBoundary(globalPos))
            values.setAllDirichlet();
        return values;
    }

    //! Evaluate the boundary conditions for a dirichlet control volume.
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        if (isOnInflowBoundary(globalPos))
            return PrimaryVariables(headAtInlet_);
        else
            return PrimaryVariables(1.0);
    }

    //! Returns the temperature in the domain
    Scalar temperature() const { return 283.15; /* 10°C*/ }

    //! Returns if a given position is on inflow boundary
    bool isOnInflowBoundary(const GlobalPosition& globalPos) const { return globalPos[0] < 1.0e-6 && globalPos[2] > 90.0; }

    //! Returns if a given position is on outflow boundary
    bool isOnOutflowBoundary(const GlobalPosition& globalPos) const { return globalPos[1] < 1.0e-6 && globalPos[2] < 10.0; }

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    template<class VTKWriter>
    void addVtkFields(VTKWriter& vtk)
    {
        const auto& gridView = this->fvGridGeometry().gridView();
        kField_.resize(gridView.size(0));
        poro_.resize(gridView.size(0));

        for (const auto& e : elements(this->fvGridGeometry().gridView()))
        {
            const auto eIdx = this->fvGridGeometry().elementMapper().index(e);
            kField_[eIdx] = this->spatialParams().permeability(e);
            poro_[eIdx] = this->spatialParams().porosity(e);
        }

        vtk.addCellData(kField_, "Permeability");
        vtk.addCellData(poro_, "Porosity");
    }

private:
    Scalar headAtInlet_;
    std::vector<Scalar> kField_;
    std::vector<Scalar> poro_;
};

} // end namespace Dumux

#endif
