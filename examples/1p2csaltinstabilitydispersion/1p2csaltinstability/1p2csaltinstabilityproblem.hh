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
/**
 * \file
 * \brief Definition of the 1p2c salt instability problem
  *
 */
#ifndef DUMUX_SALTINSTABILITY1P2C_PROBLEM_HH
#define DUMUX_SALTINSTABILITY1P2C_PROBLEM_HH

#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/1pnc/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dune/grid/yaspgrid.hh>
#include "../watersaltfluidsystem.hh"
#include "1p2csaltinstabilityspatialparameters.hh"
#include "boxfickslawdispersion.hh"

namespace Dumux {

template <class TypeTag>
class SaltInstability1p2cProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct SaltInstability1p2cProblemTypeTag { using InheritsFrom = std::tuple<BoxModel, OnePNC>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::SaltInstability1p2cProblemTypeTag> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::SaltInstability1p2cProblemTypeTag> { using type = SaltInstability1p2cProblem<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::SaltInstability1p2cProblemTypeTag> { using type = WaterSaltFluidSystem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::SaltInstability1p2cProblemTypeTag>
{
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = SaltInstability1p2cSpatialParams<FVGridGeometry, Scalar>;
};

//Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::SaltInstability1p2cProblemTypeTag> { static constexpr bool value = false; };

template<class TypeTag>
struct MolecularDiffusionType<TypeTag, TTag::SaltInstability1p2cProblemTypeTag> { using type = FicksLawDispersionImplementation<TypeTag, DiscretizationMethod::box>; };

} // end namespace Properties

template <class TypeTag>
class SaltInstability1p2cProblem : public PorousMediumFlowProblem<TypeTag>
{
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementIterator = typename GridView::template Codim<0>::Iterator;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

     // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
      // indices of the equations
    static constexpr int conti0EqIdx = Indices::conti0EqIdx;

    // Grid and world dimension
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    SaltInstability1p2cProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
      FluidSystem::init();
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 36 degrees Celsius.
     */
    Scalar temperature() const
    {
        return 273.15 + 20; // in [K]
    };

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        const auto xMax = this->fvGridGeometry().bBoxMax()[0];
        const auto yMax = this->fvGridGeometry().bBoxMax()[1];

        values.setAllNeumann();

         if (globalPos[0]>xMax-eps_)
         {
         values.setAllDirichlet();
         }

         if (globalPos[1]>yMax-eps_ && globalPos[0]>0.5*yMax && globalPos[0]<1.5*yMax)
         {
           values.setAllDirichlet();
         }

         if (globalPos[0]<eps_)
         {
         values.setAllDirichlet();
         }

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        const auto yMax = this->fvGridGeometry().bBoxMax()[1];

        initial_(values, globalPos);
        if (globalPos[1]>yMax-eps_ && globalPos[0]>0.5*yMax && globalPos[0]<1.5*yMax)
        {
            values[FluidSystem::SaltIdx] =  0.5;
        }

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolumeFace& scvf) const
    {
        //const auto globalPos = element.geometry().corner(scvf.insideScvIdx());

        NumEqVector values(0.0);

        return values;
    }
    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * For this method, the \a values parameter stores the rate mass
     * of a component is generated or annihilate per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        return NumEqVector(0.0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;

        initial_(values, globalPos);

        return values;
    }
    // \}

private:
    // the internal method for the initial condition
    void initial_(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        Scalar densityW = 1025;
        values[Indices::pressureIdx] = 1.0133e5+(depthBOR_-globalPos[1])*densityW*9.81; //initial condition for the pressure
        values[FluidSystem::SaltIdx] = 0.0; //initial condition for the salt molefraction
    }

    static constexpr Scalar eps_ = 1e-6;
    static constexpr Scalar depthBOR_ = 1.0; // [m]
};

} // end namespace Dumux

#endif // DUMUX_HENRY1P2C_PROBLEM_HH
