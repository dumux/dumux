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
 * \ingroup TwoPNCTests
 * \brief Problem where air is injected under a low permeable layer in a depth of 2700m.
 */
#ifndef DUMUX_TWOPNC_DIFFUSION_PROBLEM_HH
#define DUMUX_TWOPNC_DIFFUSION_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/porousmediumflow/2pnc/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/fluidsystems/h2on2.hh>

#include "spatialparams.hh"
#include <dumux/flux/maxwellstefanslaw.hh>

#ifndef DIFFUSIONTYPE // default to Fick's law if not set through CMake
#define DIFFUSIONTYPE FicksLaw<TypeTag>
#endif

namespace Dumux {

template <class TypeTag>
class TwoPNCDiffusionProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct TwoPNCDiffusion { using InheritsFrom = std::tuple<TwoPNC>; };
struct TwoPNCDiffusionCC { using InheritsFrom = std::tuple<TwoPNCDiffusion, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPNCDiffusion> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPNCDiffusion> { using type = TwoPNCDiffusionProblem<TypeTag>; };

// // Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TwoPNCDiffusion>
{
    using type = FluidSystems::H2ON2<GetPropType<TypeTag, Properties::Scalar>,
                                     FluidSystems::H2ON2DefaultPolicy</*fastButSimplifiedRelations=*/true>>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TwoPNCDiffusion>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TwoPNCDiffusionSpatialParams<GridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::TwoPNCDiffusion> { static constexpr bool value = true; };

//! Here we set FicksLaw or TwoPNCDiffusionsLaw
template<class TypeTag>
struct MolecularDiffusionType<TypeTag, TTag::TwoPNCDiffusion> { using type = DIFFUSIONTYPE; };

//! Set the default formulation to pw-Sn: This can be over written in the problem.
template<class TypeTag>
struct Formulation<TypeTag, TTag::TwoPNCDiffusion>
{ static constexpr auto value = TwoPFormulation::p0s1; };

} // end namespace Properties


/*!
 * \ingroup TwoPNCTests
 * \brief Problem where air is injected under a low permeable layer in a depth of 2700m.
 */
template <class TypeTag>
class TwoPNCDiffusionProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;

    //! Property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();

public:
    TwoPNCDiffusionProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // initialize the tables of the fluid system
        FluidSystem::init();

        name_ = getParam<std::string>("Problem.Name");

        // stating in the console whether mole or mass fractions are used
        if(useMoles)
        {
            std::cout<<"problem uses mole-fractions"<<std::endl;
        }
        else
        {
            std::cout<<"problem uses mass-fractions"<<std::endl;
        }
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Returns the temperature [K]
     */
    Scalar temperature() const
    { return 293.15; }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllDirichlet();
        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
         PrimaryVariables priVars;
         priVars.setState(Indices::firstPhaseOnly);
         priVars[Indices::pressureIdx] = 1e5;
         priVars[Indices::switchIdx] = 1e-5 ;

         if (globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_)
             priVars[Indices::switchIdx] = 1e-3;

        return priVars;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean influx.
     * \param globalPos The global position
     *
     * \note The units must be according to either using mole or mass fractions (mole/(m^2*s) or kg/(m^2*s)).
     */
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector values(0.0);
        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial values for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        return initial_(globalPos);
    }

    // \}

private:
    /*!
     * \brief Evaluates the initial values for a control volume.
     *
     * The internal method for the initial condition
     *
     * \param globalPos The global position
     */
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars.setState(Indices::firstPhaseOnly);

        //mole-fraction formulation
        priVars[Indices::switchIdx] = 1e-5;
        priVars[Indices::pressureIdx] = 1e5;
        return priVars;
    }

    Scalar temperature_;
    static constexpr Scalar eps_ = 1e-6;

   std::string name_ ;


};

} // end namespace Dumux

#endif
