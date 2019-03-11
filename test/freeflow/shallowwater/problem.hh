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
 * \ingroup ShallowWaterTest
 * \brief Test for the Shallow water model (wet dam break).
 */
#ifndef DUMUX_DAM_BREAK_TEST_PROBLEM_HH
#define DUMUX_DAM_BREAK_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>
#include <dumux/discretization/cctpfa.hh>
#include "spatialparams.hh"

#include <dumux/freeflow/shallowwater/model.hh>
#include <dumux/freeflow/shallowwater/problem.hh>
#include <dumux/flux/shallowwater/riemannproblem.hh>


namespace Dumux
{
/*!
 * \ingroup ShallowWaterTests
 * \brief A simple dam break test for the shallow water equations
 */
template <class TypeTag>
class DamBreakProblem;


// Specify the properties for the problem
namespace Properties{

// Create new type tags
namespace TTag
{
struct ShallowWaterModel{ using InheritsFrom = std::tuple<ShallowWater>; };
struct DamBreakWet{ using InheritsFrom = std::tuple<ShallowWaterModel, CCTpfaModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::ShallowWaterModel>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::ShallowWaterModel>{ using type = Dumux::DamBreakProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::DamBreakWet>
{
private:
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = DamBreakSpatialParams<FVGridGeometry, Scalar>;
};

template<class TypeTag>
struct EnableFVGridGeometryCache<TypeTag, TTag::ShallowWaterModel> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::ShallowWaterModel> { static constexpr bool value = false; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::ShallowWaterModel> { static constexpr bool value = false; };


} // end namespace Properties



/*!
 * \ingroup Shallow water equations model
 * \ingroup ImplicitTestProblems
 *
 * \brief A simple dam break test
 *
 * This problem uses the \ref ShallowWaterModel
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_shallowwater -parameterFile test_shallowwater.input -TimeManager.TEnd 10</tt>
 *
 * where the initial time step is 0.01 seconds, and the end of the
 * simulation time is 10 seconds
 */
template <class TypeTag>
class DamBreakProblem : public ShallowWaterProblem<TypeTag>
{
    using ParentType = ShallowWaterProblem<TypeTag>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using PrimaryVariables = typename VolumeVariables::PrimaryVariables;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using NeumannFluxes = GetPropType<TypeTag, Properties::NumEqVector>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;


    enum {
        // copy some indices for convenience
        massBalanceIdx = Indices::massBalanceIdx,
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx,

        // Grid and world dimension
        dimWorld = FVGridGeometry::GridView::dimensionworld,

    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    DamBreakProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        name_ = getParam<std::string>("Problem.Name");
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }


    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the boundary type is set
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann();
        return bcTypes;
    }

    /*!
     * \brief Specifies the neumann bounday
     * \param element
     * \param fvGeometry
     * \param elemVolVars
     * \param scvf
     */

    NeumannFluxes neumann(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolVars,
                          const SubControlVolumeFace& scvf) const
    {
        NeumannFluxes values(0.0);

        //we need the Riemann invariants to compute the values depending of the boundary type
        //since we use a weak imposition we do not have a dirichlet value. We impose fluxes
        //based on q,h, etc. computed with the Riemann invariants

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& nxy = scvf.unitOuterNormal();
        const auto& ip = scvf.ipGlobal();


        auto riemannFlux = ShallowWater::riemannProblem(insideVolVars.waterDepth(),
                                                        insideVolVars.waterDepth(),
                                                        insideVolVars.velocity(0),
                                                        -insideVolVars.velocity(0),
                                                        insideVolVars.velocity(1),
                                                        -insideVolVars.velocity(1),
                                                        insideVolVars.bedSurface(),
                                                        insideVolVars.bedSurface(),
                                                        insideVolVars.gravity(),
                                                        nxy);

        values[massBalanceIdx] = riemannFlux[0];
        values[velocityXIdx]   = riemannFlux[1];
        values[velocityYIdx]   = riemannFlux[2];

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet boundary
     *        segment. For the shallow water equations we do a weak
     *        imposition of boundary conditions. You may use this for
     *        supercritical flow.
     *
     * \param globalPos The position for which the Dirichlet value is set
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {

        PrimaryVariables values(0.0);

        return values;

    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial values for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     *
     * \param globalPos The position for which the boundary type is set
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {

        PrimaryVariables values(0.0);

        values[0] = 1.0;
        values[1] = 0.0;
        values[2] = 0.0;

        // water level on the left side of the gate
        if (globalPos[0] < 10.0 + eps_)
        {
            values[0] = 4.0;
        }

        //water level on the right side of the gate
        else
        {
            values[0] = 1.0;
        }

        return values;

    };

    // \}

private:


    static constexpr Scalar eps_ = 1.0e-6;
    std::string name_;
};

} //end namespace Dumux

#endif
