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
 * \ingroup ShallowWaterTests
 * \brief A test for the Shallow water model (wet dam break).
 */
#ifndef DUMUX_DAM_BREAK_TEST_PROBLEM_HH
#define DUMUX_DAM_BREAK_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>
#include <dumux/discretization/cctpfa.hh>
#include "spatialparams.hh"

#include <dumux/common/boundarytypes.hh>
#include <dumux/freeflow/shallowwater/model.hh>
#include <dumux/freeflow/shallowwater/problem.hh>
#include <dumux/flux/shallowwater/riemannproblem.hh>
#include <dumux/flux/shallowwater/exactriemann.hh>

namespace Dumux {

/*!
 * \ingroup ShallowWaterTests
 * \brief A simple dam break test for the shallow water equations
 */
template <class TypeTag>
class DamBreakProblem;


// Specify the properties for the problem
namespace Properties {

// Create new type tags
namespace TTag {
struct DamBreakWet { using InheritsFrom = std::tuple<ShallowWater, CCTpfaModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::DamBreakWet>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DamBreakWet>
{ using type = Dumux::DamBreakProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::DamBreakWet>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = DamBreakSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::DamBreakWet>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::DamBreakWet>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::DamBreakWet>
{ static constexpr bool value = false; };

} // end namespace Properties


/*!
 * \ingroup Shallow water equations model
 * \ingroup ImplicitTestProblems
 *
 * \brief A simple dam break test (1D wet dam break).
 *
 * The domain is 20 meters long with a gate in the middle. On the left
 * side the water depth is 4 meters and on the right side the depth is 1 meter.
 * All boundaries are set to no-flow.
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
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using NeumannFluxes = GetPropType<TypeTag, Properties::NumEqVector>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    DamBreakProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        name_ = getParam<std::string>("Problem.Name");
        exactWaterDepth_.resize(gridGeometry->numDofs(), 0.0);
        exactVelocityX_.resize(gridGeometry->numDofs(), 0.0);
    }

    //! Get the analytical water depth
    const std::vector<Scalar>& getExactWaterDepth()
    {
        return exactWaterDepth_;
    }

    //! Get the analytical velocity
    const std::vector<Scalar>& getExactVelocityX()
    {
        return exactVelocityX_;
    }

    //! Udpate the analytical solution
    template<class SolutionVector, class GridVariables>
    void updateAnalyticalSolution(const SolutionVector& curSol,
                                  const GridVariables& gridVariables,
                                  const Scalar time)
    {
        //compute solution for all elements
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, curSol);

            const auto& globalPos = element.geometry().center();

            //compute the position s and the initial water depth at the gate, velocities are zero
            const Scalar s = (globalPos[0] - gatePosition_)/time;
            const Scalar waterDepthLeft =  initialWaterDepthLeft_;
            const Scalar waterDepthRight =  initialWaterDepthRight_;
            const auto gravity = this->spatialParams().gravity(globalPos);

            auto riemannResult = ShallowWater::exactRiemann(waterDepthLeft,
                                                            waterDepthRight,
                                                            0.0,
                                                            0.0,
                                                            0.0,
                                                            0.0,
                                                            gravity,
                                                            s);

            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            exactWaterDepth_[eIdx] = riemannResult.waterDepth;
            exactVelocityX_[eIdx] = riemannResult.velocityX;
        }
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
     * \param elemFluxVarsCache
     * \param scvf
     */
    NeumannFluxes neumann(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolVars,
                          const ElementFluxVariablesCache& elemFluxVarsCache,
                          const SubControlVolumeFace& scvf) const
    {
        NeumannFluxes values(0.0);

        // we need the Riemann invariants to compute the values depending of the boundary type
        // since we use a weak imposition we do not have a dirichlet value. We impose fluxes
        // based on q,h, etc. computed with the Riemann invariants
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& nxy = scvf.unitOuterNormal();
        const auto gravity = this->spatialParams().gravity(scvf.center());

        auto riemannFlux = ShallowWater::riemannProblem(insideVolVars.waterDepth(),
                                                        insideVolVars.waterDepth(),
                                                        insideVolVars.velocity(0),
                                                        -insideVolVars.velocity(0),
                                                        insideVolVars.velocity(1),
                                                        -insideVolVars.velocity(1),
                                                        insideVolVars.bedSurface(),
                                                        insideVolVars.bedSurface(),
                                                        gravity,
                                                        nxy);

        values[Indices::massBalanceIdx] = riemannFlux[0];
        values[Indices::velocityXIdx]   = riemannFlux[1];
        values[Indices::velocityYIdx]   = riemannFlux[2];

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

        values[0] = initialWaterDepthRight_;
        values[1] = 0.0;
        values[2] = 0.0;

        // water level on the left side of the gate
        if (globalPos[0] < 10.0 + eps_)
        {
            values[0] = initialWaterDepthLeft_;
        }

        //water level on the right side of the gate
        else
        {
            values[0] = initialWaterDepthRight_;
        }

        return values;
    };

    // \}

private:

    std::vector<Scalar> exactWaterDepth_;
    std::vector<Scalar> exactVelocityX_;

    static constexpr Scalar initialWaterDepthLeft_ = 4.0;
    static constexpr Scalar initialWaterDepthRight_ = 1.0;
    static constexpr Scalar channelLenght_ = 20.0;
    static constexpr Scalar gatePosition_ = 10.0;

    static constexpr Scalar eps_ = 1.0e-6;
    std::string name_;
};

} //end namespace Dumux

#endif
