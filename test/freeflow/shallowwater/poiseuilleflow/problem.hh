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
 * \brief A test for the Shallow water model (flow around a square pillar).
 */
#ifndef DUMUX_POISEUILLE_FLOW_TEST_PROBLEM_HH
#define DUMUX_POISEUILLE_FLOW_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>
#include <dumux/discretization/cctpfa.hh>
#include "spatialparams.hh"

#include <dumux/freeflow/shallowwater/model.hh>
#include <dumux/freeflow/shallowwater/problem.hh>
#include <dumux/flux/shallowwater/riemannproblem.hh>
#include <dumux/flux/shallowwater/exactriemann.hh>


namespace Dumux {

/*!
 * \ingroup ShallowWaterTests
 * \brief A simple test for the flow around a square pillar using the shallow water equations
 */
template <class TypeTag>
class PoiseuilleFlowProblem;


// Specify the properties for the problem
namespace Properties {

// Create new type tags
namespace TTag {
struct PoiseuilleFlow { using InheritsFrom = std::tuple<ShallowWater, CCTpfaModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::PoiseuilleFlow>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PoiseuilleFlow>
{ using type = Dumux::PoiseuilleFlowProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PoiseuilleFlow>
{
private:
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = PoiseuilleFlowSpatialParams<FVGridGeometry, Scalar>;
};

template<class TypeTag>
struct EnableFVGridGeometryCache<TypeTag, TTag::PoiseuilleFlow>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::PoiseuilleFlow>
{ static constexpr bool value = false; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::PoiseuilleFlow
{ static constexpr bool value = false; };

} // end namespace Properties


/*!
 * \ingroup Shallow water equations model
 * \ingroup ImplicitTestProblems
 *
 * \brief A simple test for the 2D flow through a channel with rough walls (Poiseuille flow).
 *
 * The domain has a length L =  400 meters long and a width W = 100 meters.
 * The bed level is sloped from z = -9.98 (x=0) to z = -10.0 meters (x=L)
 * The initital water depth corresponds to the analytical solution:
 * having a slope equal to ib = dh / L, where dh = -0.02 m.
 * At the west/left  (inflow)  boundary a discharge is prescribed of Q_in = -4087.5 m^3/s.
 * At the east/right (outflow) boundary a fixed water level is prescribed of theta = 0.0 m.
 * The south and north boundaries are set to roughwall type boundaries,
 * with a coefficient alphaW = 1, where:
*      alphaW = 0.0: full    slip (smooth wall)
* 0.0 <alphaW < 1.0: partial slip (partially-rough wall)
*      alphaW = 1.0: no      slip (fully-rough wall)
 * Additionally these (south and north) boundaries are (automatically) set to no-flow boundaries.
  *
 * This problem uses the \ref ShallowWaterModel
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_shallowwater -parameterFile test_shallowwater.input -TimeManager.TEnd 12000</tt>
 *
 * where the initial time step is 1.0 seconds, and the end time of the
 * simulation is 12000 seconds
 */
template <class TypeTag>
class PoiseuilleFlowProblem : public ShallowWaterProblem<TypeTag>
{
    using ParentType = ShallowWaterProblem<TypeTag>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
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
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    PoiseuilleFlowProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        name_ = getParam<std::string>("Problem.Name");
        exactWaterDepth_.resize(fvGridGeometry->numDofs(), 0.0);
        exactVelocityX_.resize(fvGridGeometry->numDofs(), 0.0);
        exactVelocityY_.resize(fvGridGeometry->numDofs(), 0.0);
    }

    //! Get the analytical water depth
    const std::vector<Scalar>& getExactWaterDepth()
    {
        return exactWaterDepth_;
    }

    //! Get the analytical U-velocity
    const std::vector<Scalar>& getExactVelocityX()
    {
        return exactVelocityX_;
    }

    //! Get the analytical V-velocity
    const std::vector<Scalar>& getExactVelocityY()
    {
        return exactVelocityY_;
    }

    //! Udpate the analytical solution
    template<class SolutionVector, class GridVariables>
    void updateAnalyticalSolution(const SolutionVector& curSol,
                                  const GridVariables& gridVariables,
                                  const Scalar time)
    {
        //compute solution for all elements
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, curSol);

            const auto& globalPos = element.geometry().center();

            //compute the position s and the initial water depth at the gate, velocities are zero
            const Scalar dh = -0.02;
            const Scalar waterDepthLeft =  initialWaterDepthLeft_;
            const Scalar waterDepthRight =  initialWaterDepthRight_;
            const auto gravity = this->spatialParams().gravity(globalPos);

            auto riemannResult = ShallowWater::exactPoiseuille(Qin,
                                                            dh,
                                                            0.0,
                                                            0.0,
                                                            0.0,
                                                            0.0,
                                                            gravity,
                                                            s);

            const auto eIdx = this->fvGridGeometry().elementMapper().index(element);
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
