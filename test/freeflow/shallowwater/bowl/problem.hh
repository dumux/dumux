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
 * \brief A test for the Shallow water model (bowl).
 */
#ifndef DUMUX_BOWL_TEST_PROBLEM_HH
#define DUMUX_BOWL_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>
#include <dumux/discretization/cctpfa.hh>
#include "spatialparams.hh"
#include <dumux/common/parameters.hh>

#include <dumux/freeflow/shallowwater/model.hh>
#include <dumux/freeflow/shallowwater/problem.hh>
#include <dumux/freeflow/shallowwater/boundaryfluxes.hh>

namespace Dumux {

template <class TypeTag>
class BowlProblem;

// Specify the properties for the problem
namespace Properties {

// Create new type tags
namespace TTag {
struct Bowl { using InheritsFrom = std::tuple<ShallowWater, CCTpfaModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::Bowl>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Bowl>
{ using type = Dumux::BowlProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Bowl>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;
public:
    using type = BowlSpatialParams<GridGeometry, Scalar, VolumeVariables>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Bowl>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Bowl>
{ static constexpr bool value = false; };
} // end namespace Properties

/*!
 * \ingroup ShallowWaterTests
 * \brief A wetting and drying test with sloshing water in a bowl.
 *
 * The domain is 4 meters long and 4 meters wide. The center of the domain is loacted at
 * x = 0 and y = 0. There is no flow over the boundaries and no friction is considered.
 *
 * This example is demanding for the implicit model if a high mesh resolution is applied
 * (e.g. 150x150 cells) in combination with a large time step size. Using the new limiting
 * (UpwindFluxLimiting = true) will help to improve the convergence for such cases.
 *
 * This test uses a low mesh resoultion and only ensures that UpwindFluxLimiting for the mobility
 * works.
 *
 * The results are checked against a analytical solution which is based on the "Thacker-Solution"
 * (William Thacker, "Some exact solutions to the nonlinear shallow-water wave equations", Journal
 * of Fluid Mechanics, 107:499–508, 1981). Further examples and details on the solution are given
 * in SWASHES (Shallow Water Analytic Solutions for Hydraulic and Environmental Studies,
 * https://www.idpoisson.fr/swashes/).
 *
 * This problem uses the \ref ShallowWaterModel
 */
template <class TypeTag>
class BowlProblem : public ShallowWaterProblem<TypeTag>
{
    using ParentType = ShallowWaterProblem<TypeTag>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using NeumannFluxes = GetPropType<TypeTag, Properties::NumEqVector>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

public:
    BowlProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        using std::sqrt;
        using std::pow;

        name_ = getParam<std::string>("Problem.Name");
        gravity_ = getParam<Scalar>("Problem.Gravity");
        bowlDepthAtCenter_ =  getParam<Scalar>("Problem.BowlDepthAtCenter");
        bowlParaboloidRadius_ =  getParam<Scalar>("Problem.BowlParaboloidRadius");
        bowlInitialWaterElevationAtCenter_ =  getParam<Scalar>("Problem.BowlInitialWaterElevationAtCenter");
        bowlAnalyticParameterOmega_ = sqrt(8.0 * gravity_ * bowlDepthAtCenter_) / bowlParaboloidRadius_;
        bowlAnalyticParameterA_ = (pow(bowlDepthAtCenter_ + bowlInitialWaterElevationAtCenter_, 2.0) -
                                  pow(bowlDepthAtCenter_, 2.0)) /
                                  (pow(bowlDepthAtCenter_ + bowlInitialWaterElevationAtCenter_, 2.0)+
                                  pow(bowlDepthAtCenter_, 2.0));
        exactWaterDepth_.resize(gridGeometry->numDofs(), 0.0);
        updateAnalyticalSolution(0.0);
    }

    //! Get the analytical water depth
    const std::vector<Scalar>& getExactWaterDepth()
    {
        return exactWaterDepth_;
    }

    //! Get the analctic water depth at
    Scalar calculateAnalyticWaterDepth(Scalar time, Scalar x, Scalar y)
    {
        using std::max;
        using std::sqrt;
        using std::pow;
        using std::cos;

        auto waterDepth = max(bowlDepthAtCenter_* ((sqrt(1.0 - pow(bowlAnalyticParameterA_, 2.0)))/
                          (1.0 - bowlAnalyticParameterA_ * cos(bowlAnalyticParameterOmega_ * time))-
                          1.0 - (pow(x, 2.0) + pow(y, 2.0)) / pow(bowlParaboloidRadius_, 2.0) *
                          ((1.0- pow(bowlAnalyticParameterA_, 2.0)) /
                          (pow(1.0 - bowlAnalyticParameterA_ *
                          cos(bowlAnalyticParameterOmega_ * time), 2.0)) - 1.0)) -
                          bowlDepthAtCenter_ *((pow(x, 2.0) + pow(y, 2.0))/
                          pow(bowlParaboloidRadius_, 2.0) - 1.0), 1.0E-5);

        return waterDepth;
    }

    //! Compute L2 error
    void computeL2error(const Scalar time,
                        const SolutionVector& curSol,
                        const GridVariables& gridVariables)
    {
        Scalar l2error = 0.0;

        //first ensure that the Analytical solution is up-to-date
        updateAnalyticalSolution(time);

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            const auto& globalPos = element.geometry().center();
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);
            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, curSol);
            exactWaterDepth_[eIdx] = calculateAnalyticWaterDepth(time, globalPos[0], globalPos[1]);

            for (auto&& scv : scvs(fvGeometry))
            {
                using std::pow;
                l2error += pow(exactWaterDepth_[eIdx] - elemVolVars[scv].waterDepth(), 2.0);
            }
        }
        using std::sqrt;
        l2error = sqrt(l2error);
        l2error = this->gridGeometry().gridView().comm().sum(l2error);

        if (this->gridGeometry().gridView().comm().rank() == 0)
        {
            std::cout << "L2 error at t =  "
                  <<  time << " seconds "
                  << " for "
                  << std::setw(6) << this->gridGeometry().gridView().size(0)
                  << " elements: "
                  << std::scientific
                  << l2error
                  << std::endl;
        }
    }

    //! Udpate the analytical solution
    void updateAnalyticalSolution(const Scalar time)
    {
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            const auto& globalPos = element.geometry().center();

            exactWaterDepth_[eIdx] = calculateAnalyticWaterDepth(time, globalPos[0], globalPos[1]);
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
    {
        return name_;
    }

     /*!
     * \brief Evaluate the source term for all balance equations within a given
     *        sub-control-volume.
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub control volume
     *
     * For this method, the \a values parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
     NumEqVector source(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolume &scv) const
    {
        NumEqVector source (0.0);

        return source;
    }

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
     * \brief Specifies the neumann boundary
     *
     *  We need the Riemann invariants to compute the values depending of the boundary type.
     *  Since we use a weak imposition we do not have a dirichlet value. We impose fluxes
     *  based on q, h, etc. computed with the Riemann invariants
     *
     * \param element
     * \param fvGeometry
     * \param elemVolVars
     * \param scvf
     */
    NeumannFluxes neumann(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolVars,
                          const ElementFluxVariablesCache& elemFluxVarsCache,
                          const SubControlVolumeFace& scvf) const
    {
        NeumannFluxes values(0.0);

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& nxy = scvf.unitOuterNormal();
        const auto gravity = this->spatialParams().gravity(scvf.center());
        std::array<Scalar, 3> boundaryStateVariables;

        //no flow with zero normal velocity and tangential velocity
        const auto vNormalGhost = -(nxy[0] * insideVolVars.velocity(0) +  nxy[1] * insideVolVars.velocity(1));
        const auto vTangentialGhost = -nxy[1] * insideVolVars.velocity(0) + nxy[0] * insideVolVars.velocity(1);

        boundaryStateVariables[0] = insideVolVars.waterDepth();
        boundaryStateVariables[1] =  nxy[0] * vNormalGhost - nxy[1] * vTangentialGhost;
        boundaryStateVariables[2] =  nxy[1] * vNormalGhost + nxy[0] * vTangentialGhost;

        auto riemannFlux = ShallowWater::riemannProblem(insideVolVars.waterDepth(),
                                                        boundaryStateVariables[0],
                                                        insideVolVars.velocity(0),
                                                        boundaryStateVariables[1],
                                                        insideVolVars.velocity(1),
                                                        boundaryStateVariables[2],
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
    PrimaryVariables initial(const Element& element) const
    {
        PrimaryVariables values(0.0);

        auto elemId = this->gridGeometry().elementMapper().index(element);
        values[0] = exactWaterDepth_[elemId];
        values[1] = 0.0;
        values[2] = 0.0;

        return values;
    };

    // \}

private:
    std::vector<Scalar> exactWaterDepth_;
    Scalar gravity_;
    Scalar bowlDepthAtCenter_;
    Scalar bowlParaboloidRadius_;
    Scalar bowlInitialWaterElevationAtCenter_;
    Scalar bowlAnalyticParameterOmega_;
    Scalar bowlAnalyticParameterA_;
    static constexpr Scalar eps_ = 1.0e-6;
    std::string name_;
};

} //end namespace Dumux

#endif
