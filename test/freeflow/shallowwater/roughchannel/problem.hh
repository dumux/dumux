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
 * \brief A test for the Shallow water model (rough channel).
 */
#ifndef DUMUX_ROUGH_CHANNEL_TEST_PROBLEM_HH
#define DUMUX_ROUGH_CHANNEL_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>
#include <dumux/discretization/cctpfa.hh>
#include "spatialparams.hh"
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>

#include <dumux/freeflow/shallowwater/model.hh>
#include <dumux/freeflow/shallowwater/problem.hh>
#include <dumux/freeflow/shallowwater/boundaryfluxes.hh>

namespace Dumux {

template <class TypeTag>
class RoughChannelProblem;

// Specify the properties for the problem
namespace Properties {

// Create new type tags
namespace TTag {
struct RoughChannel { using InheritsFrom = std::tuple<ShallowWater, CCTpfaModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::RoughChannel>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::RoughChannel>
{ using type = Dumux::RoughChannelProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RoughChannel>
{
private:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;
public:
    using type = RoughChannelSpatialParams<GridGeometry, Scalar, VolumeVariables>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::RoughChannel>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::RoughChannel>
{ static constexpr bool value = false; };
} // end namespace Properties

/*!
 * \ingroup ShallowWaterTests
 * \brief A simple flow in a rough channel with friction law after Manning.
 *
 * The domain is 500 meters long and 5 meters wide. At the left border a discharge
 * boundary condition is applied and at the right border a water depth boundary condition.
 * All other boundaries are set to no-flow. Normal flow is assumed, therefor the water depth
 * at the right border can be calculated with the formular of Gaukler-Manning-Strickler.
 *
 * \f[
 * v_m = 1/n * R_{hy}^{2/3} * I_s^{1/2}
 * \f]
 *
 * With the mean velocity
 * \f[
 * v_m = \frac{q}/{h}
 * \f]
 * the friction value n after Manning
 * the hydraulic radius R_{hy} equal to the water depth h (because normal flow is assumed)
 * the bed slope I_s and the unity inflow discharge q.
 *
 * Therefore h can be calculated with
 *
 * \f[
 * h = \left(\frac{n*q}{\sqrt{I_s}} \right)^{3/5}
 * \f]
 *
 * The formula of Gaukler Manning and Strickler is also used to calculate the analytic solution.
 *
 * This problem uses the \ref ShallowWaterModel
 */
template <class TypeTag>
class RoughChannelProblem : public ShallowWaterProblem<TypeTag>
{
    using ParentType = ShallowWaterProblem<TypeTag>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
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
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

public:
    RoughChannelProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        name_ = getParam<std::string>("Problem.Name");
        exactWaterDepth_.resize(gridGeometry->numDofs(), 0.0);
        exactVelocityX_.resize(gridGeometry->numDofs(), 0.0);
        constManningN_ = getParam<Scalar>("Problem.ManningN");
        bedSlope_ = getParam<Scalar>("Problem.BedSlope");
        discharge_ = getParam<Scalar>("Problem.Discharge");
        hBoundary_ = this->gauklerManningStrickler(discharge_,constManningN_,bedSlope_);
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

    //! Calculate the water depth with Gaukler-Manning-Strickler
    Scalar gauklerManningStrickler(Scalar discharge, Scalar manningN, Scalar bedSlope)
    {
        using std::pow;
        using std::abs;
        using std::sqrt;

        return pow(abs(discharge)*manningN/sqrt(bedSlope), 0.6);
    }

    //! Udpate the analytical solution
    void updateAnalyticalSolution()
    {
        using std::abs;

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const Scalar h = this->gauklerManningStrickler(discharge_,constManningN_,bedSlope_);
            const Scalar u = abs(discharge_)/h;

            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            exactWaterDepth_[eIdx] = h;
            exactVelocityX_[eIdx] = u;
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

        source += bottomFrictionSource(element, fvGeometry, elemVolVars, scv);

        return source;
    }

    /*!
     * \brief Compute the source term due to bottom friction
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub control volume
     *
     * \return source
     */
     NumEqVector bottomFrictionSource(const Element& element,
                                      const FVElementGeometry& fvGeometry,
                                      const ElementVolumeVariables& elemVolVars,
                                      const SubControlVolume &scv) const
     {
        NumEqVector bottomFrictionSource(0.0);

        const auto& volVars = elemVolVars[scv];
        Dune::FieldVector<Scalar, 2> bottomShearStress = this->spatialParams().frictionLaw(element, scv).shearStress(volVars);

        bottomFrictionSource[0] = 0.0;
        bottomFrictionSource[1] =bottomShearStress[0];
        bottomFrictionSource[2] =bottomShearStress[1];

        return bottomFrictionSource;
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

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& nxy = scvf.unitOuterNormal();
        const auto gravity = this->spatialParams().gravity(scvf.center());
        std::array<Scalar, 3> boundaryStateVariables;

        // impose discharge at the left side
        if (scvf.center()[0] < this->gridGeometry().bBoxMin()[0] + eps_)
        {
            boundaryStateVariables = ShallowWater::fixedDischargeBoundary(discharge_,
                                                                          insideVolVars.waterDepth(),
                                                                          insideVolVars.velocity(0),
                                                                          insideVolVars.velocity(1),
                                                                          gravity,
                                                                          nxy);
        }
        // impose water depth at the right side
        else if (scvf.center()[0] > this->gridGeometry().bBoxMax()[0] - eps_)
        {
            boundaryStateVariables =  ShallowWater::fixedWaterDepthBoundary(hBoundary_,
                                                                            insideVolVars.waterDepth(),
                                                                            insideVolVars.velocity(0),
                                                                            insideVolVars.velocity(1),
                                                                            gravity,
                                                                            nxy);
        }
        // no flow boundary
        else
        {
            boundaryStateVariables[0] = insideVolVars.waterDepth();
            boundaryStateVariables[1] = -insideVolVars.velocity(0);
            boundaryStateVariables[2] = -insideVolVars.velocity(1);
        }

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
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        values[0] = hBoundary_;
        values[1] = abs(discharge_)/hBoundary_;
        values[2] = 0.0;

        return values;
    };

    // \}

private:
    std::vector<Scalar> exactWaterDepth_;
    std::vector<Scalar> exactVelocityX_;
    Scalar hBoundary_;
    Scalar constManningN_; // analytic solution is only available for const friction.
    Scalar bedSlope_;
    Scalar discharge_; // discharge at the inflow boundary
    static constexpr Scalar eps_ = 1.0e-6;
    std::string name_;
};

} //end namespace Dumux

#endif
