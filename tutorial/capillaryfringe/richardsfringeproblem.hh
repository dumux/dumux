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
 *
 * \brief A water infiltration problem with a low-permeability lens
 *        embedded into a high-permeability domain which uses the
 *        Richards box model.
 */
#ifndef DUMUX_RICHARDS_LENSPROBLEM_HH
#define DUMUX_RICHARDS_LENSPROBLEM_HH

#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/porousmediumflow/richards/implicit/model.hh>

#include "fringespatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class RichardsFringeProblem;


// Specify the properties for the lens problem
namespace Properties
{
NEW_TYPE_TAG(RichardsFringeProblem, INHERITS_FROM(Richards));
NEW_TYPE_TAG(RichardsFringeBoxProblem, INHERITS_FROM(BoxModel, RichardsFringeProblem));
NEW_TYPE_TAG(RichardsFringeCCProblem, INHERITS_FROM(CCTpfaModel, RichardsFringeProblem));

// Set the spatial parameters
SET_TYPE_PROP(RichardsFringeProblem, SpatialParams, FringeSpatialParams<TypeTag>);

// Set the material law
SET_PROP(RichardsFringeProblem, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    // define the material law parameterized by absolute saturations
    using type = EffToAbsLaw<RegularizedVanGenuchten<Scalar>>;
};

// Use 2d YaspGrid
SET_TYPE_PROP(RichardsFringeProblem, Grid, Dune::YaspGrid<2>);

// Set the physical problem to be solved
SET_TYPE_PROP(RichardsFringeProblem, Problem, RichardsFringeProblem<TypeTag>);

// Enable gravity
SET_BOOL_PROP(RichardsFringeProblem, ProblemEnableGravity, true);
}

/*!
 * \ingroup RichardsModel
 * \ingroup ImplicitTestProblems
 *
 * \brief A water infiltration problem with a low-permeability lens
 *        embedded into a high-permeability domain which uses the
 *        Richards model.
 *
 * The domain is box shaped. Left and right boundaries are Dirichlet
 * boundaries with fixed water pressure (fixed Saturation \f$S_w = 0\f$),
 * bottom boundary is closed (Neumann 0 boundary), the top boundary
 * (Neumann 0 boundary) is also closed except for infiltration
 * section, where water is infiltrating into an initially unsaturated
 * porous medium. This problem is very similar the the LensProblem
 * which uses the TwoPBoxModel, with the main difference being that
 * the domain is initally fully saturated by gas instead of water and
 * water instead of a %DNAPL infiltrates from the top.
 *
 * This problem uses the \ref RichardsModel
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_boxrichards -parameterFile test_boxrichards.input -TimeManager.TEnd 10000000</tt>
 * <tt>./test_ccrichards -parameterFile test_ccrichards.input -TimeManager.TEnd 10000000</tt>
 *
 * where the initial time step is 100 seconds, and the end of the
 * simulation time is 10,000,000 seconds (115.7 days)
 */
template <class TypeTag>
class RichardsFringeProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum {
        // copy some indices for convenience
        pressureIdx = Indices::pressureIdx,
        conti0EqIdx = Indices::conti0EqIdx,

        // Grid and world dimension
        dimWorld = GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    /*!
     * \brief Constructor
     *
     * \param timeManager The Dumux TimeManager for simulation management.
     * \param gridView The grid view on the spatial domain of the problem
     */
    RichardsFringeProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name) + "_richards";
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

    /*!
     * \brief Returns the temperature [K] within a finite volume
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; }; // -> 10°C

    /*!
     * \brief Returns the reference pressure [Pa] of the non-wetting
     *        fluid phase within a finite volume
     *
     * This problem assumes a constant reference pressure of 1 bar.
     */
    Scalar nonWettingReferencePressure() const
    { return 1.0e5; };

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
        if (onUpperBoundary_(globalPos))
            bcTypes.setAllDirichlet();
        else
            bcTypes.setAllNeumann();
        static const bool bottomDirichlet = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Problem, BottomDirichlet);
        if (bottomDirichlet && onLowerBoundary_(globalPos))
            bcTypes.setAllDirichlet();
        return bcTypes;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param globalPos The position for which the Dirichlet value is set
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        if (onUpperBoundary_(globalPos))
        {
            static const Scalar sw = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, TopSaturation);
            const Scalar pc = MaterialLaw::pc(this->spatialParams().materialLawParamsAtPos(globalPos), sw);
            values[pressureIdx] = nonWettingReferencePressure() - pc;
            values.setState(Indices::bothPhases);
        }
        else if(onLowerBoundary_(globalPos))
        {
            static const Scalar sw = 1.0;
            const Scalar pc = MaterialLaw::pc(this->spatialParams().materialLawParamsAtPos(globalPos), sw);
            values[pressureIdx] = nonWettingReferencePressure() - pc;
            values.setState(Indices::bothPhases);
        }
        return values;
    }

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
        static const Scalar levelFactor = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, WaterLevelFactor);
        static const Scalar swBottom = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, SwBottomInitial);
        static const Scalar swTop = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, SwTopInitial);

        PrimaryVariables values(0.0);
        if (globalPos[1] > levelFactor*(this->bBoxMax()[1] - this->bBoxMin()[1]))
        {
            const Scalar sw = swTop;
            const Scalar pc = MaterialLaw::pc(this->spatialParams().materialLawParamsAtPos(globalPos), sw);
            values[pressureIdx] = nonWettingReferencePressure() - pc;
            values.setState(Indices::bothPhases);
        }
        else
        {
            const Scalar sw = swBottom;
            const Scalar pc = MaterialLaw::pc(this->spatialParams().materialLawParamsAtPos(globalPos), sw);
            values[pressureIdx] = nonWettingReferencePressure() - pc;
            values.setState(Indices::bothPhases);
        }
        return values;
    }

    bool shouldWriteRestartFile() const
    { return false; }

    // \}

private:

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->bBoxMax()[1] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->bBoxMin()[1] + eps_;
    }

    static constexpr Scalar eps_ = 1.5e-7;

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;
    std::string name_;
};

} //end namespace Dumux

#endif
