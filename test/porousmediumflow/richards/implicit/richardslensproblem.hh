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

#include <dune/grid/io/file/dgfparser.hh>

#include <dumux/porousmediumflow/richards/implicit/model.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/linear/amgbackend.hh>

#include "richardslensspatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class RichardsLensProblem;


// Specify the properties for the lens problem
namespace Properties
{
NEW_TYPE_TAG(RichardsLensProblem, INHERITS_FROM(Richards, RichardsLensSpatialParams));
NEW_TYPE_TAG(RichardsLensBoxProblem, INHERITS_FROM(BoxModel, RichardsLensProblem));
NEW_TYPE_TAG(RichardsLensCCProblem, INHERITS_FROM(CCModel, RichardsLensProblem));

// Use 2d YaspGrid
SET_TYPE_PROP(RichardsLensProblem, Grid, Dune::YaspGrid<2>);

// Set the physical problem to be solved
SET_TYPE_PROP(RichardsLensProblem, Problem, Dumux::RichardsLensProblem<TypeTag>);

// Set the wetting phase
SET_PROP(RichardsLensProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Enable gravity
SET_BOOL_PROP(RichardsLensProblem, ProblemEnableGravity, true);

// Use the AMG backend to allow parallel computation
SET_TYPE_PROP(RichardsLensProblem, LinearSolver, AMGBackend<TypeTag>);

//! Use pressure [Pa] or pressure head [cm] formulation
SET_BOOL_PROP(RichardsLensProblem, UseHead, false);
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
class RichardsLensProblem : public RichardsProblem<TypeTag>
{
    typedef RichardsProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // copy some indices for convenience
        pwIdx = Indices::pwIdx,
        hIdx = Indices::hIdx,
        contiEqIdx = Indices::contiEqIdx,

        // Grid and world dimension
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    static const bool useHead = GET_PROP_VALUE(TypeTag, UseHead);

public:
    /*!
     * \brief Constructor
     *
     * \param timeManager The Dumux TimeManager for simulation management.
     * \param gridView The grid view on the spatial domain of the problem
     */
    RichardsLensProblem(TimeManager &timeManager,
                        const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        eps_ = 3e-6;
        pnRef_ = 1e5;

        lensLowerLeft_[0] = 1.0;
        lensLowerLeft_[1] = 2.0;

        lensUpperRight_[0] = 4.0;
        lensUpperRight_[1] = 3.0;

        this->spatialParams().setLensCoords(lensLowerLeft_, lensUpperRight_);

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
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
    const std::string name() const
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
     *
     * \param element The DUNE Codim<0> entity which intersects with
     *                the finite volume in question
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The sub control volume index inside the finite
     *               volume geometry
     */
    Scalar referencePressure(const Element &element,
                             const FVElementGeometry &fvGeometry,
                             const int scvIdx) const
    { return pnRef_; };

    /*!
     * \brief Return the sources within the domain.
     *
     * \param values Stores the source values, acts as return value
     * \param globalPos The global position
     */
    void sourceAtPos(PrimaryVariables &values,
                const GlobalPosition &globalPos) const
    {
        values = 0;
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
     * \param values The boundary types for the conservation equations
     * \param globalPos The position for which the boundary type is set
     */
    void boundaryTypesAtPos(BoundaryTypes &values,
                       const GlobalPosition &globalPos) const
    {
        if (onLeftBoundary_(globalPos) ||
            onRightBoundary_(globalPos))
        {
            values.setAllDirichlet();
        }
        else
            values.setAllNeumann();
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The position for which the Dirichlet value is set
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichletAtPos(PrimaryVariables &values,
                   const GlobalPosition &globalPos) const
    {
        // use initial values as boundary conditions
        if (onInlet_(globalPos))
            // inflow of water
            values[contiEqIdx] = -0.;
        else
            initial_(values, globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     *
     * \param values The neumann values for the conservation equations
     * \param globalPos The position for which the Neumann value is set
     */
    void neumannAtPos(PrimaryVariables &values,
                 const GlobalPosition &globalPos) const
    {
        values = 0.0;
        if (onInlet_(globalPos)) {
            // inflow of water
            if (useHead)
                values[contiEqIdx] = -0.04/1000.; // kg/(m*s) / density
            else
                values[contiEqIdx] = -0.04; // kg/(m*s)
        }
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
     * \param values Storage for all primary variables of the initial condition
     * \param globalPos The position for which the boundary type is set
     */
    void initialAtPos(PrimaryVariables &values,
                 const GlobalPosition &globalPos) const
    { initial_(values, globalPos); };

    // \}

private:
    void initial_(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        Scalar sw = 0.0;
        Scalar pc =
            MaterialLaw::pc(this->spatialParams().materialLawParams(globalPos),
                            sw);
        Scalar g = this->gravity().two_norm();

        if (useHead)
            values[hIdx] = (-pc)/1000./ g *(100.);
        else
            values[pwIdx] = pnRef_ - pc;
    }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->bBoxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->bBoxMax()[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->bBoxMin()[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->bBoxMax()[1] - eps_;
    }

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        Scalar width = this->bBoxMax()[0] - this->bBoxMin()[0];
        Scalar lambda = (this->bBoxMax()[0] - globalPos[0])/width;
        return onUpperBoundary_(globalPos) && 0.5 < lambda && lambda < 2.0/3.0;
    }

    Scalar eps_;
    Scalar pnRef_;

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;
    std::string name_;
};
} //end namespace

#endif
