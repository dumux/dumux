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
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 */

#ifndef DUMUX_GENERALLENSPROBLEM_HH
#define DUMUX_GENERALLENSPROBLEM_HH

//common includes
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/dnapl.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

//box model
#include <dumux/implicit/common/implicitporousmediaproblem.hh>
#include <dumux/implicit/2p/2pmodel.hh>

//decoupled model
#include <dumux/decoupled/2p/diffusion/fv/fvpressureproperties2p.hh>
#include <dumux/decoupled/2p/transport/fv/fvtransportproperties2p.hh>
#include <dumux/decoupled/2p/impes/impesproblem2p.hh>

#include "generallensspatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class GeneralLensProblem;

//////////
// Specify the properties for the lens problem
//////////
namespace Properties
{

//Set the general problem TypeTag which does not depend on the model
NEW_TYPE_TAG(GeneralLensProblem, INHERITS_FROM(GeneralLensSpatialParams));

// Property for defining the model specific problem base class
NEW_PROP_TAG(ProblemBaseClass);

// Set the grid type
SET_TYPE_PROP(GeneralLensProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(GeneralLensProblem, Problem, Dumux::GeneralLensProblem<TypeTag>);

// Set the wetting phase
SET_PROP(GeneralLensProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(GeneralLensProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::DNAPL<Scalar> > type;
};

///////////////////////////////////////////////////
// Box model TypeTag
//////////////////////////////////////////////////

NEW_TYPE_TAG(BoxGeneralLensProblem, INHERITS_FROM(BoxTwoP, GeneralLensProblem));

// Set the problem property
SET_TYPE_PROP(BoxGeneralLensProblem, ProblemBaseClass, Dumux::ImplicitPorousMediaProblem<TypeTag>);

// Set the problem property
SET_TYPE_PROP(BoxGeneralLensProblem, SpatialParamsBaseClass,Dumux::ImplicitSpatialParams<TypeTag>);


///////////////////////////////////////////////////
// CC model TypeTag
//////////////////////////////////////////////////

NEW_TYPE_TAG(CCGeneralLensProblem, INHERITS_FROM(CCTwoP, GeneralLensProblem));

// Set the problem property
SET_TYPE_PROP(CCGeneralLensProblem, ProblemBaseClass, Dumux::ImplicitPorousMediaProblem<TypeTag>);

// Set the problem property
SET_TYPE_PROP(CCGeneralLensProblem, SpatialParamsBaseClass,Dumux::ImplicitSpatialParams<TypeTag>);


///////////////////////////////////////////////////
// Deoupled model TypeTag
//////////////////////////////////////////////////

NEW_TYPE_TAG(DecoupledGeneralLensProblem, INHERITS_FROM(FVPressureTwoP, FVTransportTwoP, IMPESTwoP, GeneralLensProblem));

// Set the problem property
SET_TYPE_PROP(DecoupledGeneralLensProblem, ProblemBaseClass, Dumux::IMPESProblem2P<TypeTag>);

// Set the problem property
SET_TYPE_PROP(DecoupledGeneralLensProblem, SpatialParamsBaseClass, Dumux::FVSpatialParams<TypeTag>);

SET_INT_PROP(DecoupledGeneralLensProblem, Formulation,
        DecoupledTwoPCommonIndices::pwsn);
}

/*!
 * \ingroup TwoPModel
 * \ingroup IMPETtests
 *
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 *
 * The domain is sized 6m times 4m and features a rectangular lens
 * with low permeablility which spans from (1 m , 2 m) to (4 m, 3 m)
 * and is surrounded by a medium with higher permability. Note that
 * this problem is discretized using only two dimensions, so from the
 * point of view of the two-phase model, the depth of the domain
 * implicitly is 1 m everywhere.
 *
 * On the top and the bottom of the domain neumann boundary conditions
 * are used, while dirichlet conditions apply on the left and right
 * boundaries.
 *
 * DNAPL is injected at the top boundary from 3m to 4m at a rate of
 * 0.04 kg/(s m^2), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * The dirichlet boundaries on the left boundary is the hydrostatic
 * pressure scaled by a factor of 1.125, while on the right side it is
 * just the hydrostatic pressure. The DNAPL saturation on both sides
 * is zero.
 *
 * This problem uses the \ref TwoPModel.
 *
 * This problem should typically be simulated until \f$t_{\text{end}}
 * \approx 20\,000\;s\f$ is reached. A good choice for the initial time step
 * size is \f$t_{\text{inital}} = 250\;s\f$.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_generalproblem2p [-ModelType Box/Decoupled]</tt>
 */
template <class TypeTag >
class GeneralLensProblem : public GET_PROP_TYPE(TypeTag, ProblemBaseClass)
{
    typedef typename GET_PROP_TYPE(TypeTag, ProblemBaseClass) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;

    enum {
        // primary variable indices
        pwIdx = Indices::pwIdx,
        snIdx = Indices::snIdx,

        // phase indices
        nPhaseIdx = Indices::nPhaseIdx,

        // world dimension
        dimWorld = GridView::dimensionworld
    };


    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    GeneralLensProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        eps_ = 3e-6;
        temperature_ = 273.15 + 20; // -> 20°C
        this->timeManager().startNextEpisode(500);
    }

    /*!
     * \name Problem parameters
     */
    // \{


    bool shouldWriteRestartFile()
    {
        return false;
    }

    bool shouldWriteOutput() const
    {
        if (this->timeManager().time() < eps_ ||
            this->timeManager().willBeFinished() ||
            this->timeManager().episodeWillBeOver())
        {
            return true;
        }
        return false;
    }

    void episodeEnd()
    {
        if (!this->timeManager().willBeFinished())
        this->timeManager().startNextEpisode(500);
    }

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string name() const
    {
        return "generallens_" + GET_RUNTIME_PARAM(TypeTag, std::string, ModelType);
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * \param globalPos The global coordinates
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    {
        return temperature_;
    }

    //! Returns the reference pressure for evaluation of constitutive relations
    /*!
    * \param globalPos The global coordinates
    *
    * Only for decoupled model -> incompressible.
    */
    Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
    {
        return 1.0e5; // -> 1 bar
    }

    /*! \brief Returns the source term
     * \param values return values
     * \param globalPos The global coordinates
     */
    void sourceAtPos(PrimaryVariables &values,const GlobalPosition& globalPos) const
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
     *        used for which equation on a given boundary control volume.
     *
     * \param values The boundary types for the conservation equations
    * \param globalPos The global coordinates of the boundary
     */
    void boundaryTypesAtPos(BoundaryTypes &values,
                            const GlobalPosition &globalPos) const
    {
        if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos)) {
            values.setAllDirichlet();
        }
        else {
            values.setAllNeumann();
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary.
     *
     * \param values The dirichlet values for the primary variables
    * \param globalPos The global coordinates of the boundary
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichletAtPos(PrimaryVariables &values,
                        const GlobalPosition &globalPos) const
    {
        Scalar densityW = WettingPhase::density(temperature_, /*pressure=*/1e5);

        if (onLeftBoundary_(globalPos))
        {
            Scalar height = this->bBoxMax()[1] - this->bBoxMin()[1];
            Scalar depth = this->bBoxMax()[1] - globalPos[1];
            Scalar alpha = (1 + 1.5/height);

            // hydrostatic pressure scaled by alpha
            values[pwIdx] = 1e5 - alpha*densityW*this->gravity()[1]*depth;
            values[snIdx] = 0.0;
        }
        else if (onRightBoundary_(globalPos))
        {
            Scalar depth = this->bBoxMax()[1] - globalPos[1];

            // hydrostatic pressure
            values[pwIdx] = 1e5 - densityW*this->gravity()[1]*depth;
            values[snIdx] = 0.0;
        }
        else
            values = 0.0;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values The neumann values for the conservation equations [kg / (m^2 *s )]
     * \param globalPos The position of the boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumannAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        values = 0.0;
        if (onInlet_(globalPos)) {
            values[nPhaseIdx] = -0.04; // kg / (m * s)
        }
    }
    // \}

    /*!
     * \name Volume terms
     */
    // \{


    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initialAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        Scalar depth = this->bBoxMax()[1] - globalPos[1];
        Scalar densityW = WettingPhase::density(temperature_, /*pressure=*/1.0e5);

        // hydrostatic pressure
        values[pwIdx] = 1e5 - densityW*this->gravity()[1]*depth;
        values[snIdx] = 0.0;

    }
    // \}

private:

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

    Scalar temperature_;
    Scalar eps_;
};
} //end namespace

#endif
