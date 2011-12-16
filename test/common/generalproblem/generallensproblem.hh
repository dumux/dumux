/*****************************************************************************
 *   Copyright (C) 2011 by Markus Wolff                                      *
 *   Copyright (C) 2007-2008 by Klaus Mosthaf                                *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
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
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/simplednapl.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

//box model
#include <dumux/boxmodels/2p/2pmodel.hh>

//decoupled model
#include <dumux/decoupled/2p/impes/impesproblem2p.hh>
#include <dumux/decoupled/2p/diffusion/fv/fvvelocity2p.hh>
#include <dumux/decoupled/2p/transport/fv/fvsaturation2p.hh>
#include<dumux/decoupled/2p/transport/fv/evalcflflux_coats.hh>

#include "generallensspatialparameters.hh"

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
NEW_TYPE_TAG(GeneralLensProblem, INHERITS_FROM(GeneralLensSpatialParameters));

// Property for defining the model specific problem base class
NEW_PROP_TAG(ProblemBaseClass);

// Set the grid type
SET_PROP(GeneralLensProblem, Grid)
{
//    typedef Dune::UGGrid<2> type;
    typedef Dune::YaspGrid<2> type;
};

// Set the problem property
SET_PROP(GeneralLensProblem, Problem)
{
    typedef Dumux::GeneralLensProblem<TypeTag> type;
};

// Set the wetting phase
SET_PROP(GeneralLensProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(GeneralLensProblem, NonWettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleDNAPL<Scalar> > type;
};

// Enable gravity
SET_BOOL_PROP(GeneralLensProblem, EnableGravity, true);

///////////////////////////////////////////////////
// Box model TypeTag
//////////////////////////////////////////////////

NEW_TYPE_TAG(BoxGeneralLensProblem, INHERITS_FROM(GeneralLensProblem, BoxTwoP));


// Set the problem property
SET_PROP(BoxGeneralLensProblem, ProblemBaseClass)
{
    typedef Dumux::TwoPProblem<TypeTag> type;
};

// Set the problem property
SET_PROP(BoxGeneralLensProblem, SpatialParamsBaseClass)
{
    typedef Dumux::BoxSpatialParameters<TypeTag> type;
};


///////////////////////////////////////////////////
// Deoupled model TypeTag
//////////////////////////////////////////////////

NEW_TYPE_TAG(DecoupledGeneralLensProblem, INHERITS_FROM(DecoupledTwoP, Transport, GeneralLensProblem));


// Set the problem property
SET_PROP(DecoupledGeneralLensProblem, ProblemBaseClass)
{
    typedef Dumux::IMPESProblem2P<TypeTag> type;
};

// Set the problem property
SET_PROP(DecoupledGeneralLensProblem, SpatialParamsBaseClass)
{
    typedef Dumux::FVSpatialParameters<TypeTag> type;
};

// Set the model properties
SET_TYPE_PROP(DecoupledGeneralLensProblem, TransportModel, Dumux::FVSaturation2P<TTAG(DecoupledGeneralLensProblem)>);

SET_PROP(DecoupledGeneralLensProblem, PressureModel)
{
    typedef Dumux::FVVelocity2P<TTAG(DecoupledGeneralLensProblem)> type;
};

SET_INT_PROP(DecoupledGeneralLensProblem, Formulation,
        DecoupledTwoPCommonIndices::pwSn);

SET_TYPE_PROP(DecoupledGeneralLensProblem, EvalCflFluxFunction, Dumux::EvalCflFluxCoats<TypeTag>);

SET_SCALAR_PROP(DecoupledGeneralLensProblem, CFLFactor, 1.0);
}

/*!
 * \ingroup TwoPBoxModel
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
 * This problem uses the \ref TwoPBoxModel.
 *
 * This problem should typically be simulated until \f$t_{\text{end}}
 * \approx 20\,000\;s\f$ is reached. A good choice for the initial time step
 * size is \f$t_{\text{inital}} = 250\;s\f$.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_2p 20000 250</tt>
 */
template <class TypeTag >
class GeneralLensProblem : public GET_PROP_TYPE(TypeTag, PTAG(ProblemBaseClass))
{
    typedef GeneralLensProblem<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ProblemBaseClass)) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(WettingPhase)) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NonWettingPhase)) NonWettingPhase;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),

        // primary variable indices
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        pwIdx = Indices::pwIdx,
        SnIdx = Indices::SnIdx,

        // equation indices
        contiWEqIdx = Indices::contiWEqIdx,
        contiNEqIdx = Indices::contiNEqIdx,

        // phase indices
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,


        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };


    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     * \param lensLowerLeft Global position of the lenses lower left corner
     * \param lensUpperRight Global position of the lenses upper right corner
     */
    GeneralLensProblem(TimeManager &timeManager,
                       const GridView &gridView,
                       const GlobalPosition &lensLowerLeft,
                       const GlobalPosition &lensUpperRight)
        : ParentType(timeManager, gridView)
    {
        temperature_ = 273.15 + 20; // -> 20°C
        this->spatialParameters().setLensCoords(lensLowerLeft, lensUpperRight);
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
    };

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
    };

    //! Returns the reference pressure for evaluation of constitutive relations
    /*!
    * \param globalPos The global coordinates
    *
    * Only for decoupled model -> incompressible.
    */
    Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
    {
        return 1e5; // -> 10°C
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
            Scalar height = this->bboxMax()[1] - this->bboxMin()[1];
            Scalar depth = this->bboxMax()[1] - globalPos[1];
            Scalar alpha = (1 + 1.5/height);

            // hydrostatic pressure scaled by alpha
            values[pwIdx] = 1e5 - alpha*densityW*this->gravity()[1]*depth;
            values[SnIdx] = 0.0;
        }
        else if (onRightBoundary_(globalPos))
        {
            Scalar depth = this->bboxMax()[1] - globalPos[1];

            // hydrostatic pressure
            values[pwIdx] = 1e5 - densityW*this->gravity()[1]*depth;
            values[SnIdx] = 0.0;
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
        Scalar depth = this->bboxMax()[1] - globalPos[1];
        Scalar densityW = WettingPhase::density(temperature_, /*pressure=*/1e5);

        // hydrostatic pressure
        values[pwIdx] = 1e5 - densityW*this->gravity()[1]*depth;
        values[SnIdx] = 0.0;

    }
    // \}

private:

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->bboxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->bboxMax()[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->bboxMin()[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->bboxMax()[1] - eps_;
    }

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        Scalar width = this->bboxMax()[0] - this->bboxMin()[0];
        Scalar lambda = (this->bboxMax()[0] - globalPos[0])/width;
        return onUpperBoundary_(globalPos) && 0.5 < lambda && lambda < 2.0/3.0;
    }

    Scalar temperature_;
    static constexpr Scalar eps_ = 3e-6;
};
} //end namespace

#endif
