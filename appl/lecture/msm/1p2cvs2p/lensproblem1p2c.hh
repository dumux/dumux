/*****************************************************************************
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
#ifndef DUMUX_LENSPROBLEM1P2C_HH
#define DUMUX_LENSPROBLEM1P2C_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/simplednapl.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include <dumux/boxmodels/1p2c/1p2cmodel.hh>

#include <water_contaminant.hh>
#include "lensspatialparameters1p2c.hh"

namespace Dumux
{

template <class TypeTag>
class LensProblem;

//////////
// Specify the properties for the lens problem
//////////
namespace Properties
{
NEW_TYPE_TAG(LensProblem, INHERITS_FROM(BoxOnePTwoC));

// Set the grid type
SET_PROP(LensProblem, Grid)
{
#if HAVE_UG
    typedef Dune::UGGrid<2> type;
#else
    typedef Dune::YaspGrid<2> type;
    //typedef Dune::SGrid<2, 2> type;
#endif
};

// Set the problem property
SET_PROP(LensProblem, Problem)
{
    typedef Dumux::LensProblem<TypeTag> type;
};

// Set fluid configuration
SET_PROP(LensProblem, FluidSystem)
{
    typedef Dumux::WaterContaminant<TypeTag> type;
};

// Set the spatial parameters
SET_PROP(LensProblem, SpatialParameters)
{
    typedef Dumux::LensSpatialParameters1p2c<TypeTag> type;
};

// Define whether mole(true) or mass(false) fractions are used
SET_BOOL_PROP(LensProblem, UseMoles, true);
// Disable gravity
SET_BOOL_PROP(LensProblem, EnableGravity, false);

// Disable Jacobian recycling
SET_BOOL_PROP(LensProblem, EnableJacobianRecycling, false);
}

/*!
 * \ingroup TwoPBoxProblems
 * \brief Soil decontamination problem where a contaminant infiltrates a fully
 *        water saturated medium.
 *
 * The whole problem is symmetric.
 *
 * The domain is sized 5m times 4m and features a rectangular lens
 * with low permeablility which spans from (0.8 m , 2 m) to (4 m, 3 m)
 * and is surrounded by a medium with higher permability.
 *
 * On the top and the bottom of the domain Dirichlet boundary conditions
 * are used, while Neumann conditions apply on the left and right
 * boundaries.
 *
 * The contaminant is injected at the top boundary from 2.25m to 2.75m at a variable rate.
 *
 * The Dirichlet boundaries on the top boundary is different to the bottom pressure.
 *
 * This problem uses the \ref TwoPBoxModel.
 *
 * This problem should typically simulated until \f$t_{\text{end}} =
 * 30\,000\;s\f$ is reached. A good choice for the initial time step size
 * is \f$t_{\text{inital}} = 100\;s\f$.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./lens_1p2c 30000 100</tt>
 */
template <class TypeTag >
class LensProblem : public OnePTwoCBoxProblem<TypeTag>
{
    typedef LensProblem<TypeTag> ThisType;
    typedef OnePTwoCBoxProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePTwoCIndices)) Indices;
    enum
    {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        // indices of the equations
        contiEqIdx = Indices::contiEqIdx,
        transEqIdx = Indices::transEqIdx,
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        x1Idx = Indices::x1Idx
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef typename GET_PROP(TypeTag, PTAG(ParameterTree)) Params;

public:
    LensProblem(TimeManager &timeManager,
                const GridView &gridView,
                const GlobalPosition &lowerLeft,
                const GlobalPosition &upperRight,
                const GlobalPosition &lensLowerLeft,
                const GlobalPosition &lensUpperRight)
        : ParentType(timeManager, gridView)
    {
        this->spatialParameters().setLensCoords(lensLowerLeft, lensUpperRight);

        bboxMin_ = lowerLeft;
        bboxMax_ = upperRight;

        upperPressure_ = Params::tree().template get<double>("Boundary.UpperPressure");
        lowerPressure_ = Params::tree().template get<double>("Boundary.LowerPressure");
        infiltrationRate_ = Params::tree().template get<double>("Boundary.InfiltrationRate");
        infiltrationStartTime_= 1.0e-9;//The infiltrations starts always after the first time step!
        infiltrationEndTime_= Params::tree().template get<double>("Boundary.InfiltrationEndTime");
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const char *name() const
    {
        std::string simName = "lens-1p2c_run";
        return (str(boost::format("%s")
                    %simName).c_str());
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    {
        return 273.15 + 10; // -> 10��C
    };

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
     * \param vertex The vertex on the boundary for which the
     *               conditions needs to be specified
     */
    void boundaryTypesAtPos(BoundaryTypes &values, const GlobalPosition& globalPos) const
    {
        if (onUpperBoundary_(globalPos) || onLowerBoundary_(globalPos))
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        if (onInlet_(globalPos))
            values.setNeumann(transEqIdx);
//        if (onLowerBoundary_(globalPos))
//            values.setNeumann(transEqIdx);
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param vertex The vertex representing the "half volume on the boundary"
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
//    	if (onInlet_(globalPos) && this->timeManager().time() <= infiltrationEndTime_)
//    	{
//    		values[pressureIdx] = upperPressure_;
//    		values[x1Idx] =0.8;
//    	}
//    	else
    		if (onUpperBoundary_(globalPos))
        {
            values[pressureIdx] = upperPressure_;
            values[x1Idx] = 0.0;
        }
        else if (onLowerBoundary_(globalPos))
        {
            values[pressureIdx] = lowerPressure_;
            values[x1Idx] = 0.0;
        }
        else
            values = 0.0;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
        values = 0.0;
        const Scalar& time = this->timeManager().time();
        if (time <= infiltrationEndTime_ && infiltrationStartTime_ <= time)
        {
            if (onInlet_(globalPos))
                values[transEqIdx] = -infiltrationRate_;
        }
    }
    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void sourceAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
        values = Scalar(0.0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initialAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
        // no contaminant, hydrostatic pressure
        const Scalar depth = this->bboxMax()[1] - globalPos[1];
        const Scalar height = this->bboxMax()[1] - this->bboxMin()[1];

        values[contiEqIdx] = upperPressure_ - depth/height*(upperPressure_-lowerPressure_);
        //if (globalPos[1]>3.5 && globalPos[1]<3.75 && globalPos[0]>2.25 && globalPos[0]<2.75 )
        //{ values[transEqIdx] = 0.001;}
        //else
        values[transEqIdx] = 0.0;
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
        return onUpperBoundary_(globalPos)
            && (bboxMax_[0] - 0.45*width)/width > lambda
            && lambda > (bboxMax_[0] - 0.55*width)/width;
    }

    static const Scalar eps_ = 3e-6;
    GlobalPosition bboxMin_;
    GlobalPosition bboxMax_;

    Scalar upperPressure_;
    Scalar lowerPressure_;
    Scalar infiltrationRate_;
    Scalar infiltrationStartTime_;
    Scalar infiltrationEndTime_;
};
} //end namespace

#endif
