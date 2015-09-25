// $Id: lensproblem.hh 2539 2009-10-13 14:49:40Z bernd $
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
#ifndef DUMUX_LENSPROBLEM2P_HH
#define DUMUX_LENSPROBLEM2P_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>
#include "external_interface.hh"

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/simplednapl.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include <dumux/boxmodels/2p/2pmodel.hh>

#include "lensspatialparameters2p.hh"

namespace Dumux
{

template <class TypeTag>
class LensProblem;

//////////
// Specify the properties for the lens problem
//////////
namespace Properties
{
NEW_TYPE_TAG(LensProblem, INHERITS_FROM(BoxTwoP));

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

#if HAVE_DUNE_PDELAB
SET_PROP(LensProblem, LocalFEMSpace)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    enum{dim = GridView::dimension};

public:
    typedef Dune::PDELab::Q1LocalFiniteElementMap<Scalar,Scalar,dim>  type; // for cubes
//    typedef Dune::PDELab::P1LocalFiniteElementMap<Scalar,Scalar,dim>  type; // for simplices
};
#endif // HAVE_DUNE_PDELAB

// Set the problem property
SET_PROP(LensProblem, Problem)
{
    typedef Dumux::LensProblem<TypeTag> type;
};

// Set the wetting phase
SET_PROP(LensProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(LensProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleDNAPL<Scalar> > type;
};

// Set the spatial parameters
SET_PROP(LensProblem, SpatialParameters)
{
    typedef Dumux::LensSpatialParameters2p<TypeTag> type;
};

// Enable gravity
SET_BOOL_PROP(LensProblem, EnableGravity, false);
}

/*!
 * \ingroup TwoPBoxProblems
 * \brief Soil decontamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 *
 * The domain is sized 5m times 4m and features a rectangular lens
 * with low permeablility which spans from (1 m , 2 m) to (4 m, 3 m)
 * and is surrounded by a medium with higher permability.
 *
 * On the top and the bottom of the domain neumann boundary conditions
 * are used, while dirichlet conditions apply on the left and right
 * boundaries.
 *
 * DNAPL is injected at the top boundary from 3m to 4m at a rate of
 * 0.04 kg/(s m), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * The dirichlet boundaries on the left boundary is the hydrostatic
 * pressure scaled by a factor of 1.125, while on the right side it is
 * just the hydrostatic pressure. The DNAPL saturation on both sides
 * is zero.
 *
 * This problem uses the \ref TwoPBoxModel.
 *
 * This problem should typically simulated until \f$t_{\text{end}} =
 * 50\,000\;s\f$ is reached. A good choice for the initial time step size
 * is \f$t_{\text{inital}} = 1\,000\;s\f$.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./lens_2p 50000 100</tt>
 */
template <class TypeTag >
class LensProblem : public TwoPProblem<TypeTag>
{
    typedef LensProblem<TypeTag> ThisType;
    typedef TwoPProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef TwoPFluidState<TypeTag> FluidState;

    enum
    {
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

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

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
        bboxMin_[0] = lowerLeft[0];
        bboxMin_[1] = lowerLeft[1];
        bboxMax_[0] = upperRight[0];
        bboxMax_[1] = upperRight[1];

        //load interface-file
        Dumux::InterfaceProblemProperties interfaceProbProps("interface2p.xml");

        lowerPressure_ = interfaceProbProps.IPP_LowerPressure;
        upperPressure_ = interfaceProbProps.IPP_UpperPressure;
        infiltrationRate_ = interfaceProbProps.IPP_InfiltrationRate;
        //infiltrationStartTime_= interfaceProbProps.IPP_InfiltrationStartTime;
        infiltrationStartTime_= 1.0e-9;//The infiltrations starts always after the first time step!
        infiltrationEndTime_= interfaceProbProps.IPP_InfiltrationEndTime;
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
        std::string simName = "lens-2p_run";
        Dumux::InterfaceProblemProperties interfaceProbProps("interface2p.xml");
        Scalar simNum =  interfaceProbProps.IPP_SimulationNumber;

        return (str(boost::format("%s-%02d")
                    %simName%simNum).c_str());
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature(const Element &element,
                       const FVElementGeometry &fvElemGeom,
                       int scvIdx) const
    {
        return 273.15 + 10; // -> 10°C
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
    void boundaryTypes(BoundaryTypes &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();


        if (onUpperBoundary_(globalPos) || onLowerBoundary_(globalPos))
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        if (onInlet_(globalPos))
            values.setNeumann(contiNEqIdx);
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
    void dirichlet(PrimaryVariables &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        if (onUpperBoundary_(globalPos))
        {
            values[pwIdx] = upperPressure_;
            values[SnIdx] = 0.0;
        }
        else if (onLowerBoundary_(globalPos))
        {
            values[pwIdx] = lowerPressure_;
            values[SnIdx] = 0.0;
        }
        else
            values = 0.0;

//        Scalar densityW = this->wettingPhase().density();
//
//        if (onLeftBoundary_(globalPos))
//        {
//            Scalar depth = bboxMax_[1] - globalPos[1];
//
//            // hydrostatic pressure scaled by alpha
//            values[pW] = lowerPressure_ - densityW*this->gravity()[1]*depth;
//            values[sN] = 0.0;
//        }
//        else if (onRightBoundary_(globalPos))
//        {
//            Scalar depth = bboxMax_[1] - globalPos[1];
//
//            // hydrostatic pressure
//            values[pW] = lowerPressure_ - densityW*this->gravity()[1]*depth;
//            values[sN] = 0.0;
//        }
//        else
//            values = 0.0;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 const Intersection &isIt,
                 int scvIdx,
                 int boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);

        values = 0.0;

        const Scalar& time = this->timeManager().time();

        if (time >= infiltrationStartTime_ && time <= infiltrationEndTime_)
        {
            if (onInlet_(globalPos))
                values[contiNEqIdx] = -infiltrationRate_; // kg / (m * s)
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
    void source(PrimaryVariables &values,
                const Element &element,
                const FVElementGeometry &,
                int subControlVolumeIdx) const
    {
        values = Scalar(0.0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initial(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 int scvIdx) const
    {
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);

        // no DNAPL, hydrostatic pressure
        const Scalar depth = this->bboxMax()[1] - globalPos[1];
        const Scalar height = this->bboxMax()[1] - this->bboxMin()[1];

        values[pwIdx] = upperPressure_ - depth/height*(upperPressure_-lowerPressure_);
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
        return onUpperBoundary_(globalPos) && (bboxMax_[0]-0.35*width)/width > lambda && lambda > (bboxMax_[0]-0.55*width)/width;
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
