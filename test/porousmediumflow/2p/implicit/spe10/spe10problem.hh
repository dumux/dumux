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
#ifndef DUMUX_SPE10_PROBLEM_HH
#define DUMUX_SPE10_PROBLEM_HH

#include <dumux/material/components/simpleh2o.hh>
#include "spe10oil.hh"

#include <dumux/porousmediumflow/2p/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/implicit/cellcentered/propertydefaults.hh>

#include "spe10spatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class CCTwoPSpe10Problem;

namespace Properties
{
NEW_TYPE_TAG(TwoPSpe10Problem, INHERITS_FROM(TwoP, Spe10SpatialParams));
NEW_TYPE_TAG(CCTwoPSpe10Problem, INHERITS_FROM(CCModel, TwoPSpe10Problem));

// Set the grid type
SET_TYPE_PROP(TwoPSpe10Problem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(TwoPSpe10Problem, Problem, Dumux::CCTwoPSpe10Problem<TypeTag>);


SET_TYPE_PROP(TwoPSpe10Problem, LinearSolver, Dumux::SuperLUBackend<TypeTag> );

// Set the wetting phase
SET_PROP(TwoPSpe10Problem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(TwoPSpe10Problem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, Dumux::Oil<Scalar> > type;
};

// Enable gravity
SET_BOOL_PROP(TwoPSpe10Problem, ProblemEnableGravity, false);

}

template <class TypeTag >
class CCTwoPSpe10Problem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum
    {
        // primary variable indices
        pwIdx = Indices::pwIdx,
        snIdx = Indices::snIdx,

        // equation indices
        contiNEqIdx = Indices::contiNEqIdx,
        contiWEqIdx = Indices::contiWEqIdx,

        // phase indices
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        // world dimension
        dimWorld = GridView::dimensionworld,
        dim = GridView::dimension
    };


    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dimWorld,dimWorld> DimWorldMatrix;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    CCTwoPSpe10Problem(TimeManager &timeManager,
                const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        eps_ = 3e-6;
        temperature_ = 273.15 + 20; // -> 20°C

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
        episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, EpisodeLength);

        this->timeManager().startNextEpisode(episodeLength_);

        pIn_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, PressureIn);
        pOut_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, PressureOut);
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    {
        return name_;
    }

    /*!
     * \brief User defined output after the time integration
     *
     * Will be called diretly after the time integration.
     */
    void postTimeStep()
    {
        // Calculate storage terms
        PrimaryVariables storage;
        this->model().globalStorage(storage);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout<<"Storage: " << storage << std::endl;
        }
    }

    /*!
     * \brief Returns the temperature \f$ K \f$
     *
     * This problem assumes a uniform temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return temperature_; }

    /*!
     * \brief Returns the source term
     *
     * \param values Stores the source values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} / (m^\textrm{dim} \cdot s )] \f$
     * \param globalPos The global position
     */
    void sourceAtPos(PrimaryVariables &values,
                const GlobalPosition &globalPos) const
    {
        values = 0.0;
    }
    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param values Stores the value of the boundary type
     * \param globalPos The global position
     */
    void boundaryTypesAtPos(BoundaryTypes &bcTypes,
            const GlobalPosition &globalPos) const
    {
        if (globalPos[1] <  eps_ || globalPos[1] > this->bBoxMax()[1] - eps_)
        {
            bcTypes.setAllDirichlet();
        }
        else
        {
            bcTypes.setAllNeumann();
        }
    }


    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param values Stores the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param globalPos The global position
     */
    void dirichletAtPos(PrimaryVariables &values,
                        const GlobalPosition &globalPos) const
    {
        values = 0;

        if (globalPos[1] > this->bBoxMax()[1] - eps_)
        {
            values[pwIdx] = pOut_;
            values[snIdx] = 1.0;
        }
        else
        {
            values[pwIdx] = pIn_;
            values[snIdx] = 0.0;
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumannAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        values = 0.0;
    }
    // \}

    /*!
     * \name Volume terms
     */
    // \{


    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    void initialAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        typename GET_PROP_TYPE(TypeTag, FluidState) fluidState;
        fluidState.setTemperature(temperature_);
        fluidState.setPressure(FluidSystem::wPhaseIdx, 1e5);
        fluidState.setPressure(FluidSystem::nPhaseIdx, 1e5);

        Scalar densityW = FluidSystem::density(fluidState, FluidSystem::wPhaseIdx);

        Scalar depth = this->bBoxMax()[dim-1] - globalPos[dim-1];

        values[pwIdx] = pOut_ + (pIn_-pOut_)/this->bBoxMax()[dim-1] * depth - densityW*this->gravity()[dim-1]*depth;
        values[snIdx] = 1.0;
    }
    // \}

    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     */
    void addOutputVtkFields()
    {
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;

        unsigned numElements = this->gridView().size(0);

        //create required scalar fields
        ScalarField *permX = this->resultWriter().allocateManagedBuffer(numElements);
        ScalarField *permY = this->resultWriter().allocateManagedBuffer(numElements);
        ScalarField *eIdxGlobal = this->resultWriter().allocateManagedBuffer(numElements);
        ScalarField *volume = this->resultWriter().allocateManagedBuffer(numElements);

        for (const auto& element : elements(this->gridView()))
        {
            int eIdx = this->elementMapper().index(element);
            FVElementGeometry fvGeometry;
            fvGeometry.update(this->gridView(), element);

            const DimWorldMatrix K = this->spatialParams().intrinsicPermeability(element, fvGeometry, /*element data*/ 0);
            (*permX)[eIdx] = K[0][0];
            (*permY)[eIdx] = K[1][1];
            (*eIdxGlobal)[eIdx] = eIdx;
            (*volume)[eIdx] = element.geometry().volume();
        }

        //pass the scalar fields to the vtkwriter
        this->resultWriter().attachDofData(*permX, "PERMX", false); //element data
        this->resultWriter().attachDofData(*permY, "PERMY", false); //element data
        this->resultWriter().attachDofData(*eIdxGlobal, "eIdx", false); //element data
        this->resultWriter().attachDofData(*volume, "volume", false); //element data
    }

    bool shouldWriteOutput() const
    {
        return this->timeManager().timeStepIndex() == 0 ||
               this->timeManager().episodeWillBeFinished() ||
               this->timeManager().willBeFinished();
    }

    void episodeEnd()
    {
        this->timeManager().startNextEpisode(episodeLength_);
    }

    bool shouldWriteRestartFile() const
    {
        return false;
    }

private:
    Scalar temperature_;
    Scalar eps_;
    std::string name_;
    Scalar pIn_;
    Scalar pOut_;
    Scalar episodeLength_;
};
} //end namespace

#endif
