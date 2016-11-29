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
 * \brief Definition of a problem for water management in PEM fuel cells.
 */
#ifndef DUMUX_FUELCELL_PROBLEM_HH
#define DUMUX_FUELCELL_PROBLEM_HH

#include <dumux/porousmediumflow/2pnc/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/material/fluidsystems/h2on2o2.hh>
#include <dumux/material/constants.hh>
#include <dumux/material/chemistry/electrochemistry/electrochemistry.hh>

#include "fuelcellspatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class FuelCellProblem;

namespace Properties
{
NEW_TYPE_TAG(FuelCellProblem, INHERITS_FROM(TwoPNC, FuelCellSpatialParams));
NEW_TYPE_TAG(FuelCellBoxProblem, INHERITS_FROM(BoxModel, FuelCellProblem));
NEW_TYPE_TAG(FuelCellCCProblem, INHERITS_FROM(CCModel, FuelCellProblem));

// Set the grid type
SET_TYPE_PROP(FuelCellProblem, Grid, Dune::YaspGrid<2>);
// Set the problem property
SET_TYPE_PROP(FuelCellProblem, Problem, FuelCellProblem<TypeTag>);
// Set the primary variable combination for the 2pnc model
SET_INT_PROP(FuelCellProblem, Formulation, TwoPNCFormulation::pnsw);

// Set fluid configuration
SET_PROP(FuelCellProblem, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    static const bool useComplexRelations = true;
 public:
    typedef FluidSystems::H2ON2O2<Scalar, useComplexRelations> type;
};

// Set the transport equation that is replaced by the total mass balance
SET_INT_PROP(FuelCellProblem, ReplaceCompEqIdx, 3);
}


/*!
 * \ingroup TwoPNCModel
 * \ingroup ImplicitTestProblems
 * \brief Problem or water management in PEM fuel cells.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2pnc</tt>
 */
template <class TypeTag>
class FuelCellProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    enum
    {
        numComponents = FluidSystem::numComponents,
        numSecComponents = FluidSystem::numSecComponents,
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };
    enum
    {
        wCompIdx = FluidSystem::wCompIdx, //major component of the liquid phase
        nCompIdx = FluidSystem::nCompIdx, //major component of the gas phase
    };
    enum
    {
        pressureIdx = Indices::pressureIdx, //gas-phase pressure
        switchIdx = Indices::switchIdx, //liquid saturation or mole fraction
        conti0EqIdx = Indices::conti0EqIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;

    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    // Select the electrochemistry method
    typedef Dumux::ElectroChemistry<TypeTag, ElectroChemistryModel::Ochs> ElectroChemistry;
    typedef Constants<Scalar> Constant;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };
public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    FuelCellProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        nTemperature_       = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, NTemperature);
        nPressure_          = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, NPressure);
        pressureLow_        = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, PressureLow);
        pressureHigh_       = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, PressureHigh);
        temperatureLow_     = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, TemperatureLow);
        temperatureHigh_    = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, TemperatureHigh);
        temperature_        = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FluidSystem, InitialTemperature);

        name_               = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);

        pO2Inlet_            = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, ElectroChemistry, pO2Inlet);

        eps_ = 1e-6;

        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);
    }

    /*!
     * \name Problem parameters
     */

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string name() const
    { return name_; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return temperature_; }

    //! \copydoc ImplicitProblem::solDependentSource()
    void solDependentSource(PrimaryVariables &values,
                            const Element &element,
                            const FVElementGeometry &fvGeometry,
                            const int scvIdx,
                            const ElementVolumeVariables &elemVolVars) const
    {
        values = 0.0;

        const auto& globalPos = isBox ? element.geometry().corner(scvIdx)
                                      : element.geometry().center();

        //reaction sources from electro chemistry
        if(inReactionLayer_(globalPos))
        {
            const auto& volVars = elemVolVars[scvIdx];
            auto currentDensity = ElectroChemistry::calculateCurrentDensity(volVars);
            ElectroChemistry::reactionSource(values, currentDensity);
        }
    }


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
    void boundaryTypesAtPos(BoundaryTypes &values,
                            const GlobalPosition &globalPos) const
    {
        values.setAllNeumann();

        if (onUpperBoundary_(globalPos)){
            values.setAllDirichlet();
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
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);

        if(onUpperBoundary_(globalPos))
        {
            Scalar pg = 1.0e5;
            values[pressureIdx] = pg;
            values[switchIdx] = 0.3;//Sl for bothPhases
            values[switchIdx+1] = pO2Inlet_/4.315e9; //moleFraction xlO2 for bothPhases
        }
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann
     *        boundary segment in dependency on the current solution.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param intersection The intersection between element and boundary
     * \param scvIdx The local index of the sub-control volume
     * \param boundaryFaceIdx The index of the boundary face
     * \param elemVolVars All volume variables for the element
     *
     * This method is used for cases, when the Neumann condition depends on the
     * solution and requires some quantities that are specific to the fully-implicit method.
     * The \a values store the mass flux of each phase normal to the boundary.
     * Negative values indicate an inflow.
     */
     void solDependentNeumann(PrimaryVariables &values,
                      const Element &element,
                      const FVElementGeometry &fvGeometry,
                      const Intersection &intersection,
                      const int scvIdx,
                      const int boundaryFaceIdx,
                      const ElementVolumeVariables &elemVolVars) const
    { values = 0;}

    /*!
     * \name Volume terms
     */


    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    void initialAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }

    /*!
     * \brief Return the initial phase state inside a sub control volume.
     *
     * \param element The element of the sub control volume
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The sub control volume index
     */

    int initialPhasePresence(const Element &element,
                             const FVElementGeometry &fvGeometry,
                             int scvIdx) const
    {
        return Indices::bothPhases;
    }

    /*!
     * \brief Add problem specific vtk output for the electrochemistry
     */
    void addOutputVtkFields()
    {
        // add the output field specific to the electrochemistry
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // get the number of degrees of freedom
        auto numDofs = this->model().numDofs();

        // create the required scalar fields
        ScalarField *currentDensity = this->resultWriter().allocateManagedBuffer (numDofs);
        ScalarField *reactionSourceH2O = this->resultWriter().allocateManagedBuffer (numDofs);
        ScalarField *reactionSourceO2 = this->resultWriter().allocateManagedBuffer (numDofs);

        for (const auto& element : elements(this->gridView()))
        {
            FVElementGeometry fvGeometry;
            fvGeometry.update(this->gridView(), element);

            ElementVolumeVariables elemVolVars;
            elemVolVars.update(*this,
                               element,
                               fvGeometry,
                               false /* oldSol? */);

            for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
            {
                const auto& globalPos = isBox ? element.geometry().corner(scvIdx)
                                              : element.geometry().center();

                auto dofIdxGlobal = this->model().dofMapper().subIndex(element, scvIdx, dofCodim);

                //reaction sources from electro chemistry
                if(inReactionLayer_(globalPos))
                {
                    //reactionSource Output
                    PrimaryVariables source;
                    auto i = ElectroChemistry::calculateCurrentDensity(elemVolVars[scvIdx]);
                    ElectroChemistry::reactionSource(source, i);

                    (*reactionSourceH2O)[dofIdxGlobal] = source[wPhaseIdx];
                    (*reactionSourceO2)[dofIdxGlobal] = source[numComponents-1];

                    //Current Output in A/cm^2
                    (*currentDensity)[dofIdxGlobal] = i/10000;
                }
                else
                {
                    (*reactionSourceH2O)[dofIdxGlobal] = 0.0;
                    (*reactionSourceO2)[dofIdxGlobal] = 0.0;
                    (*currentDensity)[dofIdxGlobal] = 0.0;
                }
            }
        }

        this->resultWriter().attachDofData(*reactionSourceH2O, "reactionSourceH2O [mol/(sm^2)]", isBox);
        this->resultWriter().attachDofData(*reactionSourceO2, "reactionSourceO2 [mol/(sm^2)]", isBox);
        this->resultWriter().attachDofData(*currentDensity, "currentDensity [A/cm^2]", isBox);
    }

private:

    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        Scalar pg = 1.0e5;
        values[pressureIdx] = pg;
        values[switchIdx] = 0.3;//Sl for bothPhases
        values[switchIdx+1] = pO2Inlet_/4.315e9; //moleFraction xlO2 for bothPhases
    }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->bBoxMax()[1] - eps_; }

    bool inReactionLayer_(const GlobalPosition& globalPos) const
    { return globalPos[1] < 0.1*(this->bBoxMax()[1] - this->bBoxMin()[1]) + eps_; }

    Scalar temperature_;
    Scalar eps_;
    int nTemperature_;
    int nPressure_;
    std::string name_ ;
    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;
    Scalar pO2Inlet_;
};

} //end namespace

#endif
