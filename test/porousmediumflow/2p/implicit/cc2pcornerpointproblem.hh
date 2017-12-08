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
#ifndef DUMUX_CC2P_CORNERPOINT_PROBLEM_HH
#define DUMUX_CC2P_CORNERPOINT_PROBLEM_HH

#if HAVE_OPM_GRID
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/dnapl.hh>
#include <dumux/porousmediumflow/2p/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/linear/amgbackend.hh>

#include <dumux/io/cpgridcreator.hh>

#include "cc2pcornerpointspatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class CC2PCornerPointProblem;

namespace Properties
{
NEW_TYPE_TAG(CC2PCornerPointProblem, INHERITS_FROM(CCTpfaModel, TwoP, CC2PCornerPointSpatialParams));

// Set the grid type
SET_TYPE_PROP(CC2PCornerPointProblem, Grid, Dune::CpGrid);

// Set the problem property
SET_TYPE_PROP(CC2PCornerPointProblem, Problem, CC2PCornerPointProblem<TypeTag>);

// Set the grid creator
SET_TYPE_PROP(CC2PCornerPointProblem, GridCreator, CpGridCreator<TypeTag>);

// Set the fluid system
SET_PROP(CC2PCornerPointProblem, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using WettingPhase = FluidSystems::LiquidPhase<Scalar, SimpleH2O<Scalar> >;
    using NonwettingPhase = FluidSystems::LiquidPhase<Scalar, DNAPL<Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};
} // end namespace Properties

template <class TypeTag >
class CC2PCornerPointProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    // Grid and world dimension
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum {
        // primary variable indices
        pwIdx = Indices::pwIdx,
        snIdx = Indices::snIdx,

        // equation indices
        contiNEqIdx = Indices::contiNEqIdx,
    };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    CC2PCornerPointProblem(TimeManager &timeManager,
                const GridView &gridView)
    : ParentType(timeManager, gridView), gravity_(0)
    {
        gravity_[dimWorld-1] = 9.81;
        temperature_ = 273.15 + 20; // -> 20°C

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);

        injectionElement_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Problem, InjectionElement);
        injectionRate_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, InjectionRate);
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

    const GlobalPosition& gravity() const
    {
        return gravity_;
    }

    /*!
     * \brief Returns the temperature \f$ K \f$
     *
     * This problem assumes a uniform temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return temperature_; }

    /*!
     * \brief Evaluate the source term for all phases within a given
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
    PrimaryVariables source(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume &scv) const
    {
        PrimaryVariables values(0.0);

        int eIdx = GridCreator::grid().leafGridView().indexSet().index(element);
        if (eIdx == injectionElement_)
            values[contiNEqIdx] = injectionRate_/element.geometry().volume();

        return values;
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     */
    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;

        // set no-flux on top and bottom, hydrostatic on the rest
        // use the intersection normal to decide
        const GlobalPosition normal = scvf.unitOuterNormal();
        using std::abs;
        if (abs(normal[dimWorld-1]) > 0.5 - eps_)
            values.setAllNeumann();
        else
            values.setAllDirichlet();

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param values Stores the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0);

        // hydrostatic pressure
        Scalar densityW = 1000;
        values[pwIdx] = 1e5 + densityW*(this->gravity()*globalPos);
        values[snIdx] = 0.0;

        return values;
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
    PrimaryVariables neumannAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        return values;
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
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;

        // hydrostatic pressure
        Scalar densityW = 1000;
        values[pwIdx] = 1e5 + densityW*(this->gravity()*globalPos);
        values[snIdx] = 0.0;

        return values;
    }
    // \}


    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    template<class VtkOutputModule>
    void addVtkOutputFields(VtkOutputModule& outputModule) const
    {
        auto& permX = outputModule.createScalarField("PERMX [mD]", 0);
        auto& permZ = outputModule.createScalarField("PERMZ [mD]", 0);

        for (const auto& element : elements(this->gridView()))
        {
            int eIdx = this->elementMapper().index(element);
            auto fvGeometry = localView(this->model().fvGridGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                auto K = this->spatialParams().permeability(element, scv,
                        this->model().elementSolution(element, this->model().curSol()));

                // transfer output to mD = 9.86923e-16 m^2
                permX[eIdx] = K[0][0]/9.86923e-16;
                permZ[eIdx] = K[dimWorld-1][dimWorld-1]/9.86923e-16;
            }
        }
    }


private:
    Scalar temperature_;
    static constexpr Scalar eps_ = 3e-6;
    std::string name_;
    GlobalPosition gravity_;
    int injectionElement_;
    Scalar injectionRate_;
};
} //end namespace

#endif // HAVE_OPM_GRID

#endif
