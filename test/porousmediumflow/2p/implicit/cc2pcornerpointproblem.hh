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

#if HAVE_DUNE_CORNERPOINT
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/dnapl.hh>
#include <dumux/porousmediumflow/2p/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/implicit/cellcentered/propertydefaults.hh>
#include <dumux/linear/amgbackend.hh>

#include <dune/grid/CpGrid.hpp>
#include <dumux/io/cpgridcreator.hh>
#include <dumux/implicit/cornerpoint/elementvolumevariables.hh>
#include <dumux/implicit/cornerpoint/fvelementgeometry.hh>
#include <dumux/porousmediumflow/implicit/cpdarcyfluxvariables.hh>

#include "cc2pcornerpointspatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class CC2PCornerPointProblem;

namespace Properties
{
NEW_TYPE_TAG(CC2PCornerPointProblem, INHERITS_FROM(CCTwoP, CC2PCornerPointSpatialParams));

// Set the grid type
SET_TYPE_PROP(CC2PCornerPointProblem, Grid, Dune::CpGrid);

// Set the problem property
SET_TYPE_PROP(CC2PCornerPointProblem, Problem, Dumux::CC2PCornerPointProblem<TypeTag>);

// Set the grid creator
SET_TYPE_PROP(CC2PCornerPointProblem, GridCreator, Dumux::CpGridCreator<TypeTag>);

// Set properties that are specific for CpGrid
SET_TYPE_PROP(CC2PCornerPointProblem, ElementVolumeVariables, CpElementVolumeVariables<TypeTag>);
SET_TYPE_PROP(CC2PCornerPointProblem, FVElementGeometry, CpFVElementGeometry<TypeTag>);
SET_TYPE_PROP(CC2PCornerPointProblem, FluxVariables, CpDarcyFluxVariables<TypeTag>);

// Set the wetting phase
SET_PROP(CC2PCornerPointProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(CC2PCornerPointProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::DNAPL<Scalar> > type;
};

// Linear solver settings
SET_TYPE_PROP(CC2PCornerPointProblem, LinearSolver, Dumux::AMGBackend<TypeTag> );
}

template <class TypeTag >
class CC2PCornerPointProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;

    enum {

        // primary variable indices
        pwIdx = Indices::pwIdx,
        snIdx = Indices::snIdx,

        // equation indices
        contiNEqIdx = Indices::contiNEqIdx,

        // phase indices
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,


        // world dimension
        dimWorld = GridView::dimensionworld
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
    CC2PCornerPointProblem(TimeManager &timeManager,
                const GridView &gridView)
    : ParentType(timeManager, gridView), gravity_(0)
    {
        gravity_[dimWorld-1] = 9.81;
        eps_ = 3e-6;
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
    const char *name() const
    {
        return name_.c_str();
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
     * \param values The source and sink values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^3 \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scvIdx The local subcontrolvolume index
     *
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void source(PrimaryVariables &values,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int scvIdx) const
    {
        values = 0.0;

        int eIdx = GridCreator::grid().leafGridView().indexSet().index(element);
        if (eIdx == injectionElement_)
            values[contiNEqIdx] = injectionRate_/element.geometry().volume();
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
     * \param intersection The intersection
     */
    void boundaryTypes(BoundaryTypes &values,
                       const Intersection &intersection) const
    {
        // set no-flux on top and bottom, hydrostatic on the rest
        // use the intersection normal to decide
        const GlobalPosition normal = intersection.centerUnitOuterNormal();
        if (std::abs(normal[dimWorld-1]) > 0.5)
            values.setAllNeumann();
        else
            values.setAllDirichlet();
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
        // hydrostatic pressure
        Scalar densityW = 1000;
        values[pwIdx] = 1e5 + densityW*(this->gravity()*globalPos);
        values[snIdx] = 0.0;
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
        // hydrostatic pressure
        Scalar densityW = 1000;
        values[pwIdx] = 1e5 + densityW*(this->gravity()*globalPos);
        values[snIdx] = 0.0;
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
        ScalarField *permZ = this->resultWriter().allocateManagedBuffer(numElements);

        for (const auto& element : Dune::elements(this->gridView()))
        {
            FVElementGeometry fvGeometry;
            fvGeometry.update(this->gridView(), element);
            int eIdx = this->elementMapper().index(element);

            const DimWorldMatrix K = this->spatialParams().intrinsicPermeability(element, fvGeometry, /*element data*/ 0);
            // transfer output to mD = 9.86923e-16 m^2
            (*permX)[eIdx] = K[0][0]/9.86923e-16;
            (*permZ)[eIdx] = K[2][2]/9.86923e-16;
        }

        //pass the scalar fields to the vtkwriter
        this->resultWriter().attachDofData(*permX, "PERMX [mD]", false); //element data
        this->resultWriter().attachDofData(*permZ, "PERMZ [mD]", false); //element data
    }

private:
    Scalar temperature_;
    Scalar eps_;
    std::string name_;
    GlobalPosition gravity_;
    int injectionElement_;
    Scalar injectionRate_;
};
} //end namespace

#endif // HAVE_DUNE_CORNERPOINT

#endif
