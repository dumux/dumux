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
 * \brief Definition of a problem, where air is injected under a low permeable layer.
 */
#ifndef DUMUX_INJECTION_PROBLEM_HH
#define DUMUX_INJECTION_PROBLEM_HH

#include <dune/grid/io/file/dgfparser/dgfs.hh>

#include <dumux/implicit/2p2cni/2p2cnimodel.hh>
#include <dumux/implicit/common/implicitporousmediaproblem.hh>
#include <dumux/material/fluidsystems/cs2n2fluidsystem.hh>

#include "injectionspatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class InjectionProblem;

namespace Properties
{
NEW_TYPE_TAG(InjectionProblem, INHERITS_FROM(TwoPTwoCNI, InjectionSpatialParams));
NEW_TYPE_TAG(InjectionBoxProblem, INHERITS_FROM(BoxModel, InjectionProblem));
NEW_TYPE_TAG(InjectionCCProblem, INHERITS_FROM(CCModel, InjectionProblem));
    
// Set the grid type
SET_PROP(InjectionProblem, Grid)
{
    typedef Dune::SGrid<2,2> type;
};

// Set the problem property
SET_PROP(InjectionProblem, Problem)
{
    typedef Dumux::InjectionProblem<TypeTag> type;
};

// Set fluid configuration
SET_PROP(InjectionProblem, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    static const bool useComplexRelations = false;
 public:
    typedef Dumux::FluidSystems::CS2N2<Scalar, useComplexRelations> type;
};

// Enable gravity
SET_BOOL_PROP(InjectionProblem, ProblemEnableGravity, true);

SET_BOOL_PROP(InjectionProblem, ImplicitEnableJacobianRecycling, true);
SET_BOOL_PROP(InjectionProblem, VtkAddVelocity, false);
}


/*!
 * \ingroup TwoPTwoCModel
 * \ingroup ImplicitTestProblems
 * \brief Problem where air is injected under a low permeable layer in a depth of 2700m.
 *
 * The domain is sized 60m times 40m and consists of two layers, a moderately
 * permeable spatial parameters (\f$ K=10e-12\f$) for \f$ y>22m\f$ and one with a lower permeablility (\f$ K=10e-13\f$)
 * in the rest of the domain.
 *
 * Air enters a water-filled aquifer, which is situated 2700m below sea level, at the right boundary
 * (\f$ 5m<y<15m\f$) and migrates upwards due to buoyancy. It accumulates and
 * partially enters the lower permeable aquitard.
 * 
 * This problem uses the \ref TwoPTwoCModel.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2p2c -parameterFile ./test_box2p2c.input</tt> or
 * <tt>./test_cc2p2c -parameterFile ./test_cc2p2c.input</tt>
 */
template <class TypeTag>
class InjectionProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,


        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,

        conti0EqIdx = Indices::conti0EqIdx,
        contiN2EqIdx = Indices::contiNEqIdx,
        energyEqIdx = Indices::energyEqIdx
    };


    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    
    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    InjectionProblem(TimeManager &timeManager,
                     const GridView &gridView)
        : ParentType(timeManager, GridCreator::grid().leafView())
    {
        try
        {
            nTemperature_       = GET_RUNTIME_PARAM(TypeTag, int, Problem.NTemperature);
            nPressure_          = GET_RUNTIME_PARAM(TypeTag, int, Problem.NPressure);
            pressureLow_        = GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.PressureLow);
            pressureHigh_       = GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.PressureHigh);
            temperatureLow_     = GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.TemperatureLow);
            temperatureHigh_    = GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.TemperatureHigh);
            temperature_        = GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.InitialTemperature);
            name_               = GET_RUNTIME_PARAM(TypeTag, std::string, Problem.Name);
        }
        catch (Dumux::ParameterException &e) {
            std::cerr << e << ". Abort!\n";
            exit(1) ;
        }
        catch (...) {
            std::cerr << "Unknown exception thrown!\n";
            exit(1);
        }

        /* Alternative syntax:
         * typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
         * const Dune::ParameterTree &tree = ParameterTree::tree();
         * nTemperature_       = tree.template get<int>("Problem.NTemperature");
         *
         * + We see what we do
         * - Reporting whether it was used does not work
         * - Overwriting on command line not possible
         */


        eps_ = 1e-6;

        // initialize the tables of the fluid system
        //FluidSystem::init();
        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);


       // calculate the injection volume
             injVolume_ = 0.0;
             FVElementGeometry fvElemGeom;
             ElementIterator elemIt = this->gridView().template begin<0>();
             const ElementIterator endIt = this->gridView().template end<0>();
               for (; elemIt != endIt; ++ elemIt) {
                   fvElemGeom.update(this->gridView(), *elemIt);
                   const GlobalPosition &pos = (*elemIt).geometry().center();

                               if (pos[1] >3.95 && pos[1] <4.05)
                                {
                                         if (elemIt->partitionType() != Dune::InteriorEntity)
                                                continue;
                                             injVolume_ += (*elemIt).geometry().volume();
                                }

                }

    }

    /*!
     * \brief Called directly after the time integration.
     */
    void postTimeStep()
    {
        // Calculate storage terms
        PrimaryVariables storageW, storageN;
        this->model().globalPhaseStorage(storageW, wPhaseIdx);
        this->model().globalPhaseStorage(storageN, nPhaseIdx);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout<<"Storage: wetting=[" << storageW << "]"
                     << " nonwetting=[" << storageN << "]\n";
        }
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
    const std::string name() const
    { return name_; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return temperature_; };

    void sourceAtPos(PrimaryVariables &values,
                     const GlobalPosition &globalPos) const
    {
          values = Scalar(0.0);

            Scalar m1, m2, m3;
            m1=5.13e-6/injVolume_;
            m2=8.5e-6/injVolume_;
            m3=0.15/injVolume_;

            if (this->timeManager().time()<4500.)
            {
                if (globalPos[1] >3.95 && globalPos[1] <4.05)
                        {
                         values[contiN2EqIdx]=m1;  // N2
                         values[conti0EqIdx]=m2;  // CS2
                         values[energyEqIdx]=m3;   // Enthalpy
                    }
            }
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
     * \param globalPos The position for which the bc type should be evaluated
     */
    void boundaryTypesAtPos(BoundaryTypes &values, 
                            const GlobalPosition &globalPos) const
    {
        if (globalPos[1] < eps_)
            values.setAllDirichlet();
        else if (globalPos[1] > (6.-eps_))
            values.setAllDirichlet();
        else
            values.setAllNeumann();
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values The neumann values for the conservation equations
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the box scheme
     * \param is The intersection between element and boundary
     * \param scvIdx The local vertex index
     * \param boundaryFaceIdx The index of the boundary face
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &is,
                 int scvIdx,
                 int boundaryFaceIdx) const
    {
        values = 0;
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
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initialAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }

    /*!
     * \brief Return the initial phase state inside a control volume.
     *
     * \param vertex The vertex
     * \param globalIdx The index of the global vertex
     * \param globalPos The global position
     */
    int initialPhasePresence(const Vertex &vertex,
                             int &globalIdx,
                             const GlobalPosition &globalPos) const
    { return Indices::nPhaseOnly; }

    // \}

private:
    // the internal method for the initial condition
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {

        values[Indices::pressureIdx] = 1.013e5;
        values[Indices::switchIdx] = 0.00001;
        values[Indices::temperatureIdx] = temperature_;
    }

    Scalar temperature_;
    Scalar depthBOR_;
    Scalar eps_;

    Scalar injVolume_;

    int nTemperature_;
    int nPressure_;

    std::string name_ ;

    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;




};
} //end namespace

#endif
