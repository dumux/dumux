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
 * \brief Non-isothermal gas injection problem where a gas (e.g. steam/air)
 *        is injected into a unsaturated porous medium with a residually
 *        trapped NAPL contamination.
 */
#ifndef DUMUX_KUEVETTE3P3CNIPROBLEM_HH
#define DUMUX_KUEVETTE3P3CNIPROBLEM_HH

#include <dune/common/float_cmp.hh>

#include <dumux/material/fluidsystems/h2oairmesitylene.hh>

#include <dumux/porousmediumflow/3p3c/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include "kuevettespatialparams.hh"

#define ISOTHERMAL 0

namespace Dumux
{
template <class TypeTag>
class KuevetteProblem;

namespace Properties
{
NEW_TYPE_TAG(KuevetteProblem, INHERITS_FROM(ThreePThreeCNI, KuevetteSpatialParams));
NEW_TYPE_TAG(KuevetteBoxProblem, INHERITS_FROM(BoxModel, KuevetteProblem));
NEW_TYPE_TAG(KuevetteCCProblem, INHERITS_FROM(CCModel, KuevetteProblem));

// Set the grid type
SET_TYPE_PROP(KuevetteProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(KuevetteProblem, Problem, KuevetteProblem<TypeTag>);

// Set the fluid system
SET_TYPE_PROP(KuevetteProblem,
              FluidSystem,
              FluidSystems::H2OAirMesitylene<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// set newton relative tolerance
SET_SCALAR_PROP(KuevetteProblem, NewtonMaxRelativeShift, 1e-6);
}


/*!
 * \ingroup ThreePThreeCModel
 * \ingroup ImplicitTestProblems
 * \brief Non-isothermal gas injection problem where a gas (e.g. steam/air)
 *        is injected into a unsaturated porous medium with a residually
 *        trapped NAPL contamination.
 *
 * The domain is a quasi-two-dimensional container (kuevette). Its dimensions
 * are 1.5 m x 0.74 m. The top and bottom boundaries are closed, the right
 * boundary is a Dirichlet boundary allowing fluids to escape. From the left,
 * an injection of a hot water-air mixture is applied (Neumann boundary condition
 * for the mass components and the enthalpy), aimed at remediating an initial
 * NAPL (Non-Aquoeus Phase Liquid) contamination in the heterogeneous domain.
 * The contamination is initially placed partly into the coarse sand
 * and partly into a fine sand lense.
 *
 * This simulation can be varied through assigning different boundary conditions
 * at the left boundary as described in Class (2001):
 * Theorie und numerische Modellierung nichtisothermer Mehrphasenprozesse in
 * NAPL-kontaminierten por"osen Medien, Dissertation, Eigenverlag des Instituts
 * f"ur Wasserbau
 *
 * This problem uses the \ref ThreePThreeCModel and \ref NIModel model.
 *
 * To see the basic effect and the differences to scenarios with pure steam or
 * pure air injection, it is sufficient to simulated for about 2-3 hours (10000 s).
 * Complete remediation of the domain requires much longer (about 10 days simulated time).
 * To adjust the simulation time it is necessary to edit the input file.
 *
 * To run the simulation execute:
 * <tt>./test_box3p3cnikuevette test_box3p3cnikuevette.input</tt> or
 * <tt>./test_cc3p3cnikuevette test_cc3p3cnikuevette.input</tt>
 */
template <class TypeTag >
class KuevetteProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::Grid Grid;

    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {

        pressureIdx = Indices::pressureIdx,
        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,

        // Phase State
        threePhases = Indices::threePhases,
        wgPhaseOnly = Indices::wgPhaseOnly,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };


    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    /*!
     * \brief The constructor.
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    KuevetteProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        FluidSystem::init();

        name_ = GET_RUNTIME_PARAM(TypeTag, std::string, Problem.Name);
        episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, EpisodeLength);
        this->timeManager().startNextEpisode(episodeLength_);
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
     * \brief Returns the source term at specific position in the domain.
     *
     * \param values The source values for the primary variables
     * \param globalPos The position of the center of the finite volume
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
     * \param globalPos The position for which the bc type should be evaluated
     */
    void boundaryTypesAtPos(BoundaryTypes &values,
                            const GlobalPosition &globalPos) const
    {
        if(globalPos[0] > this->bBoxMax()[0] - eps_)
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
     * \param values The neumann values for the primary variables
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumannAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        values = 0;

        // negative values for injection
        if (globalPos[0] < eps_)
        {
            values[Indices::contiWEqIdx] = -0.1435; // 0.3435 [mol/(s m)] in total
            values[Indices::contiGEqIdx] = -0.2;
            values[Indices::contiNEqIdx] =  0.0;
            values[Indices::energyEqIdx] = -6929.;
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
     * \param vIdxGlobal The global index of the vertex
     * \param globalPos The global position
     */
    int initialPhasePresence(const Vertex &vertex,
                             const int &vIdxGlobal,
                             const GlobalPosition &globalPos) const
    {
        if (isInContaminationZone(globalPos))
            return threePhases;
        else
            return wgPhaseOnly;
    }

       /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer. Adjust this in case of anisotropic permeabilities.
     */
    void addOutputVtkFields()
    {
        // get the number of degrees of freedom
        unsigned numDofs = this->model().numDofs();

        // create the scalar field required for the permeabilities
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
        ScalarField *Kxx = this->resultWriter().allocateManagedBuffer(numDofs);

        FVElementGeometry fvGeometry;

        for (const auto& element : elements(this->gridView()))
        {
            fvGeometry.update(this->gridView(), element);

            for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
            {
                int dofIdxGlobal = this->model().dofMapper().subIndex(element, scvIdx, dofCodim);
                (*Kxx)[dofIdxGlobal] = this->spatialParams().intrinsicPermeability(element, fvGeometry, scvIdx);
            }
        }

        this->resultWriter().attachDofData(*Kxx, "permeability", isBox);
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
    // checks, whether a point is located inside the contamination zone
    bool isInContaminationZone(const GlobalPosition &globalPos) const
    {
        return (Dune::FloatCmp::ge<Scalar>(globalPos[0],0.2)
               && Dune::FloatCmp::le<Scalar>(globalPos[0],0.8)
               && Dune::FloatCmp::ge<Scalar>(globalPos[1],0.4)
               && Dune::FloatCmp::le<Scalar>(globalPos[1],0.65));
    }

    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        values[pressureIdx] = 1e5 ;
        values[switch1Idx] = 0.12;
        values[switch2Idx] = 1.e-6;
        values[temperatureIdx] = 293.0;

        if (isInContaminationZone(globalPos))
        {
            values[switch2Idx] = 0.07;
        }
    }

    static constexpr Scalar eps_ = 1e-6;
    std::string name_;
    Scalar episodeLength_;
};
} //end namespace

#endif
