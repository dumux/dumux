// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Copyright (C) 2007-2009 by Jochen Fritz                                 *
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
 * \brief Model for the pressure equation discretized by mimetic FD.
 */
#ifndef DUMUX_MIMETICPRESSURE2P_HH
#define DUMUX_MIMETICPRESSURE2P_HH

// dune environent:
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

// dumux environment
#include <dumux/decoupled/2p/diffusion/diffusionproperties2p.hh>
#include <dumux/decoupled/common/mimetic/mimeticproperties.hh>
#include <dumux/decoupled/2p/diffusion/mimetic/mimeticoperator.hh>

namespace Dumux
{

//! \ingroup MimeticPressure2p
/*! \brief Mimetic finite differences discretization of a two-phase pressure equation of the sequential IMPES model.
 *
 * This class provides a mimetic finite differences implementation for solving equations of the form
 * \f[
 * \text{div}\, \boldsymbol{v}_{total} = q.
 * \f]
 * The total velocity \f$\boldsymbol{v}_{total}\f$ is defined using a global pressure approach. This leads to
 * \f[
 * - \text{div}\, \left(\lambda \boldsymbol K \text{grad}\, p_{global}\right) = q.
 * \f]
 * Here, \f$ p_{global} \f$ is the global pressure of a classical fractional flow formulation
 * (see e.g. ﻿P. Binning and M. A. Celia, “Practical implementation of the fractional flow approach to multi-phase flow simulation,” Advances in water resources, vol. 22, no. 5, pp. 461-478, 1999.),
 * \f$ \boldsymbol K \f$ the absolute permeability, \f$ \lambda = \lambda_w + \lambda_n \f$ the total mobility depending on the
 * saturation (\f$ \lambda_\alpha = k_{r_\alpha} / \mu_\alpha \f$) and \f$ q \f$ the source term. Gravity is neglected in this implementation.
 *
 * \f$ p = p_D \f$ on \f$ \Gamma_{Dirichlet} \f$, and \f$ \boldsymbol v_{total} \cdot  \boldsymbol n  = q_N \f$ on \f$ \Gamma_{Neumann} \f$.
 *
 * Remark: gravity is neglected!
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag> class MimeticPressure2P
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Variables) Variables;

    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, TwoPIndices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNW,
        pglobal = Indices::pressureGlobal,
        Sw = Indices::saturationW,
        Sn = Indices::saturationNW
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, LocalStiffness) LocalStiffness;
    typedef Dune::BlockVector< Dune::FieldVector<Scalar, 1> > TraceType;
    typedef Dune::BlockVector< Dune::FieldVector<Scalar, 2*dim> > NormalVelType;
    typedef MimeticOperatorAssembler<Scalar,GridView> OperatorAssembler;
    typedef typename GET_PROP(TypeTag, SolutionTypes)::ScalarSolution ScalarSolution;

    typedef typename GET_PROP_TYPE(TypeTag, PressureCoefficientMatrix) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, PressureRHSVector) Vector;

    //initializes the matrix to store the system of equations
    void initializeMatrix();

    //function which assembles the system of equations to be solved
    void assemble(bool first)
    {
        LocalStiffness lstiff(problem_, false, problem_.gridView());
        A_.assemble(lstiff, pressTrace_, f_);
        return;
    }

    //solves the system of equations to get the spatial distribution of the pressure
    void solve();

    void postprocess()
    {
        LocalStiffness lstiff(problem_, false, problem_.gridView());
        A_.calculatePressure(lstiff, pressTrace_, normalVelocity_, pressure_);
        return;
    }

protected:
    //! \cond \private
    Problem& problem()
    {
        return problem_;
    }

    const Problem& problem() const
    {
        return problem_;
    }
    //! \endcond

public:
    //constitutive functions are initialized and stored in the variables object
    void updateMaterialLaws();


    /*! \brief Initializes the pressure model
     *
     * Initializes pressure and velocity field.
     *
     * \param solveTwice indicates if more than one iteration is allowed to get an initial pressure solution
     */
    void initialize(bool solveTwice = true)
    {
        updateMaterialLaws();
        A_.initializeMatrix();
        f_.resize(problem_.gridView().size(1));//resize to make sure the final grid size (after the problem was completely built) is used!
        pressure_.resize(problem_.gridView().size(0));
        pressTrace_.resize(problem_.gridView().size(1));
        pressure_ = 0;
        pressTrace_ = 0;
        f_ = 0;
        assemble(true);
        solve();
        postprocess();

        storePressureSolution();
        storeVelocity();

        return;
    }

    /*! \brief Pressure and velocity update */
    void update()
    {
        assemble(false);
        solve();
        postprocess();

        storePressureSolution();
        storeVelocity();

        return;
    }

    /*! \brief Globally stores the pressure solution*/
    void storePressureSolution()
    {
        int size = problem_.gridView().size(0);
        for (int i = 0; i < size; i++)
        {
            CellData& cellData = problem_.variables().cellData(i);
            storePressureSolution(i, cellData);
            cellData.fluxData().resetVelocity();
        }
    }

    /*! \brief Stores the pressure solution of a cell
     *
     * \param globalIdx Global cell index
     * \param cellData A CellData object
     */
    void storePressureSolution(int globalIdx, CellData& cellData)
    {
        cellData.setGlobalPressure(pressure_[globalIdx]);
    }

    // Globally stores the velocity solution at cell-cell interfaces
    void storeVelocity();


    /*! \brief  Function for serialization of the pressure field.
     *
     *  Function needed for restart option. Writes the pressure of a grid element to a restart file.
     *
     *  \param outstream Stream into the restart file.
     *  \param element Grid element
     */
    void serializeEntity(std::ostream &outstream, const Element &element)
    {
        int globalIdx = problem_.variables().index(element);
        outstream << pressure_[globalIdx][0];
    }

    /*! \brief  Function for deserialization of the pressure field.
     *
     *  Function needed for restart option. Reads the pressure of a grid element from a restart file.
     *
     *  \param instream Stream from the restart file.
     *  \param element Grid element
     */
    void deserializeEntity(std::istream &instream, const Element &element)
    {
        int globalIdx = problem_.variables().index(element);
        instream >> pressure_[globalIdx][0];
    }

    /*! \brief Adds pressure output to the output file
     *
     * Adds the global pressure to the output.
     *
     * \tparam MultiWriter Class defining the output writer
     * \param writer The output writer (usually a <tt>VTKMultiWriter</tt> object)
     *
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        ScalarSolution *pressure = writer.allocateManagedBuffer (problem_.gridView().size(0));

        *pressure = pressure_;

        writer.attachCellData(*pressure, "global pressure");

        return;
    }

    //! Constructs a MimeticPressure2P object
    /**
     * \param problem A problem class object
     */
    MimeticPressure2P(Problem& problem) :
    problem_(problem),
    pressure_(problem.gridView().size(0)),
    pressTrace_(problem.gridView().size(1)),
    normalVelocity_(problem.gridView().size(0)),
    f_(problem.gridView().size(1)),
    A_(problem.gridView())
    {
        if (pressureType != pglobal)
        {
            DUNE_THROW(Dune::NotImplemented, "Pressure type not supported!");
        }
        if (saturationType != Sw)
        {
            DUNE_THROW(Dune::NotImplemented, "Saturation type not supported!");
        }

        const Element& element = *(problem_.gridView().template begin<0> ());
        FluidState fluidState;
        fluidState.setPressure(wPhaseIdx, problem_.referencePressure(element));
        fluidState.setPressure(nPhaseIdx, problem_.referencePressure(element));
        fluidState.setTemperature(problem_.temperature(element));
        fluidState.setSaturation(wPhaseIdx, 1.);
        fluidState.setSaturation(nPhaseIdx, 0.);
        density_[wPhaseIdx] = FluidSystem::density(fluidState, wPhaseIdx);
        density_[nPhaseIdx] = FluidSystem::density(fluidState, nPhaseIdx);
        viscosity_[wPhaseIdx] = FluidSystem::viscosity(fluidState, wPhaseIdx);
        viscosity_[nPhaseIdx] = FluidSystem::viscosity(fluidState, nPhaseIdx);
    }

private:
    Problem& problem_;
    ScalarSolution pressure_;
    TraceType pressTrace_; //!< vector of pressure traces
    NormalVelType normalVelocity_;
    TraceType f_;
    OperatorAssembler A_;

    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];
    static const int pressureType = GET_PROP_VALUE(TypeTag, PressureFormulation); //!< gives kind of pressure used (\f$ 0 = p_w \f$, \f$ 1 = p_n \f$, \f$ 2 = p_{global} \f$)
    static const int saturationType = GET_PROP_VALUE(TypeTag, SaturationFormulation); //!< gives kind of saturation used (\f$ 0 = S_w \f$, \f$ 1 = S_n \f$)
};

//solves the system of equations to get the spatial distribution of the pressure
template<class TypeTag>
void MimeticPressure2P<TypeTag>::solve()
{
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolver) Solver;

    int verboseLevelSolver = GET_PARAM_FROM_GROUP(TypeTag, int, LinearSolver, Verbosity);

    if (verboseLevelSolver)
    std::cout << "MimeticPressure2P: solve for pressure" << std::endl;

    Solver solver(problem_);
    solver.solve(*A_, pressTrace_, f_);
    return;
}

/*! \brief Updates constitutive relations and stores them in the variable class
 *
 * Stores mobility, fractional flow function and capillary pressure for all grid cells.
 */
template<class TypeTag>
void MimeticPressure2P<TypeTag>::updateMaterialLaws()
{
    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = problem_.gridView().template end<0>();
    for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eItEnd; ++eIt)
    {
        int globalIdx = problem_.variables().index(*eIt);

        CellData& cellData = problem_.variables().cellData(globalIdx);

        //determine phase saturations from primary saturation variable

        Scalar satW = cellData.saturation(wPhaseIdx);
        Scalar pc = MaterialLaw::pC(problem_.spatialParameters().materialLawParams(*eIt), satW);

        cellData.setCapillaryPressure(pc);

        // initialize mobilities
        Scalar mobilityW = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(*eIt), satW) / viscosity_[wPhaseIdx];
        Scalar mobilityNW = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(*eIt), satW) / viscosity_[nPhaseIdx];

        // initialize mobilities
        cellData.setMobility(wPhaseIdx, mobilityW);
        cellData.setMobility(nPhaseIdx, mobilityNW);

        //initialize fractional flow functions
        cellData.setFracFlowFunc(wPhaseIdx, mobilityW / (mobilityW + mobilityNW));
        cellData.setFracFlowFunc(nPhaseIdx, mobilityNW / (mobilityW + mobilityNW));
    }
    return;
}

/*! \brief Globally stores the velocity solution at cell-cell interfaces*/
template<class TypeTag>
void MimeticPressure2P<TypeTag>::storeVelocity()
{
    // iterate through leaf grid an evaluate c0 at cell center
    const ElementIterator &eItEnd = problem_.gridView().template end<0>();
    for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eItEnd; ++eIt)
    {
        int globalIdx = problem_.variables().index(*eIt);

        CellData& cellData = problem_.variables().cellData(globalIdx);

        IntersectionIterator isIt = problem_.gridView().template ibegin(*eIt);
        const IntersectionIterator &isItEnd = problem_.gridView().template iend(*eIt);
        for (; isIt != isItEnd; ++isIt)
        {
            int idxInInside = isIt->indexInInside();
            Dune::FieldVector<Scalar, dimWorld> velocity(isIt->centerUnitOuterNormal());
            velocity*= normalVelocity_[globalIdx][idxInInside];
            cellData.fluxData().setVelocity(wPhaseIdx, idxInInside, velocity);
            cellData.fluxData().setPotential(wPhaseIdx, idxInInside, normalVelocity_[globalIdx][idxInInside]);
            cellData.fluxData().setPotential(nPhaseIdx, idxInInside, normalVelocity_[globalIdx][idxInInside]);
        }
    }
//    printvector(std::cout, problem_.variables().velocity(), "velocity", "row", 4, 1, 3);
    return;
}

}
#endif
