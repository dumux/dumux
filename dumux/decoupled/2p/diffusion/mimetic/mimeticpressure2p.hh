// $Id$
/*****************************************************************************
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Copyright (C) 2007-2009 by Jochen Fritz                                 *
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
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
#include "dumux/common/pardiso.hh"
#include <dumux/decoupled/2p/2pproperties.hh>

namespace Dumux
{

/*! \ingroup Mimetic2p
 *
 * \brief mimetic method for the pressure equation
 *
 * Provides a mimetic implementation for the evaluation
 * of equations of the form
 * \f[\text{div}\, \boldsymbol{v}_{total} = q.\f]
 * The definition of the total velocity \f$\boldsymbol{v}_total\f$ depends on the kind of pressure chosen. This could be a wetting (w) phase pressure leading to
 * \f[ - \text{div}\,  \left[\lambda \boldsymbol{K} \left(\text{grad}\, p_w + f_n \text{grad}\, p_c + \sum f_\alpha \rho_\alpha g  \text{grad}\, z\right)\right] = q, \f]
 * a non-wetting (n) phase pressure yielding
 * \f[ - \text{div}\,  \left[\lambda \boldsymbol{K}  \left(\text{grad}\, p_n - f_w \text{grad}\, p_c + \sum f_\alpha \rho_\alpha g  \text{grad}\, z\right)\right] = q, \f]
 * or a global pressure leading to
 * \f[ - \text{div}\, \left[\lambda \boldsymbol{K} \left(\text{grad}\, p_{global} + \sum f_\alpha \rho_\alpha g  \text{grad}\, z\right)\right] = q.\f]
 *  Here, \f$p\f$ denotes a pressure, \f$\boldsymbol{K}\f$ the absolute permeability, \f$\lambda\f$ the total mobility, possibly depending on the
 * saturation,\f$f\f$ the fractional flow function of a phase, \f$\rho\f$ a phase density, \f$g\f$ the gravity constant and \f$q\f$ the source term.
 * For all cases, \f$p = p_D\f$ on \f$\Gamma_{Neumann}\f$, and \f$\boldsymbol{v}_{total}  = q_N\f$
 * on \f$\Gamma_{Dirichlet}\f$.
 *
 *\tparam TypeTag The Type Tag
 */
template<class TypeTag> class MimeticPressure2P
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Variables)) Variables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

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
        Sn = Indices::saturationNW,
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalStiffness)) LocalStiffness;
    typedef Dune::BlockVector< Dune::FieldVector<Scalar, 1> > TraceType;
    typedef Dune::BlockVector< Dune::FieldVector<Scalar, 2*dim> > NormalVelType;
    typedef MimeticOperatorAssembler<Scalar,GridView> OperatorAssembler;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureCoefficientMatrix)) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureRHSVector)) Vector;

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
        A_.calculatePressure(lstiff, pressTrace_, normalVelocity_, problem_.variables().pressure());
        return;
    }

protected:
    Problem& problem()
    {
        return problem_;
    }

    const Problem& problem() const
    {
        return problem_;
    }

public:
    //constitutive functions are initialized and stored in the variables object
    void updateMaterialLaws();

    void initialize(bool solveTwice = true)
    {
        updateMaterialLaws();
        assemble(true);
        solve();
        postprocess();

        return;
    }

    void pressure(bool solveTwice = true)
    {
        assemble(false);
        solve();
        postprocess();

        return;
    }

    void calculateVelocity();

    // serialization methods
    template<class Restarter>
    void serialize(Restarter &res)
    {
       return;
    }

    template<class Restarter>
    void deserialize(Restarter &res)
    {
        return;
    }

    //! \brief Write data files
     /*  \param name file name */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        typename Variables::ScalarSolutionType *pressure = writer.template createField<Scalar, 1> (problem_.gridView().size(0));

        *pressure = problem_.variables().pressure();

        writer.addCellData(pressure, "global pressure");

        return;
    }

    //! Constructs a MimeticPressure2P object
    /**
     * \param problem The Dumux problem
     */
    MimeticPressure2P(Problem& problem) :
    problem_(problem),
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
    }

private:
    Problem& problem_;
    TraceType pressTrace_; //!< vector of pressure traces
    NormalVelType normalVelocity_;
    TraceType f_;
    OperatorAssembler A_;
protected:
    static const int pressureType = GET_PROP_VALUE(TypeTag, PTAG(PressureFormulation)); //!< gives kind of pressure used (\f$ 0 = p_w\f$, \f$ 1 = p_n\f$, \f$ 2 = p_{global}\f$)
    static const int saturationType = GET_PROP_VALUE(TypeTag, PTAG(SaturationFormulation)); //!< gives kind of saturation used (\f$ 0 = S_w\f$, \f$ 1 = S_n\f$)
};

//solves the system of equations to get the spatial distribution of the pressure
template<class TypeTag>
void MimeticPressure2P<TypeTag>::solve()
{
    typedef typename GET_PROP(TypeTag, PTAG(SolverParameters)) SolverParameters;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressurePreconditioner)) Preconditioner;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureSolver)) Solver;

    typedef TraceType Vector;
    typedef typename CROperatorAssembler<Scalar,GridView>::RepresentationType Matrix;
    typedef Dune::MatrixAdapter<Matrix,Vector,Vector> Operator;

    Operator op(*A_);
    Dune::InverseOperatorResult result;

    double reduction = SolverParameters::reductionSolver;
    int maxItSolver = SolverParameters::maxIterationNumberSolver;
    int iterPreconditioner = SolverParameters::iterationNumberPreconditioner;
    int verboseLevelSolver = SolverParameters::verboseLevelSolver;
    double relaxation = SolverParameters::relaxationPreconditioner;

    if (verboseLevelSolver)
    std::cout << "MimeticPressure2P: solve for pressure" << std::endl;

    Preconditioner preconditioner(*A_, iterPreconditioner, relaxation);
    Solver solver(op, preconditioner, reduction, maxItSolver, verboseLevelSolver);
    solver.apply(pressTrace_, f_, result);

    return;
}

//constitutive functions are updated once if new saturations are calculated and stored in the variables object
template<class TypeTag>
void MimeticPressure2P<TypeTag>::updateMaterialLaws()
{
    FluidState fluidState;

    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = problem_.gridView().template end<0>();
    for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // get global coordinate of cell center
        GlobalPosition globalPos = eIt->geometry().center();

        int globalIdx = problem_.variables().index(*eIt);

        Scalar sat = problem_.variables().saturation()[globalIdx];

        std::vector<Scalar> mobilities(2);

        problem_.variables().capillaryPressure(globalIdx)= MaterialLaw::pC(
                problem_.spatialParameters().materialLawParams(globalPos, *eIt), sat);

        Scalar temperature = problem_.temperature(globalPos, *eIt);
        Scalar pressW =  problem_.referencePressure(globalPos, *eIt);
        Scalar pressN = pressW;
        fluidState.update(sat, pressW, pressN, temperature);
        Scalar viscosityW = FluidSystem::phaseViscosity(wPhaseIdx, temperature, pressW, fluidState);
        Scalar viscosityN = FluidSystem::phaseViscosity(nPhaseIdx, temperature, pressN, fluidState);
        Scalar mobilityW = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(globalPos, *eIt), sat)
                / viscosityW;
        Scalar mobilityN = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(globalPos, *eIt), sat)
                / viscosityN;

        // initialize mobilities
        problem_.variables().mobilityWetting(globalIdx)= mobilityW;
        problem_.variables().mobilityNonwetting(globalIdx)= mobilityN;
        Scalar mobilityT = mobilityW + mobilityN;
        problem_.variables().fracFlowFuncWetting(globalIdx)= mobilityW/mobilityT;
        problem_.variables().fracFlowFuncNonwetting(globalIdx)= mobilityN/mobilityT;
    }
    return;
}

template<class TypeTag>
void MimeticPressure2P<TypeTag>::calculateVelocity()
{
    // ASSUMES axiparallel grids in 2D
    for (int i = 0; i < problem_.gridView().size(0); i++)
    {
        problem_.variables().velocity()[i][0][0] = -normalVelocity_[i][0];
        problem_.variables().velocity()[i][0][1] = 0;
        problem_.variables().velocity()[i][1][0] = normalVelocity_[i][1];
        problem_.variables().velocity()[i][1][1] = 0;
        problem_.variables().velocity()[i][2][0] = 0;
        problem_.variables().velocity()[i][2][1] = -normalVelocity_[i][2];
        problem_.variables().velocity()[i][3][0] = 0;
        problem_.variables().velocity()[i][3][1] = normalVelocity_[i][3];
    }
//    printvector(std::cout, problem_.variables().velocity(), "velocity", "row", 4, 1, 3);
    return;
}

}
#endif
