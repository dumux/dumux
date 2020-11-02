// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup SequentialTwoPModel
 * \brief Model for the pressure equation discretized by mimetic FD.
 */
#ifndef DUMUX_MIMETICPRESSURE2P_HH
#define DUMUX_MIMETICPRESSURE2P_HH

#include <dune/common/exceptions.hh>

// dumux environment
#include <dumux/porousmediumflow/sequential/mimetic/properties.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/pressure.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mimetic/operator.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mimetic/mimetic.hh>

#include <dumux/common/deprecated.hh>

namespace Dumux {

/*!
 * \ingroup SequentialTwoPModel
 * \brief Mimetic method for the pressure equation.
 *
 * Provides a mimetic implementation for the evaluation
 * of equations of the form
 * \f[\text{div}\, \boldsymbol{v}_{total} = q.\f]
 * The definition of the total velocity \f$\boldsymbol{v}_total\f$ depends on the kind of pressure chosen.
 * This could be a wetting (w) phase pressure leading to
 * \f[ - \text{div}\,  \left[\lambda \boldsymbol{K} \left(\text{grad}\, p_w + f_n \text{grad}\, p_c
 *     + \sum f_\alpha \rho_\alpha g  \text{grad}\, z\right)\right] = q, \f]
 * a nonwetting (n) phase pressure yielding
 * \f[ - \text{div}\,  \left[\lambda \boldsymbol{K}  \left(\text{grad}\, p_n - f_w \text{grad}\, p_c
 *     + \sum f_\alpha \rho_\alpha g  \text{grad}\, z\right)\right] = q, \f]
 * or a global pressure leading to
 * \f[ - \text{div}\, \left[\lambda \boldsymbol{K} \left(\text{grad}\, p_{global}
 *     + \sum f_\alpha \rho_\alpha g  \text{grad}\, z\right)\right] = q.\f]
 * Here, \f$p\f$ denotes a pressure, \f$\boldsymbol{K}\f$ the absolute permeability,
 * \f$\lambda\f$ the total mobility, possibly depending on the
 * saturation,\f$f\f$ the fractional flow function of a phase, \f$\rho\f$ a phase density,
 * \f$g\f$ the gravity constant and \f$q\f$ the source term.
 * For all cases, \f$p = p_D\f$ on \f$\Gamma_{Neumann}\f$, and \f$\boldsymbol{v}_{total}  = q_N\f$
 * on \f$\Gamma_{Dirichlet}\f$.
 *
 *\tparam TypeTag The problem Type Tag
 */
template<class TypeTag> class MimeticPressure2P
{
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNw,
        pGlobal = Indices::pressureGlobal,
        sw = Indices::saturationW,
        sn = Indices::saturationNw,
        vw = Indices::velocityW,
        vn = Indices::velocityNw,
        //! gives kind of pressure used (\f$ 0 = p_w\f$, \f$ 1 = p_n\f$, \f$ 2 = p_{global}\f$)
        pressureType = getPropValue<TypeTag, Properties::PressureFormulation>(),
        //! gives kind of saturation used (\f$ 0 = S_w\f$, \f$ 1 = S_n\f$)
        saturationType = getPropValue<TypeTag, Properties::SaturationFormulation>(),
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        numPhases = getPropValue<TypeTag, Properties::NumPhases>()
    };

    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using Grid = typename GridView::Grid;

    using Geometry = typename Element::Geometry;
    using JacobianTransposed = typename Geometry::JacobianTransposed ;

    using LocalStiffness = GetPropType<TypeTag, Properties::LocalStiffness>;
    using TraceType = Dune::BlockVector<Dune::FieldVector<Scalar, 1> >;
    using OperatorAssembler = MimeticOperatorAssemblerTwoP<TypeTag>;

    using CellData = GetPropType<TypeTag, Properties::CellData>;
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using ScalarSolutionType = typename SolutionTypes::ScalarSolution;

    using Matrix = GetPropType<TypeTag, Properties::PressureCoefficientMatrix>;
    using Vector = GetPropType<TypeTag, Properties::PressureRHSVector>;

    using DimVector = Dune::FieldVector<Scalar, dim>;

    //! Initializes the matrix to store the system of equations
    void initializeMatrix();

    //! Function which assembles the system of equations to be solved
    void assemble(bool first)
    {
        Scalar timeStep = problem_.timeManager().timeStepSize();
        Scalar maxError = 0.0;
        int size = problem_.gridView().size(0);
        for (int i = 0; i < size; i++)
        {
            Scalar sat = 0;
            using std::max;
            switch (saturationType)
            {
            case sw:
                sat = problem_.variables().cellData(i).saturation(wPhaseIdx);
                break;
            case sn:
                sat = problem_.variables().cellData(i).saturation(nPhaseIdx);
                break;
            default:
                DUNE_THROW(Dune::NotImplemented, "Only saturation formulation sw and sn are implemented!");
            }
            if (sat > 1.0)
            {
                maxError = max(maxError, (sat - 1.0) / timeStep);
            }
            if (sat < 0.0)
            {
                maxError = max(maxError, (-sat) / timeStep);
            }
        }

        lstiff_.setErrorInfo(maxError, timeStep);
        A_.assemble(lstiff_, pressTrace_, f_);
        return;
    }

    //! Solves the system of equations to get the spatial distribution of the pressure
    void solve();

    void postprocess()
    {
        A_.calculatePressure(lstiff_, pressTrace_, problem_);

        return;
    }

public:
    //! Constitutive functions are initialized and stored in the variables object
    void updateMaterialLaws();

    /*!
     * \brief Initializes the model
     *
     * Function initializes the sparse matrix to solve the global system of
     * equations and sets/calculates the initial pressure.
     *
     * \param solveTwice indicates if more than one iteration is allowed to get an initial pressure solution
     */
    void initialize(bool solveTwice = true)
    {
        const auto element = *problem_.gridView().template begin<0>();
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

        updateMaterialLaws();
        A_.initialize();
        pressTrace_.resize(problem_.gridView().size(1));
        f_.resize(problem_.gridView().size(1));
        lstiff_.initialize();
        lstiff_.reset();

        pressTrace_ = 0.0;
        f_ = 0;

        assemble(true);
        solve();
        postprocess();

        return;
    }

    /*!
     * \brief Velocity update
     *
     * Reset the velocities in the cellData.
     */
    void updateVelocity()
    {
        updateMaterialLaws();
        postprocess();
    }

    // TODO doc me!
    void update()
    {
        lstiff_.reset();

        assemble(false);
        solve();
        postprocess();

        return;
    }

    /*!
     * \brief  Write data file
     *
     * \tparam MultiWriter Class defining the output writer
     * \param writer The output writer (usually a <tt>VTKMultiWriter</tt> object)
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        int size = problem_.gridView().size(0);
        ScalarSolutionType *potential = writer.allocateManagedBuffer(size);
        ScalarSolutionType *pressure = 0;
        ScalarSolutionType *pressureSecond = 0;
        ScalarSolutionType *potentialSecond = 0;
        Dune::BlockVector < DimVector > *velocityWetting = 0;
        Dune::BlockVector < DimVector > *velocityNonwetting = 0;

        if (vtkOutputLevel_ > 0)
        {
            pressure = writer.allocateManagedBuffer(size);
            pressureSecond = writer.allocateManagedBuffer(size);
            potentialSecond = writer.allocateManagedBuffer(size);
            velocityWetting = writer.template allocateManagedBuffer<Scalar, dim>(size);
            velocityNonwetting = writer.template allocateManagedBuffer<Scalar, dim>(size);
        }


            for (const auto& element : elements(problem_.gridView()))
            {
                int eIdxGlobal = problem_.variables().index(element);
                CellData& cellData = problem_.variables().cellData(eIdxGlobal);

                if (pressureType == pw)
                {
                    (*potential)[eIdxGlobal] = cellData.potential(wPhaseIdx);
                }

                if (pressureType == pn)
                {
                    (*potential)[eIdxGlobal] = cellData.potential(nPhaseIdx);
                }

                if (vtkOutputLevel_ > 0)
                {

                if (pressureType == pw)
                {
                    (*pressure)[eIdxGlobal] = cellData.pressure(wPhaseIdx);
                    (*potentialSecond)[eIdxGlobal] = cellData.potential(nPhaseIdx);
                    (*pressureSecond)[eIdxGlobal] = cellData.pressure(nPhaseIdx);
                }

                if (pressureType == pn)
                {
                    (*pressure)[eIdxGlobal] = cellData.pressure(nPhaseIdx);
                    (*potentialSecond)[eIdxGlobal] = cellData.potential(wPhaseIdx);
                    (*pressureSecond)[eIdxGlobal] = cellData.pressure(wPhaseIdx);
                }

                const typename Element::Geometry& geometry = element.geometry();

                // get corresponding reference element
                const auto refElement = referenceElement(geometry);

                const int numberOfFaces=refElement.size(1);
                std::vector<Scalar> fluxW(numberOfFaces,0);
                std::vector<Scalar> fluxNw(numberOfFaces,0);

                // run through all intersections with neighbors and boundary
                for (const auto& intersection : intersections(problem_.gridView(), element))
                {
                    int isIndex = intersection.indexInInside();

                    fluxW[isIndex] += intersection.geometry().volume()
                        * (intersection.centerUnitOuterNormal() * cellData.fluxData().velocity(wPhaseIdx, isIndex));
                    fluxNw[isIndex] += intersection.geometry().volume()
                        * (intersection.centerUnitOuterNormal() * cellData.fluxData().velocity(nPhaseIdx, isIndex));
                }

                // calculate velocity on reference element as the Raviart-Thomas-0
                // interpolant of the fluxes
                Dune::FieldVector<Scalar, dim> refVelocity;
                // simplices
                if (refElement.type().isSimplex()) {
                    for (int dimIdx = 0; dimIdx < dim; dimIdx++)
                    {
                        refVelocity[dimIdx] = -fluxW[dim - 1 - dimIdx];
                        for (int fIdx = 0; fIdx < dim + 1; fIdx++)
                        {
                            refVelocity[dimIdx] += fluxW[fIdx]/(dim + 1);
                        }
                    }
                }
                // cubes
                else if (refElement.type().isCube()){
                    for (int i = 0; i < dim; i++)
                        refVelocity[i] = 0.5 * (fluxW[2*i + 1] - fluxW[2*i]);
                }
                // 3D prism and pyramids
                else {
                    DUNE_THROW(Dune::NotImplemented, "velocity output for prism/pyramid not implemented");
                }

                const DimVector& localPos = refElement.position(0, 0);

                // get the transposed Jacobian of the element mapping
                const JacobianTransposed jacobianT = geometry.jacobianTransposed(localPos);

                // calculate the element velocity by the Piola transformation
                DimVector elementVelocity(0);
                jacobianT.umtv(refVelocity, elementVelocity);
                elementVelocity /= geometry.integrationElement(localPos);

                (*velocityWetting)[eIdxGlobal] = elementVelocity;

                // calculate velocity on reference element as the Raviart-Thomas-0
                // interpolant of the fluxes
                // simplices
                if (refElement.type().isSimplex()) {
                    for (int dimIdx = 0; dimIdx < dim; dimIdx++)
                    {
                        refVelocity[dimIdx] = -fluxNw[dim - 1 - dimIdx];
                        for (int fIdx = 0; fIdx < dim + 1; fIdx++)
                        {
                            refVelocity[dimIdx] += fluxNw[fIdx]/(dim + 1);
                        }
                    }
                }
                // cubes
                else if (refElement.type().isCube()){
                    for (int i = 0; i < dim; i++)
                        refVelocity[i] = 0.5 * (fluxNw[2*i + 1] - fluxNw[2*i]);
                }
                // 3D prism and pyramids
                else {
                    DUNE_THROW(Dune::NotImplemented, "velocity output for prism/pyramid not implemented");
                }

                // calculate the element velocity by the Piola transformation
                elementVelocity = 0;
                jacobianT.umtv(refVelocity, elementVelocity);
                elementVelocity /= geometry.integrationElement(localPos);

                (*velocityNonwetting)[eIdxGlobal] = elementVelocity;
                }
            }

            if (pressureType == pw)
            {
                writer.attachCellData(*potential, "wetting potential");
            }

            if (pressureType == pn)
            {
                writer.attachCellData(*potential, "nonwetting potential");
            }

            if (vtkOutputLevel_ > 0)
            {
            if (pressureType == pw)
            {
                writer.attachCellData(*pressure, "wetting pressure");
                writer.attachCellData(*pressureSecond, "nonwetting pressure");
                writer.attachCellData(*potentialSecond, "nonwetting potential");
            }

            if (pressureType == pn)
            {
                writer.attachCellData(*pressure, "nonwetting pressure");
                writer.attachCellData(*pressureSecond, "wetting pressure");
                writer.attachCellData(*potentialSecond, "wetting potential");
            }

            writer.attachCellData(*velocityWetting, "wetting-velocity", dim);
            writer.attachCellData(*velocityNonwetting, "nonwetting-velocity", dim);
            }
    }

    /*!
     * \brief General method for serialization, output
     *
     * Function needed for restart option.
     *
     * \param outstream The output stream
     * \param element The grid element
     */
    void serializeEntity(std::ostream &outstream, const Element &element)
    {
        int numFaces = element.subEntities(1);
        for (int i=0; i < numFaces; i++)
        {
            int fIdxGlobal = A_.faceMapper().subIndex(element, i, 1);
            outstream << pressTrace_[fIdxGlobal][0];
        }
    }

    /*!
     * \brief General method for deserialization
     *
     * \param instream The input stream
     * \param element The grid element
     */
    void deserializeEntity(std::istream &instream, const Element &element)
    {
        int numFaces = element.subEntities(1);
        for (int i=0; i < numFaces; i++)
        {
            int fIdxGlobal = A_.faceMapper().subIndex(element, i, 1);
            instream >> pressTrace_[fIdxGlobal][0];
        }
    }

    /*!
     * \brief Constructs a MimeticPressure2P object
     *
     * \param problem The Dumux problem
     */
    MimeticPressure2P(Problem& problem) :
    problem_(problem),
    A_(problem.gridView()), lstiff_(problem_, false, problem_.gridView())
    {
        if (pressureType != pw && pressureType != pn)
        {
            DUNE_THROW(Dune::NotImplemented, "Pressure type not supported!");
        }
        if (saturationType != sw && saturationType != sn)
        {
            DUNE_THROW(Dune::NotImplemented, "Saturation type not supported!");
        }
        if (getPropValue<TypeTag, Properties::EnableCompressibility>())
        {
            DUNE_THROW(Dune::NotImplemented, "Compressibility not supported!");
        }

        density_[wPhaseIdx] = 0.0;
        density_[nPhaseIdx] = 0.0;
        viscosity_[wPhaseIdx] = 0.0;
        viscosity_[nPhaseIdx] = 0.0;

        vtkOutputLevel_ = getParam<int>("Vtk.OutputLevel");
    }

private:
    Problem& problem_;
    TraceType pressTrace_; // vector of pressure traces
    TraceType f_;
    OperatorAssembler A_;
    LocalStiffness lstiff_;

    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];

    int vtkOutputLevel_;
};

//! Solves the system of equations to get the spatial distribution of the pressure.
template<class TypeTag>
void MimeticPressure2P<TypeTag>::solve()
{
    using Solver = GetPropType<TypeTag, Properties::LinearSolver>;

    auto verboseLevelSolver = getParam<int>("LinearSolver.Verbosity", 0);

    if (verboseLevelSolver)
    std::cout << "MimeticPressure2P: solve for pressure" << std::endl;

    auto solver = getSolver<Solver>(problem_);
    solver.solve(*A_, pressTrace_, f_);

    return;
}

//! Constitutive functions are updated once if new saturations are calculated and stored in the variables object.
template<class TypeTag>
void MimeticPressure2P<TypeTag>::updateMaterialLaws()
{
        // iterate through leaf grid an evaluate c0 at cell center
        for (const auto& element : elements(problem_.gridView()))
        {
            int eIdxGlobal = problem_.variables().index(element);

            CellData& cellData = problem_.variables().cellData(eIdxGlobal);

            const Scalar satW = cellData.saturation(wPhaseIdx);

            // old material law interface is deprecated: Replace this by
            // const auto& fluidMatrixInteraction = spatialParams.fluidMatrixInteractionAtPos(element.geometry().center());
            // after the release of 3.3, when the deprecated interface is no longer supported
            const auto fluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem_.spatialParams(), element);

            // initialize mobilities
            const Scalar mobilityW = fluidMatrixInteraction.krw(satW) / viscosity_[wPhaseIdx];
            const Scalar mobilityNw = fluidMatrixInteraction.krn(satW) / viscosity_[nPhaseIdx];

            // initialize mobilities
            cellData.setMobility(wPhaseIdx, mobilityW);
            cellData.setMobility(nPhaseIdx, mobilityNw);

            // initialize fractional flow functions
            cellData.setFracFlowFunc(wPhaseIdx, mobilityW / (mobilityW + mobilityNw));
            cellData.setFracFlowFunc(nPhaseIdx, mobilityNw / (mobilityW + mobilityNw));
        }
        return;
}

} // end namespace Dumux
#endif
