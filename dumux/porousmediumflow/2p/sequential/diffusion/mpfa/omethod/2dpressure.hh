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
 * \brief  Finite volume MPFA O-method discretization of a two-phase pressure equation of the sequential IMPES model.
 */
#ifndef DUMUX_FVMPFAO2DPRESSURE2P_HH
#define DUMUX_FVMPFAO2DPRESSURE2P_HH

// dumux environment
#include <dumux/porousmediumflow/sequential/cellcentered/pressure.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/mpfa/ointeractionvolume.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/properties.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/mpfa/properties.hh>

#include <dumux/common/deprecated.hh>

namespace Dumux {
/*!
 * \ingroup SequentialTwoPModel
 * \brief Finite volume MPFA O-method discretization of a two-phase flow pressure equation of the sequential IMPES model.
 *
 * Finite volume MPFA O-method discretization of the equations
 * \f[ - \text{div}\, \boldsymbol v_t = - \text{div}\, (\lambda_t \boldsymbol K \textbf{grad}\,
 * \Phi_w + f_n \lambda_t \boldsymbol K \textbf{grad}\, \Phi_{cap}   ) = 0, \f]
 * or
 * \f[ - \text{div}\, \boldsymbol v_t = - \text{div}\, (\lambda_t \boldsymbol K \textbf{grad}\,
 * \Phi_n - f_w \lambda_t \boldsymbol K \textbf{grad}\, \Phi_{cap}   ) = 0. \f]
 *  At Dirichlet boundaries a two-point flux approximation is used.
 * \f[ \Phi = g \;  \text{on} \; \Gamma_1, \quad \text{and} \quad
 * -\text{div}\, \boldsymbol v_t \cdot \mathbf{n} = J \;  \text{on}  \; \Gamma_2. \f]
 *  Here, \f$ \Phi_\alpha \f$ denotes the potential of phase \f$ \alpha \f$, \f$ \boldsymbol K \f$ the intrinsic permeability,
 * \f$ \lambda_t \f$ the total mobility, \f$ f_\alpha \f$ the phase fractional flow function.
 *
 * More details on the equations can be found in
 *
 * M. Wolff, Y. Cao, B. Flemisch, R. Helmig, and B. Wohlmuth (2013a). Multi-point flux
 * approximation L-method in 3D: numerical convergence and application to two-phase
 * flow through porous media. In P. Bastian, J. Kraus, R. Scheichl, and M. Wheeler,
 * editors, Simulation of Flow in Porous Media - Applications in Energy and Environment. De Gruyter.
 *
 * M. Wolff, B. Flemisch, R. Helmig, I. Aavatsmark.
 * Treatment of tensorial relative permeabilities with multipoint flux approximation.
 * International Journal of Numerical Analysis and Modeling (9), pp. 725-744, 2012.
 *
 *  Remark1: only for 2-D quadrilateral grid
 *
 *  Remark2: implemented for UGGrid, ALUGrid, or YaspGrid
 *
 *\tparam TypeTag The problem Type Tag
 */
template<class TypeTag>
class FvMpfaO2dPressure2p: public FVPressure<TypeTag>
{
    using ParentType = FVPressure<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;
    using CellData = GetPropType<TypeTag, Properties::CellData>;
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using PrimaryVariables = typename SolutionTypes::PrimaryVariables;
    using ScalarSolutionType = typename SolutionTypes::ScalarSolution;

    using GridTypeIndices = GetPropType<TypeTag, Properties::GridTypeIndices>;

    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNw,
        sw = Indices::saturationW,
        sn = Indices::saturationNw
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        pressEqIdx = Indices::pressureEqIdx,
        satEqIdx = Indices::satEqIdx,
        numPhases = getPropValue<TypeTag, Properties::NumPhases>()
    };
    enum
    {
        globalCorner = 2,
        globalEdge = 3,
        neumannNeumann = 0,
        dirichletDirichlet = 1,
        dirichletNeumann = 2,
        neumannDirichlet = 3
    };

    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using IntersectionIterator = typename GridView::IntersectionIterator;
    using Intersection = typename GridView::Intersection;

    using LocalPosition = Dune::FieldVector<Scalar, dim>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GravityVector = Dune::FieldVector<Scalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

    using DimVector = Dune::FieldVector<Scalar, dim>;

    using InteractionVolume = FVMPFAOInteractionVolume<TypeTag>;

    using GlobalInteractionVolumeVector = std::vector<InteractionVolume>;
    using InnerBoundaryVolumeFaces = std::vector<Dune::FieldVector<bool, 2*dim> >;

    //! Helper function that finds the correct neighboring intersections
    Intersection getNextIntersection_(const Element&, const IntersectionIterator&);

    //! Initializes the matrix to store the system of equations
    friend class FVPressure<TypeTag>;
    void initializeMatrix();

    void storeInteractionVolumeInfo();

    //! Function which assembles the system of equations to be solved
    void assemble();
public:

    //! Constitutive functions are initialized and stored in the variables object
    void updateMaterialLaws();

    /*!
     * \brief Updates interaction volumes
     *
     * Globally rebuilds the MPFA interaction volumes.
     */
    void updateInteractionVolumeInfo()
    {
        interactionVolumes_.clear();
        innerBoundaryVolumeFaces_.clear();

        interactionVolumes_.resize(problem_.gridView().size(dim), InteractionVolume(problem_.gridView().grid()));
        innerBoundaryVolumeFaces_.resize(problem_.gridView().size(0), Dune::FieldVector<bool, 2 * dim>(false));

        storeInteractionVolumeInfo();
    }

    /*!
     * \brief Initializes the pressure model
     *
     * \copydetails FVPressure::initialize()
     */
    void initialize()
    {
        ParentType::initialize();

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

        interactionVolumes_.resize(problem_.gridView().size(dim), InteractionVolume(problem_.gridView().grid()));
        innerBoundaryVolumeFaces_.resize(problem_.gridView().size(0), Dune::FieldVector<bool, 2 * dim>(false));

        storeInteractionVolumeInfo();

        assemble();
        this->solve();

        storePressureSolution();

        return;
    }

    /*!
     * \brief Pressure update
     *
     * \copydetails FVPressure::update()
     */
    void update()
    {
        //error bounds for error term for incompressible models to correct unphysical saturation
        //over/undershoots due to saturation transport
        timeStep_ = problem_.timeManager().timeStepSize();
        maxError_ = 0.0;
        int size = problem_.gridView().size(0);
        for (int i = 0; i < size; i++)
        {
            Scalar sat = 0;
            using std::max;
            switch (saturationType_)
            {
            case sw:
                sat = problem_.variables().cellData(i).saturation(wPhaseIdx);
                break;
            case sn:
                sat = problem_.variables().cellData(i).saturation(nPhaseIdx);
                break;
            }
            if (sat > 1.0)
            {
                maxError_ = max(maxError_, (sat - 1.0) / timeStep_);
            }
            if (sat < 0.0)
            {
                maxError_ = max(maxError_, (-sat) / timeStep_);
            }
        }

        assemble();
        this->solve();

        storePressureSolution();

        return;
    }

    /*!
     * \brief Globally stores the pressure solution
     */
    void storePressureSolution()
    {
        // iterate through leaf grid an evaluate c0 at cell center
        for (const auto& element : elements(problem_.gridView()))
        {
            storePressureSolution(element);
        }
    }

    /*!
     * \brief Stores the pressure solution of a cell
     *
     * \param element Dune grid element
     */
    void storePressureSolution(const Element& element)
    {
        int eIdxGlobal = problem_.variables().index(element);
        CellData& cellData = problem_.variables().cellData(eIdxGlobal);

        switch (pressureType_)
        {
        case pw:
        {
            Scalar potW = this->pressure()[eIdxGlobal];

            Scalar gravityDiff = (problem_.bBoxMax() - element.geometry().center()) * gravity_;
            Scalar potPc = cellData.capillaryPressure()
                    + gravityDiff * (density_[nPhaseIdx] - density_[wPhaseIdx]);

            cellData.setPotential(wPhaseIdx, potW);
            cellData.setPotential(nPhaseIdx, potW + potPc);

            Scalar pressW = potW - gravityDiff * density_[wPhaseIdx];

            cellData.setPressure(wPhaseIdx, pressW);
            cellData.setPressure(nPhaseIdx, pressW + cellData.capillaryPressure());

            break;
        }
        case pn:
        {
            Scalar potNw = this->pressure()[eIdxGlobal];

            Scalar gravityDiff = (problem_.bBoxMax() - element.geometry().center()) * gravity_;
            Scalar potPc = cellData.capillaryPressure()
                    + gravityDiff * (density_[nPhaseIdx] - density_[wPhaseIdx]);

            cellData.setPotential(nPhaseIdx, potNw);
            cellData.setPotential(wPhaseIdx, potNw - potPc);

            Scalar pressNw = potNw - gravityDiff * density_[nPhaseIdx];

            cellData.setPressure(wPhaseIdx, pressNw - cellData.capillaryPressure());
            cellData.setPressure(nPhaseIdx, pressNw);

            break;
        }
        }

        cellData.fluxData().resetVelocity();
    }

    /*!
     * \brief Adds pressure output to the output file
     *
     * Adds the pressure, the potential and the capillary pressure to the output.
     * If the VtkOutputLevel is equal to zero (default) only primary variables are written,
     * if it is larger than zero also secondary variables are written.
     *
     * \tparam MultiWriter Class defining the output writer
     * \param writer The output writer (usually a <tt>VTKMultiWriter</tt> object)
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        int size = problem_.gridView().size(0);
        ScalarSolutionType *potential = writer.allocateManagedBuffer(size);

        (*potential) = this->pressure();

        if (pressureType_ == pw)
        {
            writer.attachCellData(*potential, "wetting potential");
        }

        if (pressureType_ == pn)
        {
            writer.attachCellData(*potential, "nonwetting potential");
        }

        if (vtkOutputLevel_ > 0)
        {
            ScalarSolutionType *pressure = writer.allocateManagedBuffer(size);
            ScalarSolutionType *pressureSecond = writer.allocateManagedBuffer(size);
            ScalarSolutionType *potentialSecond = writer.allocateManagedBuffer(size);
            ScalarSolutionType *pc = writer.allocateManagedBuffer(size);

            for (const auto& element : elements(problem_.gridView()))
            {
                int idx = problem_.variables().index(element);
                CellData& cellData = problem_.variables().cellData(idx);

                (*pc)[idx] = cellData.capillaryPressure();

                if (pressureType_ == pw)
                {
                    (*pressure)[idx] = cellData.pressure(wPhaseIdx);
                    (*potentialSecond)[idx] = cellData.potential(nPhaseIdx);
                    (*pressureSecond)[idx] = cellData.pressure(nPhaseIdx);
                }

                if (pressureType_ == pn)
                {
                    (*pressure)[idx] = cellData.pressure(nPhaseIdx);
                    (*potentialSecond)[idx] = cellData.potential(wPhaseIdx);
                    (*pressureSecond)[idx] = cellData.pressure(wPhaseIdx);
                }
            }

            if (pressureType_ == pw)
            {
                writer.attachCellData(*pressure, "wetting pressure");
                writer.attachCellData(*pressureSecond, "nonwetting pressure");
                writer.attachCellData(*potentialSecond, "nonwetting potential");
            }

            if (pressureType_ == pn)
            {
                writer.attachCellData(*pressure, "nonwetting pressure");
                writer.attachCellData(*pressureSecond, "wetting pressure");
                writer.attachCellData(*potentialSecond, "wetting potential");
            }

            writer.attachCellData(*pc, "capillary pressure");
        }

        return;
    }

    /*!
     * \brief Constructs a FvMpfaO2dPressure2p object
     * \param problem A problem class object
     */
    FvMpfaO2dPressure2p(Problem& problem) :
            ParentType(problem), problem_(problem), gravity_(problem.gravity()), maxError_(
                    0.), timeStep_(1.)
    {
        if (pressureType_ != pw && pressureType_ != pn)
        {
            DUNE_THROW(Dune::NotImplemented, "Pressure type not supported!");
        }
        if (saturationType_ != sw && saturationType_ != sn)
        {
            DUNE_THROW(Dune::NotImplemented, "Saturation type not supported!");
        }
        if (getPropValue<TypeTag, Properties::EnableCompressibility>())
        {
            DUNE_THROW(Dune::NotImplemented, "Compressibility not supported!");
        }
        if (dim != 2)
        {
            DUNE_THROW(Dune::NotImplemented, "Dimension not supported!");
        }

        ErrorTermFactor_ = getParam<Scalar>("Impet.ErrorTermFactor");
        ErrorTermLowerBound_ = getParam<Scalar>("Impet.ErrorTermLowerBound");
        ErrorTermUpperBound_ = getParam<Scalar>("Impet.ErrorTermUpperBound");

        density_[wPhaseIdx] = 0.;
        density_[nPhaseIdx] = 0.;
        viscosity_[wPhaseIdx] = 0.;
        viscosity_[nPhaseIdx] = 0.;

        vtkOutputLevel_ = getParam<int>("Vtk.OutputLevel");
    }

private:
    Problem& problem_;

    /*!
     * \brief Volume correction term
     *
     * Volume correction term to correct for unphysical saturation overshoots/undershoots.
     * These can occur if the estimated time step for the explicit transport was too large.
     * Correction by an artificial source term allows to correct this errors due to wrong
     * time-stepping without losing mass conservation. The error term looks as follows:
     * \f[
     *  q_{error} = \begin{cases}
     *          S < 0 & a_{error} \frac{S}{\Delta t} V \\
     *          S > 1 & a_{error} \frac{(S - 1)}{\Delta t} V \\
     *          0 \le S \le 1 & 0
     *      \end{cases}
     *  \f]
     *  where \f$a_{error}\f$ is a weighting factor (default: \f$a_{error} = 0.5\f$)
    */
    Scalar evaluateErrorTerm_(CellData& cellData)
    {
        //error term for incompressible models to correct unphysical saturation over/undershoots due to saturation transport
        // error reduction routine: volumetric error is damped and inserted to right hand side
        Scalar sat = 0;
        switch (saturationType_)
        {
        case sw:
            sat = cellData.saturation(wPhaseIdx);
            break;
        case sn:
            sat = cellData.saturation(nPhaseIdx);
            break;
        }

        Scalar error = (sat > 1.0) ? sat - 1.0 : 0.0;
        if (sat < 0.0)
        {
            error = sat;
        }
        error /= timeStep_;

        using std::abs;
        Scalar errorAbs = abs(error);

        if ((errorAbs * timeStep_ > 1e-6) && (errorAbs > ErrorTermLowerBound_ * maxError_)
                && (!problem_.timeManager().willBeFinished()))
        {
            return ErrorTermFactor_ * error;
        }
        return 0.0;
    }

protected:
    GlobalInteractionVolumeVector interactionVolumes_;
    InnerBoundaryVolumeFaces innerBoundaryVolumeFaces_;

private:
    const GravityVector& gravity_; //!< vector including the gravity constant

    Scalar maxError_;
    Scalar timeStep_;
    Scalar ErrorTermFactor_; //!< Handling of error term: relaxation factor
    Scalar ErrorTermLowerBound_; //!< Handling of error term: lower bound for error dampening
    Scalar ErrorTermUpperBound_; //!< Handling of error term: upper bound for error dampening

    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];

    int vtkOutputLevel_;

    static constexpr Scalar threshold_ = 1e-15;
    //! gives kind of pressure used (\f$ 0 = p_w\f$, \f$ 1 = p_n\f$, \f$ 2 = p_{global}\f$)
    static const int pressureType_ = getPropValue<TypeTag, Properties::PressureFormulation>();
    //! gives kind of saturation used (\f$ 0 = S_w\f$, \f$ 1 = S_n\f$)
    static const int saturationType_ = getPropValue<TypeTag, Properties::SaturationFormulation>();
    //! gives kind of velocity used (\f$ 0 = v_w\f$, \f$ 1 = v_n\f$, \f$ 2 = v_t\f$)
    static const int velocityType_ = getPropValue<TypeTag, Properties::VelocityFormulation>();
};

// TODO doc me!
template<class TypeTag>
typename FvMpfaO2dPressure2p<TypeTag>::Intersection
  FvMpfaO2dPressure2p<TypeTag>::getNextIntersection_(const Element& element,
                                                     const IntersectionIterator& isIt)
{
    auto isItBegin = problem_.gridView().ibegin(element);
    const auto isEndIt = problem_.gridView().iend(element);

    auto tempIsIt = isIt;
    auto nextIsIt = ++tempIsIt;

    // get 'nextIsIt'
    switch (getPropValue<TypeTag, Properties::GridImplementation>())
    {
        // for YaspGrid
        case GridTypeIndices::yaspGrid:
        {
            if (nextIsIt == isEndIt)
            {
                nextIsIt = isItBegin;
            }
            else
            {
                nextIsIt = ++tempIsIt;

                if (nextIsIt == isEndIt)
                {
                    auto tempIsItBegin = isItBegin;
                    nextIsIt = ++tempIsItBegin;
                }
            }

            break;
        }
        // for ALUGrid and UGGrid
        case GridTypeIndices::aluGrid:
        case GridTypeIndices::ugGrid:
        {
            if (nextIsIt == isEndIt)
                nextIsIt = isItBegin;

            break;
        }
        default:
        {
            DUNE_THROW(Dune::NotImplemented, "GridType can not be used with MPFAO implementation!");
            break;
        }
    }

    return *nextIsIt;
}

// TODO doc me!
template<class TypeTag>
void FvMpfaO2dPressure2p<TypeTag>::initializeMatrix()
{
    // determine matrix row sizes
    for (const auto& element : elements(problem_.gridView()))
    {
        // cell index
        int eIdxGlobalI = problem_.variables().index(element);

        // initialize row size
        int rowSize = 1;

        // run through all intersections with neighbors
        const auto isEndIt = problem_.gridView().iend(element);
        for (auto isIt = problem_.gridView().ibegin(element); isIt != isEndIt; ++isIt)
        {
            const auto& intersection = *isIt;
            auto nextIntersection = getNextIntersection_(element, isIt);

            if (intersection.neighbor())
                rowSize++;

            // try to find the common neighbor of intersection's and nextIntersection's outside
            if (intersection.neighbor() && nextIntersection.neighbor())
            {
                for (const auto& innerIntersection
                     : intersections(problem_.gridView(), intersection.outside()))
                    for (const auto& innerNextIntersection
                         : intersections(problem_.gridView(), nextIntersection.outside()))
                    {
                        if (innerIntersection.neighbor() && innerNextIntersection.neighbor())
                        {
                            if (innerIntersection.outside() == innerNextIntersection.outside()
                                && innerIntersection.outside() != intersection.inside())
                            {
                                rowSize++;
                            }
                        }
                    }
            }
        } // end of intersection loop

        // set number of indices in row eIdxGlobalI to rowSize
        this->A_.setrowsize(eIdxGlobalI, rowSize);

    } // end of element loop

    // indicate that size of all rows is defined
    this->A_.endrowsizes();

    // determine position of matrix entries
    for (const auto& element : elements(problem_.gridView()))
    {
        // cell index
        int eIdxGlobalI = problem_.variables().index(element);

        // add diagonal index
        this->A_.addindex(eIdxGlobalI, eIdxGlobalI);

        // run through all intersections with neighbors
        const auto isEndIt = problem_.gridView().iend(element);
        for (auto isIt = problem_.gridView().ibegin(element); isIt != isEndIt; ++isIt)
        {
            const auto& intersection = *isIt;
            auto nextIntersection = getNextIntersection_(element, isIt);

            if (intersection.neighbor())
            {
                // access neighbor
                int eIdxGlobalJ = problem_.variables().index(intersection.outside());

                // add off diagonal index
                // add index (row,col) to the matrix
                this->A_.addindex(eIdxGlobalI, eIdxGlobalJ);
            }

            if (intersection.neighbor() && nextIntersection.neighbor())
            {
                for (const auto& innerIntersection
                     : intersections(problem_.gridView(), intersection.outside()))
                    for (const auto& innerNextIntersection
                         : intersections(problem_.gridView(), nextIntersection.outside()))
                    {
                        if (innerIntersection.neighbor() && innerNextIntersection.neighbor())
                        {
                            auto innerOutside = innerIntersection.outside();

                            if (innerOutside == innerNextIntersection.outside()
                                && innerOutside != intersection.inside())
                            {
                                int eIdxGlobalJ = problem_.variables().index(innerOutside);

                                this->A_.addindex(eIdxGlobalI, eIdxGlobalJ);
                            }
                        }
                    }
            }
        } // end of intersection loop
    } // end of element loop

    // indicate that all indices are defined, check consistency
    this->A_.endindices();

    return;
}

//                 Indices used in a interaction volume of the MPFA-o method
//                 ___________________________________________________
//                 |                        |                        |
//                 | nuxy: cell geometry    |       nxy: face normal |
//                 |     vectors (see MPFA) |                        |
//                 |                        |                        |
//                 |            4-----------3-----------3            |
//                 |            | --> nu43  |  nu34 <-- |            |
//                 |            | |nu41    1|--> n43   ||nu32        |
//                 |            | v   ^     |0     ^   v|            |
//                 |____________4__0__|n14__|__n23_|_1__2____________|
//                 |            |    1    0 |     0     |
//                 |            | ^         |1   nu23 ^ |            |
//                 |            | |nu14    0|--> n12  | |            |
//                 |            | -->nu12   |   nu21<-- |            |
//                 |            1-----------1-----------2            |
//                 |          elementnumber |inter-                  |
//                 |                        |face-                   |
//                 |                        |number                  |
//                 |________________________|________________________|

// only for 2-D general quadrilateral
// TODO doc me!
template<class TypeTag>
void FvMpfaO2dPressure2p<TypeTag>::storeInteractionVolumeInfo()
{
    // introduce matrix R for vector rotation and R is initialized as zero matrix
    DimMatrix R(0);

    // evaluate matrix R
    if (dim == 2)
    {
        R[0][1] = 1;
        R[1][0] = -1;
    }

    BoundaryTypes bcType;

    // run through all elements
    for (const auto& element : elements(problem_.gridView()))
    {
        // get common geometry information for the following computation

        int eIdxGlobal1 = problem_.variables().index(element);
        // get global coordinate of cell 1 center
        const GlobalPosition& globalPos1 = element.geometry().center();

        // get absolute permeability of neighbor cell 1
        DimMatrix K1(problem_.spatialParams().intrinsicPermeability(element));

        const auto isIt12End = problem_.gridView().iend(element);
        for (auto isIt12 = problem_.gridView().ibegin(element); isIt12 != isIt12End; ++isIt12)
        {
            const auto& intersection12 = *isIt12;
            auto intersection14 = getNextIntersection_(element, isIt12);

            int indexInInside12 = intersection12.indexInInside();
            int indexInInside14 = intersection14.indexInInside();

            // get the intersection node /bar^{x_3} between 'isIt12' and 'isIt14', denoted as 'corner1234'
            // initialization of corner1234

            const auto refElement = referenceElement(element);

            GlobalPosition corner1234(0);

            int globalVertIdx1234 = 0;

            // get the global coordinate and global vertex index of corner1234
            for (int i = 0; i < intersection12.geometry().corners(); ++i)
            {
                bool finished = false;

                const GlobalPosition& isIt12corner = intersection12.geometry().corner(i);

                int localVertIdx12corner = refElement.subEntity(indexInInside12, dim - 1, i, dim);

                int globalVertIdx12corner = problem_.variables().index(element.template subEntity<dim>(localVertIdx12corner));

                for (int j = 0; j < intersection14.geometry().corners(); ++j)
                {
                    int localVertIdx14corner = refElement.subEntity(indexInInside14, dim - 1, j, dim);

                    int globalVertIdx14corner = problem_.variables().index(element.template subEntity<dim>(localVertIdx14corner));

                    if (globalVertIdx12corner == globalVertIdx14corner)
                    {
                        corner1234 = isIt12corner;

                        globalVertIdx1234 = globalVertIdx12corner;

                        finished = true;
                        break;
                    }
                }

                if (finished)
                {
                    break;
                }
            }

            if (interactionVolumes_[globalVertIdx1234].isStored())
            {
                continue;
            }
            else
            {
                interactionVolumes_[globalVertIdx1234].setStored();
                //                std::cout << "vertIdx = " << globalVertIdx1234 << "\n";
            }

            //store pointer 1
            interactionVolumes_[globalVertIdx1234].setSubVolumeElement(element, 0);
            interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection12.indexInInside(), 0, 0);
            interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection14.indexInInside(), 0, 1);

            // center of face in global coordinates, i.e., the midpoint of edge 'isIt12'
            const GlobalPosition& globalPosFace12 = intersection12.geometry().center();

            // get face volume
            Scalar faceVol12 = intersection12.geometry().volume() / 2.0;

            // get outer normal vector scaled with half volume of face 'isIt12'
            DimVector unitOuterNormal12 = intersection12.centerUnitOuterNormal();

            // center of face in global coordinates, i.e., the midpoint of edge 'isIt14'
            const GlobalPosition& globalPosFace41 = intersection14.geometry().center();

            // get face volume
            Scalar faceVol41 = intersection14.geometry().volume() / 2.0;

            // get outer normal vector scaled with half volume of face 'isIt14': for numbering of n see Aavatsmark, Eigestad
            DimVector unitOuterNormal14 = intersection14.centerUnitOuterNormal();

            // compute normal vectors nu14,nu12
            DimVector nu14(0);
            R.mv(globalPos1 - globalPosFace12, nu14);

            DimVector nu12(0);
            R.mv(globalPosFace41 - globalPos1, nu12);

            interactionVolumes_[globalVertIdx1234].setPermTimesNu(nu12, K1, 0, 0);
            interactionVolumes_[globalVertIdx1234].setPermTimesNu(nu14, K1, 0, 1);
            interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal12, 0, 0);
            interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal14, 0, 1);
            interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol12, 0, 0);
            interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol41, 0, 1);

            // compute dF1, the area of quadrilateral made by normal vectors 'nu'
            DimVector Rnu12(0);
            R.umv(nu12, Rnu12);
            using std::abs;
            interactionVolumes_[globalVertIdx1234].setDF(abs(nu14 * Rnu12), 0);

            // handle interior face
            if (intersection12.neighbor())
            {
                // access neighbor cell 2 of 'isIt12'
                auto element2 = intersection12.outside();

                int eIdxGlobal2 = problem_.variables().index(element2);

                //store pointer 2
                interactionVolumes_[globalVertIdx1234].setSubVolumeElement(element2, 1);
                interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection12.indexInOutside(), 1, 1);

                interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal12, 1, 1);
                interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol12, 1, 1);

                // get global coordinate of neighbor cell 2 center
                const GlobalPosition& globalPos2 = element2.geometry().center();

                // get absolute permeability of neighbor cell 2
                DimMatrix K2(problem_.spatialParams().intrinsicPermeability(element2));

                // 'isIt14' is an interior face
                if (intersection14.neighbor())
                {
                    // neighbor cell 3
                    // access neighbor cell 3
                    auto element4 = intersection14.outside();

                    //store pointer 4
                    interactionVolumes_[globalVertIdx1234].setSubVolumeElement(element4, 3);
                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection14.indexInOutside(), 3, 0);

                    interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal14, 3, 0);
                    interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol41, 3, 0);

                    // get basic information of cell 1,2's neighbor cell 3,4

                    // get global coordinate of neighbor cell 4 center
                    const GlobalPosition& globalPos4 = element4.geometry().center();

                    // get absolute permeability of neighbor cell 2
                    DimMatrix K4(problem_.spatialParams().intrinsicPermeability(element4));

                    // cell 3
                    GlobalPosition globalPos3(0);

                    GlobalPosition globalPosFace23(0);
                    GlobalPosition globalPosFace34(0);

                    for (const auto& intersection2
                         : intersections(problem_.gridView(), element2))
                    {
                        bool finished = false;

                        for (const auto& intersection4
                             : intersections(problem_.gridView(), element4))
                        {
                            if (intersection2.neighbor() && intersection4.neighbor())
                            {
                                auto element32 = intersection2.outside();
                                auto element34 = intersection4.outside();

                                // find the common neighbor cell between cell 2 and cell 3, except cell 1
                                if (element32 == element34 && element32 != element)
                                {
                                    //store pointer 3
                                    interactionVolumes_[globalVertIdx1234].setSubVolumeElement(element32, 2);

                                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection2.indexInInside(), 1,
                                            0);
                                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection2.indexInOutside(), 2,
                                            1);
                                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection4.indexInInside(), 3,
                                            1);
                                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection4.indexInOutside(), 2,
                                            0);

                                    // get global coordinate of neighbor cell 4 center
                                    globalPos3 = element32.geometry().center();

                                    globalPosFace23 = intersection2.geometry().center();
                                    globalPosFace34 = intersection4.geometry().center();

                                    Scalar faceVol23 = intersection2.geometry().volume() / 2.0;
                                    Scalar faceVol34 = intersection4.geometry().volume() / 2.0;

                                    // get outer normal vector scaled with half volume of face : for numbering of n see Aavatsmark, Eigestad
                                    DimVector unitOuterNormal23 = intersection2.centerUnitOuterNormal();

                                    DimVector unitOuterNormal43 = intersection4.centerUnitOuterNormal();

                                    interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal23, 1, 0);
                                    interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal23, 2, 1);
                                    interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal43, 2, 0);
                                    interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal43, 3, 1);
                                    interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol23, 1, 0);
                                    interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol23, 2, 1);
                                    interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol34, 2, 0);
                                    interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol34, 3, 1);

                                    // get absolute permeability of neighbor cell 2
                                    DimMatrix K3(
                                            problem_.spatialParams().intrinsicPermeability(element32));

                                    // compute normal vectors nu23, nu21; nu32, nu34; nu41, nu43;
                                    DimVector nu23(0);
                                    R.umv(globalPosFace12 - globalPos2, nu23);

                                    DimVector nu21(0);
                                    R.umv(globalPosFace23 - globalPos2, nu21);

                                    DimVector nu32(0);
                                    R.umv(globalPosFace34 - globalPos3, nu32);

                                    DimVector nu34(0);
                                    R.umv(globalPos3 - globalPosFace23, nu34);

                                    DimVector nu41(0);
                                    R.umv(globalPos4 - globalPosFace34, nu41);

                                    DimVector nu43(0);
                                    R.umv(globalPos4 - globalPosFace41, nu43);

                                    interactionVolumes_[globalVertIdx1234].setPermTimesNu(nu23, K2, 1, 0);
                                    interactionVolumes_[globalVertIdx1234].setPermTimesNu(nu21, K2, 1, 1);
                                    interactionVolumes_[globalVertIdx1234].setPermTimesNu(nu34, K3, 2, 0);
                                    interactionVolumes_[globalVertIdx1234].setPermTimesNu(nu32, K3, 2, 1);
                                    interactionVolumes_[globalVertIdx1234].setPermTimesNu(nu41, K4, 3, 0);
                                    interactionVolumes_[globalVertIdx1234].setPermTimesNu(nu43, K4, 3, 1);

                                    // compute dF2, dF3, dF4 i.e., the area of quadrilateral made by normal vectors 'nu'
                                    DimVector Rnu21(0);
                                    R.umv(nu21, Rnu21);
                                    interactionVolumes_[globalVertIdx1234].setDF(abs(nu23 * Rnu21), 1);

                                    DimVector Rnu34(0);
                                    R.umv(nu34, Rnu34);
                                    interactionVolumes_[globalVertIdx1234].setDF(abs(nu32 * Rnu34), 2);

                                    DimVector Rnu43(0);
                                    R.umv(nu43, Rnu43);
                                    interactionVolumes_[globalVertIdx1234].setDF(abs(nu41 * Rnu43), 3);

                                    finished = true;

                                    break;
                                }
                            }
                        }
                        if (finished)
                        {
                            break;
                        }
                    }
                }

                // 'isIt14' is on the boundary
                else
                {
                    problem_.boundaryTypes(bcType, intersection14);
                    PrimaryVariables boundValues(0.0);

                    interactionVolumes_[globalVertIdx1234].setBoundary(bcType, 3);
                    if (bcType.isNeumann(pressEqIdx))
                    {
                        problem_.neumann(boundValues, intersection14);
                        boundValues *= faceVol41;
                        interactionVolumes_[globalVertIdx1234].setNeumannCondition(boundValues, 3);
                    }
                    if (bcType.hasDirichlet())
                    {
                        problem_.dirichlet(boundValues, intersection14);
                        interactionVolumes_[globalVertIdx1234].setDirichletCondition(boundValues, 3);
                    }

                    // get common geometry information for the following computation
                    // get the information of the face 'isIt23' between cell2 and cell4 (locally numbered)

                    // center of face in global coordinates, i.e., the midpoint of edge 'isIt23'
                    GlobalPosition globalPosFace23(0);

                    // get face volume
                    Scalar faceVol23 = 0;

                    // get outer normal vector scaled with half volume of face 'isIt23'
                    DimVector unitOuterNormal23(0);

                    bool finished = false;

                    for (const auto& intersection2
                         : intersections(problem_.gridView(), element2))
                    {
                        if (intersection2.boundary())
                        {
                            for (int i = 0; i < intersection2.geometry().corners(); ++i)
                            {
                                int localVertIdx2corner = refElement.subEntity(intersection2.indexInInside(), dim - 1, i,
                                        dim);

                                int globalVertIdx2corner = problem_.variables().index(element2.template subEntity<dim>(localVertIdx2corner));

                                if (globalVertIdx2corner == globalVertIdx1234)
                                {
                                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection2.indexInInside(), 1,
                                            0);

                                    globalPosFace23 = intersection2.geometry().center();

                                    faceVol23 = intersection2.geometry().volume() / 2.0;

                                    unitOuterNormal23 = intersection2.centerUnitOuterNormal();

                                    interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal23, 1, 0);
                                    interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol23, 1, 0);

                                    problem_.boundaryTypes(bcType, intersection2);
                                    boundValues = 0.0;

                                    interactionVolumes_[globalVertIdx1234].setBoundary(bcType, 1);
                                    if (bcType.isNeumann(pressEqIdx))
                                    {
                                        problem_.neumann(boundValues, intersection2);
                                        boundValues *= faceVol23;
                                        interactionVolumes_[globalVertIdx1234].setNeumannCondition(boundValues, 1);
                                    }
                                    if (bcType.hasDirichlet())
                                    {
                                        problem_.dirichlet(boundValues, intersection2);
                                        interactionVolumes_[globalVertIdx1234].setDirichletCondition(boundValues, 1);
                                    }

                                    interactionVolumes_[globalVertIdx1234].setOutsideFace(2);

                                    innerBoundaryVolumeFaces_[eIdxGlobal1][intersection12.indexInInside()] = true;
                                    innerBoundaryVolumeFaces_[eIdxGlobal2][intersection12.indexInOutside()] = true;

                                    // compute normal vectors nu23, nu21;
                                    DimVector nu23(0);
                                    R.umv(globalPosFace12 - globalPos2, nu23);

                                    DimVector nu21(0);
                                    R.umv(globalPosFace23 - globalPos2, nu21);

                                    interactionVolumes_[globalVertIdx1234].setPermTimesNu(nu23, K2, 1, 0);
                                    interactionVolumes_[globalVertIdx1234].setPermTimesNu(nu21, K2, 1, 1);

                                    // compute dF2 i.e., the area of quadrilateral made by normal vectors 'nu'
                                    DimVector Rnu21(0);
                                    R.umv(nu21, Rnu21);
                                    interactionVolumes_[globalVertIdx1234].setDF(abs(nu23 * Rnu21), 1);

                                    finished = true;

                                    break;
                                }
                            }
                        }
                        if (finished)
                        {
                            break;
                        }
                    }
                    if (!finished)
                    {
                        DUNE_THROW(
                                Dune::NotImplemented,
                                "fvmpfao2pfaboundpressure2p.hh, l. 997: boundary shape not available as interaction volume shape");
                    }
                }
            }

            // handle boundary face 'isIt12'
            else
            {
                problem_.boundaryTypes(bcType, *isIt12);
                PrimaryVariables boundValues(0.0);

                interactionVolumes_[globalVertIdx1234].setBoundary(bcType, 0);
                if (bcType.isNeumann(pressEqIdx))
                {
                    problem_.neumann(boundValues, *isIt12);
                    boundValues *= faceVol12;
                    interactionVolumes_[globalVertIdx1234].setNeumannCondition(boundValues, 0);
                }
                if (bcType.hasDirichlet())
                {
                    problem_.dirichlet(boundValues, *isIt12);
                    interactionVolumes_[globalVertIdx1234].setDirichletCondition(boundValues, 0);
                }

                // 'isIt14' is on boundary
                if (intersection14.boundary())
                {
                    problem_.boundaryTypes(bcType, intersection14);
                    boundValues = 0.0;

                    interactionVolumes_[globalVertIdx1234].setBoundary(bcType, 3);
                    if (bcType.isNeumann(pressEqIdx))
                    {
                        problem_.neumann(boundValues, intersection14);
                        boundValues *= faceVol41;
                        interactionVolumes_[globalVertIdx1234].setNeumannCondition(boundValues, 3);
                    }
                    if (bcType.hasDirichlet())
                    {
                        problem_.dirichlet(boundValues, intersection14);
                        interactionVolumes_[globalVertIdx1234].setDirichletCondition(boundValues, 3);
                    }

                    interactionVolumes_[globalVertIdx1234].setOutsideFace(1);
                    interactionVolumes_[globalVertIdx1234].setOutsideFace(2);
                }

                // 'isIt14' is inside
                else
                {
                    // neighbor cell 3
                    // access neighbor cell 3
                    auto element4 = intersection14.outside();
                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection14.indexInOutside(), 3, 0);

                    //store pointer 4
                    interactionVolumes_[globalVertIdx1234].setSubVolumeElement(element4, 3);

                    interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal14, 3, 0);
                    interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol41, 3, 0);

                    // get global coordinate of neighbor cell 3 center
                    const GlobalPosition& globalPos4 = element4.geometry().center();

                    int eIdxGlobal4 = problem_.variables().index(element4);

                    bool finished = false;

                    // get the information of the face 'isIt34' between cell3 and cell4 (locally numbered)
                    for (const auto& intersection4
                         : intersections(problem_.gridView(), element4))
                    {
                        if (intersection4.boundary())
                        {
                            for (int i = 0; i < intersection4.geometry().corners(); ++i)
                            {
                                int localVertIdx4corner = refElement.subEntity(intersection4.indexInInside(), dim - 1, i,
                                        dim);

                                int globalVertIdx4corner = problem_.variables().index(element4.template subEntity<dim>(localVertIdx4corner));

                                if (globalVertIdx4corner == globalVertIdx1234)
                                {
                                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection4.indexInInside(), 3,
                                            1);

                                    const GlobalPosition& globalPosFace34 = intersection4.geometry().center();

                                    Scalar faceVol34 = intersection4.geometry().volume() / 2.0;

                                    DimVector unitOuterNormal43 = intersection4.centerUnitOuterNormal();

                                    interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal43, 3, 1);
                                    interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol34, 3, 1);

                                    problem_.boundaryTypes(bcType, intersection4);
                                    boundValues = 0.0;

                                    interactionVolumes_[globalVertIdx1234].setBoundary(bcType, 2);
                                    if (bcType.isNeumann(pressEqIdx))
                                    {
                                        problem_.neumann(boundValues, intersection4);
                                        boundValues *= faceVol34;
                                        interactionVolumes_[globalVertIdx1234].setNeumannCondition(boundValues, 2);
                                    }
                                    if (bcType.hasDirichlet())
                                    {
                                        problem_.dirichlet(boundValues, intersection4);
                                        interactionVolumes_[globalVertIdx1234].setDirichletCondition(boundValues, 2);
                                    }

                                    interactionVolumes_[globalVertIdx1234].setOutsideFace(1);

                                    innerBoundaryVolumeFaces_[eIdxGlobal1][intersection14.indexInInside()] = true;
                                    innerBoundaryVolumeFaces_[eIdxGlobal4][intersection14.indexInOutside()] = true;

                                    // get absolute permeability of neighbor cell 2
                                    DimMatrix K4(
                                            problem_.spatialParams().intrinsicPermeability(element4));

                                    // compute normal vectors nu41, nu43;
                                    DimVector nu41(0);
                                    R.umv(globalPos4 - globalPosFace34, nu41);

                                    DimVector nu43(0);
                                    R.umv(globalPos4 - globalPosFace41, nu43);

                                    interactionVolumes_[globalVertIdx1234].setPermTimesNu(nu41, K4, 3, 0);
                                    interactionVolumes_[globalVertIdx1234].setPermTimesNu(nu43, K4, 3, 1);

                                    // compute dF1, dF3 i.e., the area of quadrilateral made by normal vectors 'nu'
                                    DimVector Rnu43(0);
                                    R.umv(nu43, Rnu43);
                                    interactionVolumes_[globalVertIdx1234].setDF(abs(nu41 * Rnu43), 3);

                                    finished = true;

                                    break;
                                }
                            }
                        }
                        if (finished)
                        {
                            break;
                        }
                    }
                    if (!finished)
                    {
                        DUNE_THROW(
                                Dune::NotImplemented,
                                "fvmpfao2pfaboundpressure2p.hh, l. 1164: boundary shape not available as interaction volume shape");
                    }
                }
            }

        } // end all intersections
    } // end grid traversal

    return;
}

// only for 2-D general quadrilateral
//TODO doc me!
template<class TypeTag>
void FvMpfaO2dPressure2p<TypeTag>::assemble()
{
    // initialization: set global matrix this->A_ to zero
    this->A_ = 0;
    this->f_ = 0;

    // run through all vertices
    for (const auto& vertex : vertices(problem_.gridView()))
    {
        int vIdxGlobal = problem_.variables().index(vertex);

        InteractionVolume& interactionVolume = interactionVolumes_[vIdxGlobal];

        if (interactionVolume.isInnerVolume())
        {

            auto element1 = interactionVolume.getSubVolumeElement(0);
            auto element2 = interactionVolume.getSubVolumeElement(1);
            auto element3 = interactionVolume.getSubVolumeElement(2);
            auto element4 = interactionVolume.getSubVolumeElement(3);

            // get global coordinate of cell centers
            const GlobalPosition& globalPos1 = element1.geometry().center();
            const GlobalPosition& globalPos2 = element2.geometry().center();
            const GlobalPosition& globalPos3 = element3.geometry().center();
            const GlobalPosition& globalPos4 = element4.geometry().center();

            // cell volumes
            Scalar volume1 = element1.geometry().volume();
            Scalar volume2 = element2.geometry().volume();
            Scalar volume3 = element3.geometry().volume();
            Scalar volume4 = element4.geometry().volume();

            // cell index
            int eIdxGlobal1 = problem_.variables().index(element1);
            int eIdxGlobal2 = problem_.variables().index(element2);
            int eIdxGlobal3 = problem_.variables().index(element3);
            int eIdxGlobal4 = problem_.variables().index(element4);

            //get the cell Data
            CellData& cellData1 = problem_.variables().cellData(eIdxGlobal1);
            CellData& cellData2 = problem_.variables().cellData(eIdxGlobal2);
            CellData& cellData3 = problem_.variables().cellData(eIdxGlobal3);
            CellData& cellData4 = problem_.variables().cellData(eIdxGlobal4);

            // evaluate right hand side
            PrimaryVariables source(0.0);
            problem_.source(source, element1);
            this->f_[eIdxGlobal1] += volume1 / (4.0)
                    * (source[wPhaseIdx] / density_[wPhaseIdx] + source[nPhaseIdx] / density_[nPhaseIdx]);
            problem_.source(source, element2);
            this->f_[eIdxGlobal2] += volume2 / (4.0)
                    * (source[wPhaseIdx] / density_[wPhaseIdx] + source[nPhaseIdx] / density_[nPhaseIdx]);
            problem_.source(source, element3);
            this->f_[eIdxGlobal3] += volume3 / (4.0)
                    * (source[wPhaseIdx] / density_[wPhaseIdx] + source[nPhaseIdx] / density_[nPhaseIdx]);
            problem_.source(source, element4);
            this->f_[eIdxGlobal4] += volume4 / (4.0)
                    * (source[wPhaseIdx] / density_[wPhaseIdx] + source[nPhaseIdx] / density_[nPhaseIdx]);

            this->f_[eIdxGlobal1] += evaluateErrorTerm_(cellData1) * volume1 / (4.0);
            this->f_[eIdxGlobal2] += evaluateErrorTerm_(cellData2) * volume2 / (4.0);
            this->f_[eIdxGlobal3] += evaluateErrorTerm_(cellData3) * volume3 / (4.0);
            this->f_[eIdxGlobal4] += evaluateErrorTerm_(cellData4) * volume4 / (4.0);

            //get mobilities of the phases
            Dune::FieldVector<Scalar, numPhases> lambda1(cellData1.mobility(wPhaseIdx));
            lambda1[nPhaseIdx] = cellData1.mobility(nPhaseIdx);

            //compute total mobility of cell 1
            Scalar lambdaTotal1 = lambda1[wPhaseIdx] + lambda1[nPhaseIdx];

            //get mobilities of the phases
            Dune::FieldVector<Scalar, numPhases> lambda2(cellData2.mobility(wPhaseIdx));
            lambda2[nPhaseIdx] = cellData2.mobility(nPhaseIdx);

            //compute total mobility of cell 1
            Scalar lambdaTotal2 = lambda2[wPhaseIdx] + lambda2[nPhaseIdx];

            //get mobilities of the phases
            Dune::FieldVector<Scalar, numPhases> lambda3(cellData3.mobility(wPhaseIdx));
            lambda3[nPhaseIdx] = cellData3.mobility(nPhaseIdx);

            //compute total mobility of cell 1
            Scalar lambdaTotal3 = lambda3[wPhaseIdx] + lambda3[nPhaseIdx];

            //get mobilities of the phases
            Dune::FieldVector<Scalar, numPhases> lambda4(cellData4.mobility(wPhaseIdx));
            lambda4[nPhaseIdx] = cellData4.mobility(nPhaseIdx);

            //compute total mobility of cell 1
            Scalar lambdaTotal4 = lambda4[wPhaseIdx] + lambda4[nPhaseIdx];

            Scalar gn12nu14 = interactionVolume.getNtkrkNu_df(lambdaTotal1, 0, 0, 1);
            Scalar gn12nu12 = interactionVolume.getNtkrkNu_df(lambdaTotal1, 0, 0, 0);
            Scalar gn14nu14 = interactionVolume.getNtkrkNu_df(lambdaTotal1, 0, 1, 1);
            Scalar gn14nu12 = interactionVolume.getNtkrkNu_df(lambdaTotal1, 0, 1, 0);
            Scalar gn12nu23 = interactionVolume.getNtkrkNu_df(lambdaTotal2, 1, 1, 0);
            Scalar gn12nu21 = interactionVolume.getNtkrkNu_df(lambdaTotal2, 1, 1, 1);
            Scalar gn23nu23 = interactionVolume.getNtkrkNu_df(lambdaTotal2, 1, 0, 0);
            Scalar gn23nu21 = interactionVolume.getNtkrkNu_df(lambdaTotal2, 1, 0, 1);
            Scalar gn43nu32 = interactionVolume.getNtkrkNu_df(lambdaTotal3, 2, 0, 1);
            Scalar gn43nu34 = interactionVolume.getNtkrkNu_df(lambdaTotal3, 2, 0, 0);
            Scalar gn23nu32 = interactionVolume.getNtkrkNu_df(lambdaTotal3, 2, 1, 1);
            Scalar gn23nu34 = interactionVolume.getNtkrkNu_df(lambdaTotal3, 2, 1, 0);
            Scalar gn43nu41 = interactionVolume.getNtkrkNu_df(lambdaTotal4, 3, 1, 0);
            Scalar gn43nu43 = interactionVolume.getNtkrkNu_df(lambdaTotal4, 3, 1, 1);
            Scalar gn14nu41 = interactionVolume.getNtkrkNu_df(lambdaTotal4, 3, 0, 0);
            Scalar gn14nu43 = interactionVolume.getNtkrkNu_df(lambdaTotal4, 3, 0, 1);

            // compute transmissibility matrix T = CA^{-1}B+F
            Dune::FieldMatrix<Scalar, 2 * dim, 2 * dim> C(0), F(0), A(0), B(0);

            // evaluate matrix C, F, A, B
            C[0][0] = -gn12nu12;
            C[0][3] = -gn12nu14;
            C[1][0] = gn23nu21;
            C[1][1] = -gn23nu23;
            C[2][1] = gn43nu32;
            C[2][2] = gn43nu34;
            C[3][2] = -gn14nu43;
            C[3][3] = gn14nu41;

            F[0][0] = gn12nu12 + gn12nu14;
            F[1][1] = -gn23nu21 + gn23nu23;
            F[2][2] = -gn43nu34 - gn43nu32;
            F[3][3] = gn14nu43 - gn14nu41;

            A[0][0] = gn12nu12 + gn12nu21;
            A[0][1] = -gn12nu23;
            A[0][3] = gn12nu14;
            A[1][0] = -gn23nu21;
            A[1][1] = gn23nu23 + gn23nu32;
            A[1][2] = gn23nu34;
            A[2][1] = -gn43nu32;
            A[2][2] = -gn43nu34 - gn43nu43;
            A[2][3] = gn43nu41;
            A[3][0] = -gn14nu12;
            A[3][2] = gn14nu43;
            A[3][3] = -gn14nu41 - gn14nu14;

            B[0][0] = gn12nu12 + gn12nu14;
            B[0][1] = gn12nu21 - gn12nu23;
            B[1][1] = -gn23nu21 + gn23nu23;
            B[1][2] = gn23nu34 + gn23nu32;
            B[2][2] = -gn43nu34 - gn43nu32;
            B[2][3] = -gn43nu43 + gn43nu41;
            B[3][0] = -gn14nu12 - gn14nu14;
            B[3][3] = gn14nu43 - gn14nu41;

            // compute T
            A.invert();
            F += C.rightmultiply(B.leftmultiply(A));
            Dune::FieldMatrix<Scalar, 2 * dim, 2 * dim> T(F);

            //                        std::cout<<"Tpress = "<<T<<"\n";

            // assemble the global matrix this->A_ and right hand side f
            this->A_[eIdxGlobal1][eIdxGlobal1] += T[0][0] + T[3][0];
            this->A_[eIdxGlobal1][eIdxGlobal2] += T[0][1] + T[3][1];
            this->A_[eIdxGlobal1][eIdxGlobal3] += T[0][2] + T[3][2];
            this->A_[eIdxGlobal1][eIdxGlobal4] += T[0][3] + T[3][3];

            this->A_[eIdxGlobal2][eIdxGlobal1] += -T[0][0] + T[1][0];
            this->A_[eIdxGlobal2][eIdxGlobal2] += -T[0][1] + T[1][1];
            this->A_[eIdxGlobal2][eIdxGlobal3] += -T[0][2] + T[1][2];
            this->A_[eIdxGlobal2][eIdxGlobal4] += -T[0][3] + T[1][3];

            this->A_[eIdxGlobal3][eIdxGlobal1] -= T[1][0] + T[2][0];
            this->A_[eIdxGlobal3][eIdxGlobal2] -= T[1][1] + T[2][1];
            this->A_[eIdxGlobal3][eIdxGlobal3] -= T[1][2] + T[2][2];
            this->A_[eIdxGlobal3][eIdxGlobal4] -= T[1][3] + T[2][3];

            this->A_[eIdxGlobal4][eIdxGlobal1] += T[2][0] - T[3][0];
            this->A_[eIdxGlobal4][eIdxGlobal2] += T[2][1] - T[3][1];
            this->A_[eIdxGlobal4][eIdxGlobal3] += T[2][2] - T[3][2];
            this->A_[eIdxGlobal4][eIdxGlobal4] += T[2][3] - T[3][3];

            if (innerBoundaryVolumeFaces_[eIdxGlobal1][interactionVolume.getIndexOnElement(0, 0)])
            {
                this->A_[eIdxGlobal1][eIdxGlobal1] += T[0][0];
                this->A_[eIdxGlobal1][eIdxGlobal2] += T[0][1];
                this->A_[eIdxGlobal1][eIdxGlobal3] += T[0][2];
                this->A_[eIdxGlobal1][eIdxGlobal4] += T[0][3];
            }
            if (innerBoundaryVolumeFaces_[eIdxGlobal1][interactionVolume.getIndexOnElement(0, 1)])
            {
                this->A_[eIdxGlobal1][eIdxGlobal1] += T[3][0];
                this->A_[eIdxGlobal1][eIdxGlobal2] += T[3][1];
                this->A_[eIdxGlobal1][eIdxGlobal3] += T[3][2];
                this->A_[eIdxGlobal1][eIdxGlobal4] += T[3][3];

            }
            if (innerBoundaryVolumeFaces_[eIdxGlobal2][interactionVolume.getIndexOnElement(1, 0)])
            {
                this->A_[eIdxGlobal2][eIdxGlobal1] += T[1][0];
                this->A_[eIdxGlobal2][eIdxGlobal2] += T[1][1];
                this->A_[eIdxGlobal2][eIdxGlobal3] += T[1][2];
                this->A_[eIdxGlobal2][eIdxGlobal4] += T[1][3];
            }
            if (innerBoundaryVolumeFaces_[eIdxGlobal2][interactionVolume.getIndexOnElement(1, 1)])
            {
                this->A_[eIdxGlobal2][eIdxGlobal1] += -T[0][0];
                this->A_[eIdxGlobal2][eIdxGlobal2] += -T[0][1];
                this->A_[eIdxGlobal2][eIdxGlobal3] += -T[0][2];
                this->A_[eIdxGlobal2][eIdxGlobal4] += -T[0][3];
            }
            if (innerBoundaryVolumeFaces_[eIdxGlobal3][interactionVolume.getIndexOnElement(2, 0)])
            {
                this->A_[eIdxGlobal3][eIdxGlobal1] -= T[2][0];
                this->A_[eIdxGlobal3][eIdxGlobal2] -= T[2][1];
                this->A_[eIdxGlobal3][eIdxGlobal3] -= T[2][2];
                this->A_[eIdxGlobal3][eIdxGlobal4] -= T[2][3];
            }
            if (innerBoundaryVolumeFaces_[eIdxGlobal3][interactionVolume.getIndexOnElement(2, 1)])
            {
                this->A_[eIdxGlobal3][eIdxGlobal1] -= T[1][0];
                this->A_[eIdxGlobal3][eIdxGlobal2] -= T[1][1];
                this->A_[eIdxGlobal3][eIdxGlobal3] -= T[1][2];
                this->A_[eIdxGlobal3][eIdxGlobal4] -= T[1][3];
            }
            if (innerBoundaryVolumeFaces_[eIdxGlobal4][interactionVolume.getIndexOnElement(3, 0)])
            {
                this->A_[eIdxGlobal4][eIdxGlobal1] += -T[3][0];
                this->A_[eIdxGlobal4][eIdxGlobal2] += -T[3][1];
                this->A_[eIdxGlobal4][eIdxGlobal3] += -T[3][2];
                this->A_[eIdxGlobal4][eIdxGlobal4] += -T[3][3];
            }
            if (innerBoundaryVolumeFaces_[eIdxGlobal4][interactionVolume.getIndexOnElement(3, 1)])
            {
                this->A_[eIdxGlobal4][eIdxGlobal1] += T[2][0];
                this->A_[eIdxGlobal4][eIdxGlobal2] += T[2][1];
                this->A_[eIdxGlobal4][eIdxGlobal3] += T[2][2];
                this->A_[eIdxGlobal4][eIdxGlobal4] += T[2][3];
            }

            //add capillary pressure and gravity terms to right-hand-side
            //calculate capillary pressure velocity
            Dune::FieldVector<Scalar, 2 * dim> pc(0);
            pc[0] = cellData1.capillaryPressure();
            pc[1] = cellData2.capillaryPressure();
            pc[2] = cellData3.capillaryPressure();
            pc[3] = cellData4.capillaryPressure();

            Dune::FieldVector<Scalar, 2 * dim> gravityDiff(0);

//            std::cout<<"maxPos = "<<problem_.bBoxMax()<<"\n";

            gravityDiff[0] = (problem_.bBoxMax() - globalPos1) * gravity_ * (density_[nPhaseIdx] - density_[wPhaseIdx]);
            gravityDiff[1] = (problem_.bBoxMax() - globalPos2) * gravity_ * (density_[nPhaseIdx] - density_[wPhaseIdx]);
            gravityDiff[2] = (problem_.bBoxMax() - globalPos3) * gravity_ * (density_[nPhaseIdx] - density_[wPhaseIdx]);
            gravityDiff[3] = (problem_.bBoxMax() - globalPos4) * gravity_ * (density_[nPhaseIdx] - density_[wPhaseIdx]);

            pc += gravityDiff;

            if (pc[0] == 0 && pc[1] == 0 && pc[2] == 0 && pc[3] == 0)
            {
                continue;
            }

            Dune::FieldVector<Scalar, 2 * dim> pcFlux(0);

            T.mv(pc, pcFlux);

            //            std::cout<<"pcFlux = "<<pcFlux<<"\n";

            Scalar pcPotential12 = pcFlux[0];
            Scalar pcPotential14 = pcFlux[3];
            Scalar pcPotential32 = -pcFlux[1];
            Scalar pcPotential34 = -pcFlux[2];

            //            std::cout<<"pcPotential12 = "<<pcPotential12<<"\n";

            //compute mobilities of face 1
            Dune::FieldVector<Scalar, numPhases> lambda12Upw(0.0);
            lambda12Upw[wPhaseIdx] = (pcPotential12 >= 0) ? lambda1[wPhaseIdx] : lambda2[wPhaseIdx];
            lambda12Upw[nPhaseIdx] = (pcPotential12 >= 0) ? lambda1[nPhaseIdx] : lambda2[nPhaseIdx];

            //compute mobilities of face 4
            Dune::FieldVector<Scalar, numPhases> lambda14Upw(0.0);
            lambda14Upw[wPhaseIdx] = (pcPotential14 >= 0) ? lambda1[wPhaseIdx] : lambda4[wPhaseIdx];
            lambda14Upw[nPhaseIdx] = (pcPotential14 >= 0) ? lambda1[nPhaseIdx] : lambda4[nPhaseIdx];

            //compute mobilities of face 2
            Dune::FieldVector<Scalar, numPhases> lambda32Upw(0.0);
            lambda32Upw[wPhaseIdx] = (pcPotential32 >= 0) ? lambda3[wPhaseIdx] : lambda2[wPhaseIdx];
            lambda32Upw[nPhaseIdx] = (pcPotential32 >= 0) ? lambda3[nPhaseIdx] : lambda2[nPhaseIdx];

            //compute mobilities of face 3
            Dune::FieldVector<Scalar, numPhases> lambda34Upw(0.0);
            lambda34Upw[wPhaseIdx] = (pcPotential34 >= 0) ? lambda3[wPhaseIdx] : lambda4[wPhaseIdx];
            lambda34Upw[nPhaseIdx] = (pcPotential34 >= 0) ? lambda3[nPhaseIdx] : lambda4[nPhaseIdx];

            for (int i = 0; i < numPhases; i++)
            {
                Scalar lambdaT12 = lambda12Upw[wPhaseIdx] + lambda12Upw[nPhaseIdx];
                Scalar lambdaT14 = lambda14Upw[wPhaseIdx] + lambda14Upw[nPhaseIdx];
                Scalar lambdaT32 = lambda32Upw[wPhaseIdx] + lambda32Upw[nPhaseIdx];
                Scalar lambdaT34 = lambda34Upw[wPhaseIdx] + lambda34Upw[nPhaseIdx];
                Scalar fracFlow12 = (lambdaT12 > threshold_) ? lambda12Upw[i] / (lambdaT12) : 0.0;
                Scalar fracFlow14 = (lambdaT14 > threshold_) ? lambda14Upw[i] / (lambdaT14) : 0.0;
                Scalar fracFlow32 = (lambdaT32 > threshold_) ? lambda32Upw[i] / (lambdaT32) : 0.0;
                Scalar fracFlow34 = (lambdaT34 > threshold_) ? lambda34Upw[i] / (lambdaT34) : 0.0;

                Dune::FieldVector<Scalar, 2 * dim> pcFluxReal(pcFlux);

                pcFluxReal[0] *= fracFlow12;
                pcFluxReal[1] *= fracFlow32;
                pcFluxReal[2] *= fracFlow34;
                pcFluxReal[3] *= fracFlow14;

//                using std::isnan;
//                if (isnan(pcFluxReal.two_norm()))
//                                std::cout<<"pcFlux = "<<pcFlux<<"\n";

                switch (pressureType_)
                {
                case pw:
                {
                    if (i == nPhaseIdx)
                    {
                        //add capillary pressure term to right hand side
                        this->f_[eIdxGlobal1] -= (pcFluxReal[0] + pcFluxReal[3]);
                        this->f_[eIdxGlobal2] -= (pcFluxReal[1] - pcFluxReal[0]);
                        this->f_[eIdxGlobal3] -= (-pcFluxReal[2] - pcFluxReal[1]);
                        this->f_[eIdxGlobal4] -= (-pcFluxReal[3] + pcFluxReal[2]);

                        if (innerBoundaryVolumeFaces_[eIdxGlobal1][interactionVolume.getIndexOnElement(0, 0)])
                        {
                            this->f_[eIdxGlobal1] -= pcFluxReal[0];
                        }
                        if (innerBoundaryVolumeFaces_[eIdxGlobal1][interactionVolume.getIndexOnElement(0, 1)])
                        {
                            this->f_[eIdxGlobal1] -= pcFluxReal[3];
                        }
                        if (innerBoundaryVolumeFaces_[eIdxGlobal2][interactionVolume.getIndexOnElement(1, 0)])
                        {
                            this->f_[eIdxGlobal2] -= pcFluxReal[1];
                        }
                        if (innerBoundaryVolumeFaces_[eIdxGlobal2][interactionVolume.getIndexOnElement(1, 1)])
                        {
                            this->f_[eIdxGlobal2] += pcFluxReal[0];
                        }
                        if (innerBoundaryVolumeFaces_[eIdxGlobal3][interactionVolume.getIndexOnElement(2, 0)])
                        {
                            this->f_[eIdxGlobal3] += pcFluxReal[2];
                        }
                        if (innerBoundaryVolumeFaces_[eIdxGlobal3][interactionVolume.getIndexOnElement(2, 1)])
                        {
                            this->f_[eIdxGlobal3] += pcFluxReal[1];
                        }
                        if (innerBoundaryVolumeFaces_[eIdxGlobal4][interactionVolume.getIndexOnElement(3, 0)])
                        {
                            this->f_[eIdxGlobal4] += pcFluxReal[3];
                        }
                        if (innerBoundaryVolumeFaces_[eIdxGlobal4][interactionVolume.getIndexOnElement(3, 1)])
                        {
                            this->f_[eIdxGlobal4] -= pcFluxReal[2];
                        }
                    }
                    break;
                }
                case pn:
                {
                    if (i == wPhaseIdx)
                    {
                        //add capillary pressure term to right hand side
                        this->f_[eIdxGlobal1] += (pcFluxReal[0] + pcFluxReal[1]);
                        this->f_[eIdxGlobal2] += (pcFluxReal[1] - pcFluxReal[0]);
                        this->f_[eIdxGlobal3] += (-pcFluxReal[2] - pcFluxReal[1]);
                        this->f_[eIdxGlobal4] += (-pcFluxReal[3] + pcFluxReal[2]);

                        if (innerBoundaryVolumeFaces_[eIdxGlobal1][interactionVolume.getIndexOnElement(0, 0)])
                        {
                            this->f_[eIdxGlobal1] += pcFluxReal[0];
                        }
                        if (innerBoundaryVolumeFaces_[eIdxGlobal1][interactionVolume.getIndexOnElement(0, 1)])
                        {
                            this->f_[eIdxGlobal1] += pcFluxReal[3];
                        }
                        if (innerBoundaryVolumeFaces_[eIdxGlobal2][interactionVolume.getIndexOnElement(1, 0)])
                        {
                            this->f_[eIdxGlobal2] += pcFluxReal[1];
                        }
                        if (innerBoundaryVolumeFaces_[eIdxGlobal2][interactionVolume.getIndexOnElement(1, 1)])
                        {
                            this->f_[eIdxGlobal2] -= pcFluxReal[0];
                        }
                        if (innerBoundaryVolumeFaces_[eIdxGlobal3][interactionVolume.getIndexOnElement(2, 0)])
                        {
                            this->f_[eIdxGlobal3] -= pcFluxReal[2];
                        }
                        if (innerBoundaryVolumeFaces_[eIdxGlobal3][interactionVolume.getIndexOnElement(2, 1)])
                        {
                            this->f_[eIdxGlobal3] -= pcFluxReal[1];
                        }
                        if (innerBoundaryVolumeFaces_[eIdxGlobal4][interactionVolume.getIndexOnElement(3, 0)])
                        {
                            this->f_[eIdxGlobal4] -= pcFluxReal[3];
                        }
                        if (innerBoundaryVolumeFaces_[eIdxGlobal4][interactionVolume.getIndexOnElement(3, 1)])
                        {
                            this->f_[eIdxGlobal4] += pcFluxReal[2];
                        }
                    }
                    break;
                }
                }
            }
        }

        // at least one face on boundary!
        else
        {
            for (int elemIdx = 0; elemIdx < 2 * dim; elemIdx++)
            {
                bool isOutside = false;
                for (int fIdx = 0; fIdx < dim; fIdx++)
                {
                    int intVolFaceIdx = interactionVolume.getFaceIndexFromSubVolume(elemIdx, fIdx);
                    if (interactionVolume.isOutsideFace(intVolFaceIdx))
                    {
                        isOutside = true;
                        break;
                    }
                }
                if (isOutside)
                {
                    continue;
                }

                auto element = interactionVolume.getSubVolumeElement(elemIdx);

                // get global coordinate of cell centers
                const GlobalPosition& globalPos = element.geometry().center();

                // cell volumes
                Scalar volume = element.geometry().volume();

                // cell index
                int eIdxGlobal = problem_.variables().index(element);

                //get the cell Data
                CellData& cellData = problem_.variables().cellData(eIdxGlobal);

                //permeability vector at boundary
                DimMatrix permeability(problem_.spatialParams().intrinsicPermeability(element));

                // evaluate right hand side
                PrimaryVariables source(0);
                problem_.source(source, element);
                this->f_[eIdxGlobal] += volume / (4.0)
                        * (source[wPhaseIdx] / density_[wPhaseIdx] + source[nPhaseIdx] / density_[nPhaseIdx]);

                this->f_[eIdxGlobal] += evaluateErrorTerm_(cellData) * volume / (4.0);

                //get mobilities of the phases
                Dune::FieldVector<Scalar, numPhases> lambda(cellData.mobility(wPhaseIdx));
                lambda[nPhaseIdx] = cellData.mobility(nPhaseIdx);

                Scalar pc = cellData.capillaryPressure();

                Scalar gravityDiff = (problem_.bBoxMax() - globalPos) * gravity_
                        * (density_[nPhaseIdx] - density_[wPhaseIdx]);

                pc += gravityDiff; //minus because of gravity definition!

                for (int fIdx = 0; fIdx < dim; fIdx++)
                {
                    int intVolFaceIdx = interactionVolume.getFaceIndexFromSubVolume(elemIdx, fIdx);

                    if (interactionVolume.isBoundaryFace(intVolFaceIdx))
                    {

                        if (interactionVolume.getBoundaryType(intVolFaceIdx).isDirichlet(pressEqIdx))
                        {
                            int boundaryFaceIdx = interactionVolume.getIndexOnElement(elemIdx, fIdx);

                            const auto refElement = referenceElement(element);

                            const LocalPosition& localPos = refElement.position(boundaryFaceIdx, 1);

                            const GlobalPosition& globalPosFace = element.geometry().global(localPos);

                            DimVector distVec(globalPosFace - globalPos);
                            Scalar dist = distVec.two_norm();
                            DimVector unitDistVec(distVec);
                            unitDistVec /= dist;

                            Scalar faceArea = interactionVolume.getFaceArea(elemIdx, fIdx);

                            // get pc and lambda at the boundary
                            Scalar satWBound = cellData.saturation(wPhaseIdx);
                            //check boundary sat at face 1
                            if (interactionVolume.getBoundaryType(intVolFaceIdx).isDirichlet(satEqIdx))
                            {
                                Scalar satBound = interactionVolume.getDirichletValues(intVolFaceIdx)[saturationIdx];
                                switch (saturationType_)
                                {
                                case sw:
                                {
                                    satWBound = satBound;
                                    break;
                                }
                                case sn:
                                {
                                    satWBound = 1 - satBound;
                                    break;
                                }
                                }

                            }

                            // old material law interface is deprecated: Replace this by
                            // const auto& fluidMatrixInteraction = spatialParams.fluidMatrixInteractionAtPos(element.geometry().center());
                            // after the release of 3.3, when the deprecated interface is no longer supported
                            const auto fluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem_.spatialParams(), element);

                            Scalar pcBound = fluidMatrixInteraction.pc(satWBound);

                            Scalar gravityDiffBound = (problem_.bBoxMax() - globalPosFace) * gravity_
                                    * (density_[nPhaseIdx] - density_[wPhaseIdx]);

                            pcBound += gravityDiffBound;

                            Dune::FieldVector<Scalar, numPhases> lambdaBound(fluidMatrixInteraction.krw(satWBound));
                            lambdaBound[nPhaseIdx] = fluidMatrixInteraction.krn(satWBound);
                            lambdaBound[wPhaseIdx] /= viscosity_[wPhaseIdx];
                            lambdaBound[nPhaseIdx] /= viscosity_[nPhaseIdx];

                            Scalar potentialBound = interactionVolume.getDirichletValues(intVolFaceIdx)[pressureIdx];
                            Scalar gdeltaZ = (problem_.bBoxMax()-globalPosFace) * gravity_;

                            //calculate potential gradients
                            Scalar potentialDiffW = 0;
                            Scalar potentialDiffNw = 0;
                            switch (pressureType_)
                            {
                            case pw:
                            {
                                potentialBound += density_[wPhaseIdx]*gdeltaZ;
                                potentialDiffW = (cellData.potential(wPhaseIdx) - potentialBound) / dist;
                                potentialDiffNw = (cellData.potential(nPhaseIdx) - potentialBound - pcBound) / dist;
                                break;
                            }
                            case pn:
                            {
                                potentialBound += density_[nPhaseIdx]*gdeltaZ;
                                potentialDiffW = (cellData.potential(wPhaseIdx) - potentialBound + pcBound) / dist;
                                potentialDiffNw = (cellData.potential(nPhaseIdx) - potentialBound) / dist;
                                break;
                            }
                            }

                            Scalar lambdaTotal = (potentialDiffW >= 0.) ? lambda[wPhaseIdx] : lambdaBound[wPhaseIdx];
                            lambdaTotal += (potentialDiffNw >= 0.) ? lambda[nPhaseIdx] : lambdaBound[nPhaseIdx];

                            DimVector permTimesNormal(0);
                            permeability.mv(unitDistVec, permTimesNormal);

                            //calculate current matrix entry
                            Scalar entry = lambdaTotal * (unitDistVec * permTimesNormal) / dist * faceArea;

                            //calculate right hand side
                            Scalar pcFlux = 0;

                            switch (pressureType_)
                            {
                            case pw:
                            {
                                // calculate capillary pressure gradient
                                DimVector pcGradient = unitDistVec;
                                pcGradient *= (pc - pcBound) / dist;

                                //add capillary pressure term to right hand side
                                pcFlux = 0.5 * (lambda[nPhaseIdx] + lambdaBound[nPhaseIdx])
                                        * (permTimesNormal * pcGradient) * faceArea;

                                break;
                            }
                            case pn:
                            {
                                // calculate capillary pressure gradient
                                DimVector pcGradient = unitDistVec;
                                pcGradient *= (pc - pcBound) / dist;

                                //add capillary pressure term to right hand side
                                pcFlux = 0.5 * (lambda[wPhaseIdx] + lambdaBound[wPhaseIdx])
                                        * (permTimesNormal * pcGradient) * faceArea;

                                break;

                            }
                            }

                            // set diagonal entry and right hand side entry
                            this->A_[eIdxGlobal][eIdxGlobal] += entry;
                            this->f_[eIdxGlobal] += entry * potentialBound;

                            if (pc == 0 && pcBound == 0)
                            {
                                continue;
                            }

                            for (int i = 0; i < numPhases; i++)
                            {
                                switch (pressureType_)
                                {
                                case pw:
                                {
                                    if (i == nPhaseIdx)
                                    {
                                        //add capillary pressure term to right hand side
                                        this->f_[eIdxGlobal] -= pcFlux;
                                    }
                                    break;
                                }
                                case pn:
                                {
                                    if (i == wPhaseIdx)
                                    {
                                        //add capillary pressure term to right hand side
                                        this->f_[eIdxGlobal] += pcFlux;
                                    }
                                    break;
                                }
                                }
                            }

                        }
                        else if (interactionVolume.getBoundaryType(intVolFaceIdx).isNeumann(pressEqIdx))
                        {
                            Scalar J = interactionVolume.getNeumannValues(intVolFaceIdx)[wPhaseIdx]
                                    / density_[wPhaseIdx];
                            J += interactionVolume.getNeumannValues(intVolFaceIdx)[nPhaseIdx] / density_[nPhaseIdx];
                            this->f_[eIdxGlobal] -= J;
                        }
                        else
                        {
                            std::cout << "interactionVolume.getBoundaryType(intVolFaceIdx).isNeumann(pressEqIdx)"
                                    << interactionVolume.getBoundaryType(intVolFaceIdx).isNeumann(pressEqIdx) << "\n";
                            DUNE_THROW(Dune::NotImplemented,
                                    "No valid boundary condition type defined for pressure equation!");
                        }
                    }
                }
            }

        } // end boundaries

    } // end vertex iterator

    // only do more if we have more than one process
    if (problem_.gridView().comm().size() > 1)
    {
        // set ghost and overlap element entries
        for (const auto& element : elements(problem_.gridView()))
        {
            if (element.partitionType() == Dune::InteriorEntity)
                continue;

            // get the global index of the cell
            int eIdxGlobalI = problem_.variables().index(element);

            this->A_[eIdxGlobalI] = 0.0;
            this->A_[eIdxGlobalI][eIdxGlobalI] = 1.0;
            this->f_[eIdxGlobalI] = this->pressure()[eIdxGlobalI];
        }
    }

    return;
}


/*!
 * \brief Updates constitutive relations and stores them in the variable class
 *
 * Stores mobility, fractional flow function and capillary pressure for all grid cells.
 */
template<class TypeTag>
void FvMpfaO2dPressure2p<TypeTag>::updateMaterialLaws()
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

        const Scalar pc = fluidMatrixInteraction.pc(satW);

        cellData.setCapillaryPressure(pc);

        // initialize mobilities
        const Scalar mobilityW = fluidMatrixInteraction.krw(satW) / viscosity_[wPhaseIdx];
        const Scalar mobilityNw = fluidMatrixInteraction.krn(satW) / viscosity_[nPhaseIdx];

        // initialize mobilities
        cellData.setMobility(wPhaseIdx, mobilityW);
        cellData.setMobility(nPhaseIdx, mobilityNw);

        //initialize fractional flow functions
        cellData.setFracFlowFunc(wPhaseIdx, mobilityW / (mobilityW + mobilityNw));
        cellData.setFracFlowFunc(nPhaseIdx, mobilityNw / (mobilityW + mobilityNw));
    }
    return;
}

} // end namespace Dumux
#endif
