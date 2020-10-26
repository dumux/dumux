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
 * \brief Grid adaptive finite volume MPFA L-method discretization of a two-phase
 * pressure equation of the sequential IMPES model.
 */
#ifndef DUMUX_FVMPFAL2DPRESSURE2P_ADAPTIVE_HH
#define DUMUX_FVMPFAL2DPRESSURE2P_ADAPTIVE_HH

// dumux environment
#include <dumux/porousmediumflow/sequential/cellcentered/pressure.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/mpfa/linteractionvolume.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/properties.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/mpfa/properties.hh>
#include "2dtransmissibilitycalculator.hh"

#include <dumux/common/deprecated.hh>

namespace Dumux {

/*!
 * \ingroup SequentialTwoPModel
 * \brief Grid adaptive finite volume MPFA L-method discretization of a two-phase flow pressure equation of the sequential IMPES model.
 *
 * Grid adaptive finite volume MPFA L-method discretization of the equations
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
 * Wolff 2013: http://elib.uni-stuttgart.de/opus/volltexte/2013/8661/
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
 *
 *  Remark1: only for 2-D quadrilateral grid
 *
 *  Remark2: implemented for UGGrid, ALUGrid
 *
 *  Remark3: Allowed difference in grid levels of two neighboring cells: 1
 *
 *\tparam TypeTag The problem Type Tag
 */
template<class TypeTag>
class FvMpfaL2dPressure2pAdaptive: public FVPressure<TypeTag>
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
    using Grid = typename GridView::Grid;

    using LocalPosition = Dune::FieldVector<Scalar, dim>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

    using DimVector = Dune::FieldVector<Scalar, dim>;

public:


    /*!
     * \brief Type of the interaction volume objects
     *
     * Type of the interaction volume objects used to store the geometric information which is needed
     * to calculated the transmissibility matrices of one MPFA interaction volume.
     */
    using InteractionVolume = FVMPFALInteractionVolume<TypeTag>;
    using TransmissibilityCalculator = FvMpfaL2dTransmissibilityCalculator<TypeTag>;
private:

    using GlobalInteractionVolumeVector = std::vector<InteractionVolume>;
    using InnerBoundaryVolumeFaces = std::vector<Dune::FieldVector<bool, 2*dim> >;

    //! helper function thats find the correct neighboring intersections
    Intersection getNextIntersection_(const Element&, const IntersectionIterator&);

    //! initializes the matrix to store the system of equations
    friend class FVPressure<TypeTag>;
    void initializeMatrix();

    void storeInteractionVolumeInfo();

    void printInteractionVolumes();

    //! function which assembles the system of equations to be solved
    void assemble();

public:

    //! constitutive functions are initialized and stored in the variables object
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
//        printInteractionVolumes();
    }

    /*!
     * \brief Initializes the pressure model
     *
     * \copydetails ParentType::initialize()
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

        updateInteractionVolumeInfo();

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
     * \brief Pressure update
     *
     * \copydetails ParentType::update()
     */
    void update()
    {
        int gridSize = problem_.gridView().size(0);

        //error bounds for error term for incompressible models to correct unphysical saturation
        //over/undershoots due to saturation transport
        timeStep_ = problem_.timeManager().timeStepSize();
        maxError_ = 0.0;
        for (int i = 0; i < gridSize; i++)
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

        if (problem_.gridAdapt().wasAdapted())
        {
            // update RHS vector, matrix
            this->A_.setSize(gridSize, gridSize); //
            this->f_.resize(gridSize);
            this->pressure().resize(gridSize);

            for (int i = 0; i < gridSize; i++)
            {
                CellData& cellData = problem_.variables().cellData(i);

                switch (pressureType_)
                {
                case pw:
                    this->pressure()[i] = cellData.pressure(wPhaseIdx);
                    break;
                case pn:
                    this->pressure()[i] = cellData.pressure(nPhaseIdx);
                    break;
                }
            }

            initializeMatrix();
            updateInteractionVolumeInfo();
        }
        assemble();
        this->solve();


        storePressureSolution();

        return;
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
     * \brief Constructs a FvMpfaL2dPressure2pAdaptive object
     * \param problem A problem class object
     */
    FvMpfaL2dPressure2pAdaptive(Problem& problem) :
        ParentType(problem), problem_(problem), transmissibilityCalculator_(problem),
                    gravity_(problem.gravity()),
                    maxError_(0.), timeStep_(1.)
    {
        if (pressureType_ != pw && pressureType_ != pn)
        {
            DUNE_THROW(Dune::NotImplemented, "Pressure type not supported!");
        }
        if (saturationType_ != sw && saturationType_ != sn)
        {
            DUNE_THROW(Dune::NotImplemented, "Saturation type not supported!");
        }
        // The methods to calculate transmissibility are also used by the compressible
        // compositional models (2p2c), although this implementation does not support it.
        // Hence warning is not thrown if the file is not used from compositional modules
        // that require 2p2cproperties.hh.
        #ifndef DUMUX_2P2CPROPERTIES_HH
        if (getPropValue<TypeTag, Properties::EnableCompressibility>())
        {
            DUNE_THROW(Dune::NotImplemented, "Compressibility not supported!");
        }
        #endif
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
    TransmissibilityCalculator transmissibilityCalculator_;

protected:
    GlobalInteractionVolumeVector interactionVolumes_;//!< Global Vector of interaction volumes
    InnerBoundaryVolumeFaces innerBoundaryVolumeFaces_;//!< Vector marking faces which intersect the boundary

private:
    const Dune::FieldVector<Scalar, dimWorld>& gravity_; //!< vector including the gravity constant

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

    // TODO doc me!
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

};

// TODO doc me!
template<class TypeTag>
typename FvMpfaL2dPressure2pAdaptive<TypeTag>::Intersection
  FvMpfaL2dPressure2pAdaptive<TypeTag>::getNextIntersection_(const Element& element,
                                                             const IntersectionIterator& isIt)
{
    auto isItBegin = problem_.gridView().ibegin(element);
    const auto isEndIt = problem_.gridView().iend(element);

    auto tempIsIt = isIt;
    auto nextIsIt = ++tempIsIt;

    // get 'nextIsIt'
    switch (getPropValue<TypeTag, Properties::GridImplementation>())
    {
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
            DUNE_THROW(Dune::NotImplemented,
                       "GridType can not be used with adaptive MPFAL!");
            break;
        }
    }

    return *nextIsIt;
}

// TODO doc me!
template<class TypeTag>
void FvMpfaL2dPressure2pAdaptive<TypeTag>::initializeMatrix()
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

            if (intersection.neighbor())
            {
                rowSize++;

                auto nextIntersection = getNextIntersection_(element, isIt);

                if (nextIntersection.neighbor())
                {
                    bool isCorner = true;

                    for (const auto& innerIntersection
                         : intersections(problem_.gridView(), intersection.outside()))
                    {
                        for (const auto& innerNextIntersection
                             : intersections(problem_.gridView(), nextIntersection.outside()))
                        {
                            if (innerIntersection.neighbor() && innerNextIntersection.neighbor())
                            {
                                if (innerIntersection.outside() == nextIntersection.outside()
                                    || innerNextIntersection.outside() == intersection.outside())
                                {
                                    isCorner = false;
                                    break;
                                }
                            }
                        }
                        if (!isCorner)
                        {
                            break;
                        }
                    }
                    if (isCorner)
                    {
                        rowSize++;
                    }
                }
            }

        } // end of intersection loop

        using std::min;
        rowSize = min(rowSize, 13); //in 2-D

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

            if (intersection.neighbor())
            {
                // access neighbor
                auto outside = intersection.outside();
                int eIdxGlobalJ = problem_.variables().index(outside);

                // add off diagonal index
                // add index (row,col) to the matrix
                this->A_.addindex(eIdxGlobalI, eIdxGlobalJ);

                if (element.level() < outside.level())
                {
                    continue;
                }

                auto nextIntersection = getNextIntersection_(element, isIt);

                if (nextIntersection.neighbor())
                {
                    // access the common neighbor of intersection's and nextIntersection's outside
                    if (element.level() < nextIntersection.outside().level())
                    {
                        continue;
                    }

                    for (const auto& innerIntersection
                         : intersections(problem_.gridView(), intersection.outside()))
                    {
                        for (const auto& innerNextIntersection
                             : intersections(problem_.gridView(), nextIntersection.outside()))
                        {
                            if (innerIntersection.neighbor() && innerNextIntersection.neighbor())
                            {
                                auto innerOutside = innerIntersection.outside();
                                auto nextOutside = nextIntersection.outside();

                                if (innerOutside == innerNextIntersection.outside() && innerOutside != element
                                        && innerOutside != nextOutside)
                                {
                                    int eIdxGlobalCorner = problem_.variables().index(innerOutside);

                                    this->A_.addindex(eIdxGlobalI, eIdxGlobalCorner);

                                    if (element.level() > outside.level())
                                    {
                                        int eIdxGlobalJCorner = problem_.variables().index(nextOutside);

                                        this->A_.addindex(eIdxGlobalJ, eIdxGlobalJCorner);
                                    }
                                    if (element.level() > nextOutside.level())
                                    {
                                        int eIdxGlobalJCorner = problem_.variables().index(nextOutside);

                                        this->A_.addindex(eIdxGlobalJCorner, eIdxGlobalJ);
                                    }
                                    if (element.level() > innerOutside.level())
                                    {
                                        this->A_.addindex(eIdxGlobalCorner, eIdxGlobalI);
                                    }
                                }
                            }
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
void FvMpfaL2dPressure2pAdaptive<TypeTag>::storeInteractionVolumeInfo()
{
    BoundaryTypes bcType;

    // run through all elements
    for (const auto& element : elements(problem_.gridView()))
    {
        // get index
        int eIdxGlobal1 = problem_.variables().index(element);

        const auto refElement = referenceElement(element);

        const auto isEndIt12 = problem_.gridView().iend(element);
        for (auto isIt12 = problem_.gridView().ibegin(element); isIt12 != isEndIt12; ++isIt12)
        {
            const auto& intersection12 = *isIt12;

            int indexInInside12 = intersection12.indexInInside();

            if (intersection12.neighbor())
            {
                // access neighbor cell 2 of 'intersection12'
                auto element2 = intersection12.outside();
                int eIdxGlobal2 = problem_.variables().index(element2);

                if (element.level() < element2.level())
                {
                    continue;
                }

                auto intersection14 = getNextIntersection_(element, isIt12);

                int indexInInside14 = intersection14.indexInInside();

                //get center vertex
                GlobalPosition corner1234(0);
                int globalVertIdx1234 = 0;
                bool finished = false;
                // get the global coordinate and global vertex index of corner1234
                for (int i = 0; i < intersection12.geometry().corners(); ++i)
                {
                    int localVertIdx12corner = refElement.subEntity(indexInInside12, dim - 1, i, dim);

                    int globalVertIdx12corner = problem_.variables().index(element.template subEntity<dim>(localVertIdx12corner));

                    for (int j = 0; j < intersection14.geometry().corners(); ++j)
                    {
                        int localVertIdx14corner = refElement.subEntity(indexInInside14, dim - 1, j, dim);

                        int globalVertIdx14corner = problem_.variables().index(element.template subEntity<dim>(localVertIdx14corner));

                        if (globalVertIdx12corner == globalVertIdx14corner)
                        {
                            corner1234 = element.geometry().corner(localVertIdx12corner);

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

                if (interactionVolumes_[globalVertIdx1234].isStored() || !finished)
                {
                    continue;
                }
                else
                {
                    interactionVolumes_[globalVertIdx1234].setStored();
                }

                interactionVolumes_[globalVertIdx1234].setCenterPosition(corner1234);

                //store pointer 1
                interactionVolumes_[globalVertIdx1234].setSubVolumeElement(element, 0);
                interactionVolumes_[globalVertIdx1234].setIndexOnElement(indexInInside12, 0, 0);
                interactionVolumes_[globalVertIdx1234].setIndexOnElement(indexInInside14, 0, 1);

                // center of face in global coordinates, i.e., the midpoint of edge 'intersection12'
                const GlobalPosition& globalPosFace12 = intersection12.geometry().center();

                // get face volume
                Scalar faceVol12 = intersection12.geometry().volume() / 2.0;

                // get outer normal vector scaled with half volume of face 'intersection12'
                DimVector unitOuterNormal12 = intersection12.centerUnitOuterNormal();

                // center of face in global coordinates, i.e., the midpoint of edge 'intersection14'
                GlobalPosition globalPosFace41 = intersection14.geometry().center();

                // get face volume
                Scalar faceVol41 = intersection14.geometry().volume() / 2.0;

                // get outer normal vector scaled with half volume of face 'intersection14': for numbering of n see Aavatsmark, Eigestad
                DimVector unitOuterNormal14 = intersection14.centerUnitOuterNormal();

                interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal12, 0, 0);
                interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal14, 0, 1);
                //get the normals of from cell 2 and 4
                unitOuterNormal14 *= -1;
                unitOuterNormal12 *= -1;
                interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol12, 0, 0);
                interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol41, 0, 1);
                interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace12, 0, 0);
                interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace41, 0, 1);

                //store pointer 2
                interactionVolumes_[globalVertIdx1234].setSubVolumeElement(element2, 1);
                interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection12.indexInOutside(), 1, 1);
                interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal12, 1, 1);
                interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol12, 1, 1);
                interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace12, 1, 1);

                // 'intersection14' is an interior face
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
                    interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace41, 3, 0);

                    // cell 3
                    GlobalPosition globalPosFace23(0);
                    GlobalPosition globalPosFace34(0);

                    if (element4.level() < element.level())
                    {
                        bool isHangingNode = false;

                        for (const auto& intersection2
                             : intersections(problem_.gridView(), element2))
                        {
                            bool breakLoop = false;
                            for (const auto& intersection4
                             : intersections(problem_.gridView(), element4))
                            {
                                if (intersection2.neighbor() && intersection4.neighbor())
                                {
                                    auto element32 = intersection2.outside();
                                    auto element34 = intersection4.outside();

                                    //hanging node!
                                    if (element32 == element4)
                                    {
                                        if (element.level() != element2.level())
                                        {
                                            breakLoop = true;
                                            isHangingNode = false;
                                            break;
                                        }

                                        isHangingNode = true;

                                        interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection2.indexInInside(),
                                                1, 0);
                                        interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection2.indexInOutside(),
                                                3, 1);

                                        globalPosFace23 = intersection2.geometry().center();
                                        interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace23, 1, 0);
                                        interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace23, 3, 1);

                                        Scalar faceVol23 = intersection2.geometry().volume() / 2.0;
                                        interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol23, 1, 0);
                                        interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol23, 3, 1);

                                        // get outer normal vector scaled with half volume of face : for numbering of n see Aavatsmark, Eigestad
                                        DimVector unitOuterNormal23 = intersection2.centerUnitOuterNormal();

                                        interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal23, 1, 0);
                                        unitOuterNormal23 *= -1;
                                        interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal23, 3, 1);
                                    }
                                    else if (element34 == element2)
                                    {
                                        if (element.level() != element2.level())
                                        {
                                            breakLoop = true;
                                            isHangingNode = false;
                                            break;
                                        }

                                        isHangingNode = true;

                                        interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection4.indexInInside(),
                                                3, 1);
                                        interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection4.indexInOutside(),
                                                1, 0);

                                        // get outer normal vector scaled with half volume of face : for numbering of n see Aavatsmark, Eigestad
                                        DimVector unitOuterNormal43 = intersection4.centerUnitOuterNormal();
                                        interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal43, 3, 1);
                                    }
                                }
                            }
                            if (breakLoop)
                            {
                                break;
                            }
                        }
                        if (!isHangingNode)
                        {
                            bool regularNode = false;

                            for (const auto& intersection2
                                 : intersections(problem_.gridView(), element2))
                            {
                                for (const auto& intersection4
                                     : intersections(problem_.gridView(), element4))
                                {
                                    if (intersection4.neighbor())
                                    {
                                        auto element41 = intersection4.outside();

                                        if (element41 == element && element41.level() > element.level())
                                        {
                                            //adjust values of intersection12 in case of hanging nodes
                                            globalPosFace41 = intersection4.geometry().center();
                                            faceVol41 = intersection4.geometry().volume() / 2.0;

                                            interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol41, 0, 1);
                                            interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace41, 0,
                                                    1);
                                            interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol41, 3, 0);
                                            interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace41, 3,
                                                    0);
                                        }
                                    }

                                    if (intersection2.neighbor() && intersection4.neighbor())
                                    {
                                        auto element32 = intersection2.outside();
                                        auto element34 = intersection4.outside();

                                        //hanging node!
                                        if (element32 == element34 && element32 != element)
                                        {
                                            //store pointer 3
                                            interactionVolumes_[globalVertIdx1234].setSubVolumeElement(element32,
                                                    2);

                                            interactionVolumes_[globalVertIdx1234].setIndexOnElement(
                                                    intersection2.indexInInside(), 1, 0);
                                            interactionVolumes_[globalVertIdx1234].setIndexOnElement(
                                                    intersection2.indexInOutside(), 2, 1);
                                            interactionVolumes_[globalVertIdx1234].setIndexOnElement(
                                                    intersection4.indexInInside(), 3, 1);
                                            interactionVolumes_[globalVertIdx1234].setIndexOnElement(
                                                    intersection4.indexInOutside(), 2, 0);

                                            globalPosFace23 = intersection2.geometry().center();
                                            globalPosFace34 = intersection4.geometry().center();

                                            Scalar faceVol23 = intersection2.geometry().volume() / 2.0;
                                            Scalar faceVol34 = intersection4.geometry().volume() / 2.0;

                                            // get outer normal vector scaled with half volume of face : for numbering of n see Aavatsmark, Eigestad
                                            DimVector unitOuterNormal23 = intersection2.centerUnitOuterNormal();
                                            DimVector unitOuterNormal43 = intersection4.centerUnitOuterNormal();

                                            interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal23, 1, 0);
                                            unitOuterNormal23 *= -1;
                                            interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal23, 2, 1);
                                            interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal43, 3, 1);
                                            unitOuterNormal43 *= -1;
                                            interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal43, 2, 0);

                                            if (element32.level() > element2.level())
                                            {
                                                for (const auto& intersection3
                                                     : intersections(problem_.gridView(), element32))
                                                {
                                                    if (intersection3.neighbor())
                                                    {
                                                        auto element23 = intersection3.outside();

                                                        if (element23 == element2)
                                                        {
                                                            globalPosFace23 = intersection3.geometry().center();
                                                            faceVol23 = intersection3.geometry().volume() / 2.0;

                                                        }
                                                    }
                                                }
                                            }
                                            interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol23, 1, 0);
                                            interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol23, 2, 1);
                                            interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace23, 1,
                                                    0);
                                            interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace23, 2,
                                                    1);

                                            if (element34.level() > element4.level())
                                            {
                                                for (const auto& intersection3
                                                     : intersections(problem_.gridView(), element34))
                                                {
                                                    if (intersection3.neighbor())
                                                    {
                                                        auto element43 = intersection3.outside();

                                                        if (element43 == element4)
                                                        {
                                                            globalPosFace34 = intersection3.geometry().center();
                                                            faceVol34 = intersection3.geometry().volume() / 2.0;

                                                        }
                                                    }
                                                }
                                            }

                                            interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol34, 2, 0);
                                            interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol34, 3, 1);
                                            interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace34, 2,
                                                    0);
                                            interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace34, 3,
                                                    1);

                                            regularNode = true;
                                            break;
                                        }
                                    }
                                    if (regularNode)
                                    {
                                        break;
                                    }
                                }
                            }
                            if (!regularNode)
                            {
                                interactionVolumes_[globalVertIdx1234].reset();
                                continue;
                            }
                        }
                    }
                    else
                    {
                        bool regularNode = false;
                        for (const auto& intersection2
                             : intersections(problem_.gridView(), element2))
                        {
                            for (const auto& intersection4
                                 : intersections(problem_.gridView(), element4))
                            {
                                if (intersection4.neighbor())
                                {
                                    auto element41 = intersection4.outside();

                                    if (element41 == element && element41.level() > element.level())
                                    {
                                        //adjust values of intersection12 in case of hanging nodes
                                        globalPosFace41 = intersection4.geometry().center();
                                        faceVol41 = intersection4.geometry().volume() / 2.0;

                                        interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol41, 0, 1);
                                        interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace41, 0, 1);
                                        interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol41, 3, 0);
                                        interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace41, 3, 0);
                                    }
                                }

                                if (intersection2.neighbor() && intersection4.neighbor())
                                {
                                    auto element32 = intersection2.outside();
                                    auto element34 = intersection4.outside();

                                    // find the common neighbor cell between cell 2 and cell 3, except cell 1
                                    if (element32 == element34 && element32 != element)
                                    {
                                        //store pointer 3
                                        interactionVolumes_[globalVertIdx1234].setSubVolumeElement(element32, 2);

                                        interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection2.indexInInside(),
                                                1, 0);
                                        interactionVolumes_[globalVertIdx1234].setIndexOnElement(
                                                intersection2.indexInOutside(), 2, 1);
                                        interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection4.indexInInside(),
                                                3, 1);
                                        interactionVolumes_[globalVertIdx1234].setIndexOnElement(
                                                intersection4.indexInOutside(), 2, 0);

                                        globalPosFace23 = intersection2.geometry().center();
                                        globalPosFace34 = intersection4.geometry().center();

                                        Scalar faceVol23 = intersection2.geometry().volume() / 2.0;
                                        Scalar faceVol34 = intersection4.geometry().volume() / 2.0;

                                        // get outer normal vector scaled with half volume of face : for numbering of n see Aavatsmark, Eigestad
                                        DimVector unitOuterNormal23 = intersection2.centerUnitOuterNormal();

                                        DimVector unitOuterNormal43 = intersection4.centerUnitOuterNormal();

                                        interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal23, 1, 0);
                                        unitOuterNormal23 *= -1;
                                        interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal23, 2, 1);
                                        interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal43, 3, 1);
                                        unitOuterNormal43 *= -1;
                                        interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal43, 2, 0);
                                        interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol23, 1, 0);
                                        interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol23, 2, 1);
                                        interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol34, 2, 0);
                                        interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol34, 3, 1);
                                        interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace23, 1, 0);
                                        interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace23, 2, 1);
                                        interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace34, 2, 0);
                                        interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace34, 3, 1);

                                        regularNode = true;
                                        break;
                                    }
                                }
                            }
                            if (regularNode)
                            {
                                break;
                            }
                        }
                        if (!regularNode)
                        {
                            interactionVolumes_[globalVertIdx1234].reset();
                            continue;
                        }
                    }

                }

                // 'intersection14' is on the boundary
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

                    finished = false;

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
                                    interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace23, 1, 0);

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


                                    if (element.level() == element2.level())
                                    {
                                        innerBoundaryVolumeFaces_[eIdxGlobal1][intersection12.indexInInside()] = true;
                                        innerBoundaryVolumeFaces_[eIdxGlobal2][intersection12.indexInOutside()] = true;
                                    }
                                    else if (element.level() < element2.level())
                                    {
                                        innerBoundaryVolumeFaces_[eIdxGlobal2][intersection12.indexInOutside()] = true;
                                    }
                                    else if (element.level() > element2.level())
                                    {
                                        innerBoundaryVolumeFaces_[eIdxGlobal1][intersection12.indexInInside()] = true;
                                    }

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

            // handle boundary face 'intersection12'
            else
            {
                auto intersection14 = getNextIntersection_(element, isIt12);

                int indexInInside14 = intersection14.indexInInside();

                //get center vertex
                GlobalPosition corner1234(0);
                int globalVertIdx1234 = 0;

                // get the global coordinate and global vertex index of corner1234
                for (int i = 0; i < intersection12.geometry().corners(); ++i)
                {
                    bool finished = false;

                    int localVertIdx12corner = refElement.subEntity(indexInInside12, dim - 1, i, dim);

                    int globalVertIdx12corner = problem_.variables().index(element.template subEntity<dim>(localVertIdx12corner));

                    for (int j = 0; j < intersection14.geometry().corners(); ++j)
                    {
                        int localVertIdx14corner = refElement.subEntity(indexInInside14, dim - 1, j, dim);

                        int globalVertIdx14corner = problem_.variables().index(element.template subEntity<dim>(localVertIdx14corner));

                        if (globalVertIdx12corner == globalVertIdx14corner)
                        {
                            corner1234 = element.geometry().corner(localVertIdx12corner);

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
                }

                interactionVolumes_[globalVertIdx1234].setCenterPosition(corner1234);

                //store pointer 1
                interactionVolumes_[globalVertIdx1234].setSubVolumeElement(element, 0);
                interactionVolumes_[globalVertIdx1234].setIndexOnElement(indexInInside12, 0, 0);
                interactionVolumes_[globalVertIdx1234].setIndexOnElement(indexInInside14, 0, 1);

                // center of face in global coordinates, i.e., the midpoint of edge 'intersection12'
                const GlobalPosition& globalPosFace12 = intersection12.geometry().center();

                // get face volume
                Scalar faceVol12 = intersection12.geometry().volume() / 2.0;

                // get outer normal vector scaled with half volume of face 'intersection12'
                DimVector unitOuterNormal12 = intersection12.centerUnitOuterNormal();

                // center of face in global coordinates, i.e., the midpoint of edge 'intersection14'
                const GlobalPosition& globalPosFace41 = intersection14.geometry().center();

                // get face volume
                Scalar faceVol41 = intersection14.geometry().volume() / 2.0;

                // get outer normal vector scaled with half volume of face 'intersection14': for numbering of n see Aavatsmark, Eigestad
                DimVector unitOuterNormal14 = intersection14.centerUnitOuterNormal();

                interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal12, 0, 0);
                interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal14, 0, 1);
                //get the normals of from cell 2 and 4
                unitOuterNormal14 *= -1;
                unitOuterNormal12 *= -1;
                interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol12, 0, 0);
                interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol41, 0, 1);
                interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace12, 0, 0);
                interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace41, 0, 1);

                problem_.boundaryTypes(bcType, intersection12);
                PrimaryVariables boundValues(0.0);

                interactionVolumes_[globalVertIdx1234].setBoundary(bcType, 0);
                if (bcType.isNeumann(pressEqIdx))
                {
                    problem_.neumann(boundValues, intersection12);
                    boundValues *= faceVol12;
                    interactionVolumes_[globalVertIdx1234].setNeumannCondition(boundValues, 0);
                }
                if (bcType.hasDirichlet())
                {
                    problem_.dirichlet(boundValues, intersection12);
                    interactionVolumes_[globalVertIdx1234].setDirichletCondition(boundValues, 0);
                }

                // 'intersection14' is on boundary
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

                // 'intersection14' is inside
                else if (intersection14.neighbor())
                {
                    // neighbor cell 3
                    // access neighbor cell 3
                    auto element4 = intersection14.outside();
                    int eIdxGlobal4 = problem_.variables().index(element4);

                    if ((element.level() == element4.level() && eIdxGlobal1 > eIdxGlobal4)
                            || element.level() < element4.level())
                    {
                        interactionVolumes_[globalVertIdx1234].reset();
                        continue;
                    }

                    interactionVolumes_[globalVertIdx1234].setIndexOnElement(intersection14.indexInOutside(), 3, 0);

                    //store pointer 4
                    interactionVolumes_[globalVertIdx1234].setSubVolumeElement(element4, 3);

                    interactionVolumes_[globalVertIdx1234].setNormal(unitOuterNormal14, 3, 0);
                    interactionVolumes_[globalVertIdx1234].setFaceArea(faceVol41, 3, 0);
                    interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace41, 3, 0);

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
                                    interactionVolumes_[globalVertIdx1234].setFacePosition(globalPosFace34, 3, 1);

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

                                    if (element.level() == element4.level())
                                    {
                                        innerBoundaryVolumeFaces_[eIdxGlobal1][intersection14.indexInInside()] = true;
                                        innerBoundaryVolumeFaces_[eIdxGlobal4][intersection14.indexInOutside()] = true;
                                    }
                                    if (element.level() < element4.level())
                                    {
                                         innerBoundaryVolumeFaces_[eIdxGlobal4][intersection14.indexInOutside()] = true;
                                    }
                                    if (element.level() > element4.level())
                                    {
                                        innerBoundaryVolumeFaces_[eIdxGlobal1][intersection14.indexInInside()] = true;
                                    }


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
                                "fvmpfao2pfaboundpressure2p_adaptive.hh: boundary shape not available as interaction volume shape");
                    }
                }
                else
                {
                    DUNE_THROW(Dune::NotImplemented,
                            "fvmpfao2pfaboundpressure2p_adaptive.hh: interface type not supported!");
                }
            }

        } // end all intersections
    } // end grid traversal

    return;
}

// TODO doc me!
template<class TypeTag>
void FvMpfaL2dPressure2pAdaptive<TypeTag>::printInteractionVolumes()
{
    for (const auto& vertex : vertices(problem_.gridView()))
    {
        int vIdxGlobal = problem_.variables().index(vertex);

        InteractionVolume& interactionVolume = interactionVolumes_[vIdxGlobal];

        if (interactionVolume.getElementNumber() > 2)
        {
            interactionVolume.printInteractionVolumeInfo();
            std::cout << "global vertex index: " << vIdxGlobal << "\n";
            if (interactionVolume.getElementNumber() == 3)
            {
                auto element1 = interactionVolume.getSubVolumeElement(0);
                auto element2 = interactionVolume.getSubVolumeElement(1);
                auto element4 = interactionVolume.getSubVolumeElement(3);

                int eIdxGlobal1 = problem_.variables().index(element1);
                int eIdxGlobal2 = problem_.variables().index(element2);
                int eIdxGlobal4 = problem_.variables().index(element4);

                std::cout << "global element index 1: " << eIdxGlobal1 << "\n";
                std::cout << "global element index 2: " << eIdxGlobal2 << "\n";
                std::cout << "global element index 4: " << eIdxGlobal4 << "\n";
            }
            if (interactionVolume.getElementNumber() == 4)
            {
                auto element1 = interactionVolume.getSubVolumeElement(0);
                auto element2 = interactionVolume.getSubVolumeElement(1);
                auto element3 = interactionVolume.getSubVolumeElement(2);
                auto element4 = interactionVolume.getSubVolumeElement(3);

                int eIdxGlobal1 = problem_.variables().index(element1);
                int eIdxGlobal2 = problem_.variables().index(element2);
                int eIdxGlobal3 = problem_.variables().index(element3);
                int eIdxGlobal4 = problem_.variables().index(element4);

                std::cout << "global element index 1: " << eIdxGlobal1 << "\n";
                std::cout << "global element index 2: " << eIdxGlobal2 << "\n";
                std::cout << "global element index 3: " << eIdxGlobal3 << "\n";
                std::cout << "global element index 4: " << eIdxGlobal4 << "\n";
            }
        }
    }
}

// only for 2-D general quadrilateral
// TODO doc me!
template<class TypeTag>
void FvMpfaL2dPressure2pAdaptive<TypeTag>::assemble()
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
            if (interactionVolume.getElementNumber() == 4)
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
                this->f_[eIdxGlobal1] += volume1 / (4.0) * (source[wPhaseIdx] / density_[wPhaseIdx] + source[nPhaseIdx] / density_[nPhaseIdx]);
                problem_.source(source, element2);
                this->f_[eIdxGlobal2] += volume2 / (4.0) * (source[wPhaseIdx] / density_[wPhaseIdx] + source[nPhaseIdx] / density_[nPhaseIdx]);
                problem_.source(source, element3);
                this->f_[eIdxGlobal3] += volume3 / (4.0) * (source[wPhaseIdx] / density_[wPhaseIdx] + source[nPhaseIdx] / density_[nPhaseIdx]);
                problem_.source(source, element4);
                this->f_[eIdxGlobal4] += volume4 / (4.0) * (source[wPhaseIdx] / density_[wPhaseIdx] + source[nPhaseIdx] / density_[nPhaseIdx]);

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

                std::vector < DimVector > lambda(2 * dim);

                lambda[0][0] = lambdaTotal1;
                lambda[0][1] = lambdaTotal1;
                lambda[1][0] = lambdaTotal2;
                lambda[1][1] = lambdaTotal2;
                lambda[2][0] = lambdaTotal3;
                lambda[2][1] = lambdaTotal3;
                lambda[3][0] = lambdaTotal4;
                lambda[3][1] = lambdaTotal4;

                //add capillary pressure and gravity terms to right-hand-side
                //calculate capillary pressure velocity
                Dune::FieldVector<Scalar, 2 * dim> pc(0);
                pc[0] = cellData1.capillaryPressure();
                pc[1] = cellData2.capillaryPressure();
                pc[2] = cellData3.capillaryPressure();
                pc[3] = cellData4.capillaryPressure();

                Dune::FieldVector<Scalar, 2 * dim> gravityDiff(0);

                gravityDiff[0] = (problem_.bBoxMax() - globalPos1) * gravity_ * (density_[nPhaseIdx] - density_[wPhaseIdx]);
                gravityDiff[1] = (problem_.bBoxMax() - globalPos2) * gravity_ * (density_[nPhaseIdx] - density_[wPhaseIdx]);
                gravityDiff[2] = (problem_.bBoxMax() - globalPos3) * gravity_ * (density_[nPhaseIdx] - density_[wPhaseIdx]);
                gravityDiff[3] = (problem_.bBoxMax() - globalPos4) * gravity_ * (density_[nPhaseIdx] - density_[wPhaseIdx]);

                pc += gravityDiff;

                Dune::FieldVector<Scalar, 2 * dim> pcFlux(0);

                Scalar pcPotential12 = 0;
                Scalar pcPotential14 = 0;
                Scalar pcPotential32 = 0;
                Scalar pcPotential34 = 0;

                DimVector Tu(0);
                Dune::FieldVector<Scalar, 2 * dim - dim + 1> u(0);
                Dune::FieldMatrix<Scalar, dim, 2 * dim - dim + 1> T(0);

                int transmissibilityType = transmissibilityCalculator_.calculateTransmissibility(T, interactionVolume, lambda, 0, 1, 2, 3);

                if (transmissibilityType == TransmissibilityCalculator::rightTriangle)
                {
                    if (innerBoundaryVolumeFaces_[eIdxGlobal1][interactionVolume.getIndexOnElement(0, 0)])
                    {
                        T *= 2;
                    }
                    this->A_[eIdxGlobal1][eIdxGlobal2] += T[1][0];
                    this->A_[eIdxGlobal1][eIdxGlobal3] += T[1][1];
                    this->A_[eIdxGlobal1][eIdxGlobal1] += T[1][2];

                    this->A_[eIdxGlobal2][eIdxGlobal2] -= T[1][0];
                    this->A_[eIdxGlobal2][eIdxGlobal3] -= T[1][1];
                    this->A_[eIdxGlobal2][eIdxGlobal1] -= T[1][2];

                    u[0] = pc[1];
                    u[1] = pc[2];
                    u[2] = pc[0];

                    T.mv(u, Tu);

                    pcFlux[0] = Tu[1];
                    pcPotential12 = Tu[1];
                }
                else if (transmissibilityType == TransmissibilityCalculator::leftTriangle)
                {
                    if (innerBoundaryVolumeFaces_[eIdxGlobal1][interactionVolume.getIndexOnElement(0, 0)])
                    {
                        T *= 2;
                    }
                    this->A_[eIdxGlobal1][eIdxGlobal1] += T[1][0];
                    this->A_[eIdxGlobal1][eIdxGlobal4] += T[1][1];
                    this->A_[eIdxGlobal1][eIdxGlobal2] += T[1][2];

                    this->A_[eIdxGlobal2][eIdxGlobal1] -= T[1][0];
                    this->A_[eIdxGlobal2][eIdxGlobal4] -= T[1][1];
                    this->A_[eIdxGlobal2][eIdxGlobal2] -= T[1][2];

                    u[0] = pc[0];
                    u[1] = pc[3];
                    u[2] = pc[1];

                    T.mv(u, Tu);

                    pcFlux[0] = Tu[1];
                    pcPotential12 = Tu[1];
                }

                transmissibilityType = transmissibilityCalculator_.calculateTransmissibility(T, interactionVolume, lambda, 1, 2, 3, 0);

                if (transmissibilityType == TransmissibilityCalculator::rightTriangle)
                {
                    if (innerBoundaryVolumeFaces_[eIdxGlobal2][interactionVolume.getIndexOnElement(1, 0)])
                    {
                        T *= 2;
                    }
                    this->A_[eIdxGlobal2][eIdxGlobal3] += T[1][0];
                    this->A_[eIdxGlobal2][eIdxGlobal4] += T[1][1];
                    this->A_[eIdxGlobal2][eIdxGlobal2] += T[1][2];

                    this->A_[eIdxGlobal3][eIdxGlobal3] -= T[1][0];
                    this->A_[eIdxGlobal3][eIdxGlobal4] -= T[1][1];
                    this->A_[eIdxGlobal3][eIdxGlobal2] -= T[1][2];

                    u[0] = pc[2];
                    u[1] = pc[3];
                    u[2] = pc[1];

                    T.mv(u, Tu);

                    pcFlux[1] = Tu[1];
                    pcPotential32 = -Tu[1];
                }
                else if (transmissibilityType == TransmissibilityCalculator::leftTriangle)
                {
                    if (innerBoundaryVolumeFaces_[eIdxGlobal2][interactionVolume.getIndexOnElement(1, 0)])
                    {
                        T *= 2;
                    }
                    this->A_[eIdxGlobal2][eIdxGlobal2] += T[1][0];
                    this->A_[eIdxGlobal2][eIdxGlobal1] += T[1][1];
                    this->A_[eIdxGlobal2][eIdxGlobal3] += T[1][2];

                    this->A_[eIdxGlobal3][eIdxGlobal2] -= T[1][0];
                    this->A_[eIdxGlobal3][eIdxGlobal1] -= T[1][1];
                    this->A_[eIdxGlobal3][eIdxGlobal3] -= T[1][2];

                    u[0] = pc[1];
                    u[1] = pc[0];
                    u[2] = pc[2];

                    T.mv(u, Tu);

                    pcFlux[1] = Tu[1];
                    pcPotential32 = -Tu[1];
                }

                transmissibilityType = transmissibilityCalculator_.calculateTransmissibility(T, interactionVolume, lambda, 2, 3, 0, 1);

                if (transmissibilityType == TransmissibilityCalculator::rightTriangle)
                {
                    if (innerBoundaryVolumeFaces_[eIdxGlobal3][interactionVolume.getIndexOnElement(2, 0)])
                    {
                        T *= 2;
                    }
                    this->A_[eIdxGlobal3][eIdxGlobal4] += T[1][0];
                    this->A_[eIdxGlobal3][eIdxGlobal1] += T[1][1];
                    this->A_[eIdxGlobal3][eIdxGlobal3] += T[1][2];

                    this->A_[eIdxGlobal4][eIdxGlobal4] -= T[1][0];
                    this->A_[eIdxGlobal4][eIdxGlobal1] -= T[1][1];
                    this->A_[eIdxGlobal4][eIdxGlobal3] -= T[1][2];

                    u[0] = pc[3];
                    u[1] = pc[0];
                    u[2] = pc[2];

                    T.mv(u, Tu);

                    pcFlux[2] = Tu[1];
                    pcPotential34 = Tu[1];
                }
                else if (transmissibilityType == TransmissibilityCalculator::leftTriangle)
                {
                    if (innerBoundaryVolumeFaces_[eIdxGlobal3][interactionVolume.getIndexOnElement(2, 0)])
                    {
                        T *= 2;
                    }
                    this->A_[eIdxGlobal3][eIdxGlobal3] += T[1][0];
                    this->A_[eIdxGlobal3][eIdxGlobal2] += T[1][1];
                    this->A_[eIdxGlobal3][eIdxGlobal4] += T[1][2];

                    this->A_[eIdxGlobal4][eIdxGlobal3] -= T[1][0];
                    this->A_[eIdxGlobal4][eIdxGlobal2] -= T[1][1];
                    this->A_[eIdxGlobal4][eIdxGlobal4] -= T[1][2];

                    u[0] = pc[2];
                    u[1] = pc[1];
                    u[2] = pc[3];

                    T.mv(u, Tu);

                    pcFlux[2] = Tu[1];
                    pcPotential34 = Tu[1];
                }

                transmissibilityType = transmissibilityCalculator_.calculateTransmissibility(T, interactionVolume, lambda, 3, 0, 1, 2);

                if (transmissibilityType == TransmissibilityCalculator::rightTriangle)
                {
                    if (innerBoundaryVolumeFaces_[eIdxGlobal4][interactionVolume.getIndexOnElement(3, 0)])
                    {
                        T *= 2;
                    }
                    this->A_[eIdxGlobal4][eIdxGlobal1] += T[1][0];
                    this->A_[eIdxGlobal4][eIdxGlobal2] += T[1][1];
                    this->A_[eIdxGlobal4][eIdxGlobal4] += T[1][2];

                    this->A_[eIdxGlobal1][eIdxGlobal1] -= T[1][0];
                    this->A_[eIdxGlobal1][eIdxGlobal2] -= T[1][1];
                    this->A_[eIdxGlobal1][eIdxGlobal4] -= T[1][2];

                    u[0] = pc[0];
                    u[1] = pc[1];
                    u[2] = pc[3];

                    T.mv(u, Tu);

                    pcFlux[3] = Tu[1];
                    pcPotential14 = -Tu[1];
                }
                else if (transmissibilityType == TransmissibilityCalculator::leftTriangle)
                {
                    if (innerBoundaryVolumeFaces_[eIdxGlobal4][interactionVolume.getIndexOnElement(3, 0)])
                    {
                        T *= 2;
                    }
                    this->A_[eIdxGlobal4][eIdxGlobal4] += T[1][0];
                    this->A_[eIdxGlobal4][eIdxGlobal3] += T[1][1];
                    this->A_[eIdxGlobal4][eIdxGlobal1] += T[1][2];

                    this->A_[eIdxGlobal1][eIdxGlobal4] -= T[1][0];
                    this->A_[eIdxGlobal1][eIdxGlobal3] -= T[1][1];
                    this->A_[eIdxGlobal1][eIdxGlobal1] -= T[1][2];

                    u[0] = pc[3];
                    u[1] = pc[2];
                    u[2] = pc[0];

                    T.mv(u, Tu);

                    pcFlux[3] = Tu[1];
                    pcPotential14 = -Tu[1];
                }

                if (pc[0] == 0 && pc[1] == 0 && pc[2] == 0 && pc[3] == 0)
                {
                    continue;
                }

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

                    switch (pressureType_)
                    {
                    case pw:
                    {
                        if (i == nPhaseIdx)
                        {
                            //add capillary pressure term to right hand side
                            this->f_[eIdxGlobal1] -= (pcFluxReal[0] - pcFluxReal[3]);
                            this->f_[eIdxGlobal2] -= (pcFluxReal[1] - pcFluxReal[0]);
                            this->f_[eIdxGlobal3] -= (pcFluxReal[2] - pcFluxReal[1]);
                            this->f_[eIdxGlobal4] -= (pcFluxReal[3] - pcFluxReal[2]);

                        }
                        break;
                    }
                    case pn:
                    {
                        if (i == wPhaseIdx)
                        {
                            //add capillary pressure term to right hand side
                            this->f_[eIdxGlobal1] += (pcFluxReal[0] - pcFluxReal[3]);
                            this->f_[eIdxGlobal2] += (pcFluxReal[1] - pcFluxReal[0]);
                            this->f_[eIdxGlobal3] += (pcFluxReal[2] - pcFluxReal[1]);
                            this->f_[eIdxGlobal4] += (pcFluxReal[3] - pcFluxReal[2]);
                        }
                        break;
                    }
                    }
                }
            }
            else if (interactionVolume.getElementNumber() == 3)
            {
                auto element1 = interactionVolume.getSubVolumeElement(0);
                auto element2 = interactionVolume.getSubVolumeElement(1);
                auto element4 = interactionVolume.getSubVolumeElement(3);

                // get global coordinate of cell centers
                const GlobalPosition& globalPos1 = element1.geometry().center();
                const GlobalPosition& globalPos2 = element2.geometry().center();
                const GlobalPosition& globalPos4 = element4.geometry().center();

                // cell volumes
                Scalar volume1 = element1.geometry().volume();
                Scalar volume2 = element2.geometry().volume();

                // cell index
                int eIdxGlobal1 = problem_.variables().index(element1);
                int eIdxGlobal2 = problem_.variables().index(element2);
                int eIdxGlobal4 = problem_.variables().index(element4);

                //get the cell Data
                CellData& cellData1 = problem_.variables().cellData(eIdxGlobal1);
                CellData& cellData2 = problem_.variables().cellData(eIdxGlobal2);
                CellData& cellData4 = problem_.variables().cellData(eIdxGlobal4);

                // evaluate right hand side -> only add source for the cells without hanging node!
                // In doing so every cell gets the source from 4 vertices and the division by 4 is correct!
                PrimaryVariables source(0.0);
                problem_.source(source, element1);
                this->f_[eIdxGlobal1] += volume1 / (4.0) * (source[wPhaseIdx] / density_[wPhaseIdx] + source[nPhaseIdx] / density_[nPhaseIdx]);
                problem_.source(source, element2);
                this->f_[eIdxGlobal2] += volume2 / (4.0) * (source[wPhaseIdx] / density_[wPhaseIdx] + source[nPhaseIdx] / density_[nPhaseIdx]);

                this->f_[eIdxGlobal1] += evaluateErrorTerm_(cellData1) * volume1 / (4.0);
                this->f_[eIdxGlobal2] += evaluateErrorTerm_(cellData2) * volume2 / (4.0);

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
                Dune::FieldVector<Scalar, numPhases> lambda4(cellData4.mobility(wPhaseIdx));
                lambda4[nPhaseIdx] = cellData4.mobility(nPhaseIdx);

                //compute total mobility of cell 1
                Scalar lambdaTotal4 = lambda4[wPhaseIdx] + lambda4[nPhaseIdx];

                std::vector < DimVector > lambda(4);

                lambda[0][0] = lambdaTotal1;
                lambda[0][1] = lambdaTotal1;
                lambda[1][0] = lambdaTotal2;
                lambda[1][1] = lambdaTotal2;
                lambda[3][0] = lambdaTotal4;
                lambda[3][1] = lambdaTotal4;

                //add capillary pressure and gravity terms to right-hand-side
                //calculate capillary pressure velocity
                Dune::FieldVector<Scalar, 3> pc(0);
                pc[0] = cellData1.capillaryPressure();
                pc[1] = cellData2.capillaryPressure();
                pc[2] = cellData4.capillaryPressure();

                Dune::FieldVector<Scalar, 3> gravityDiff(0);

                gravityDiff[0] = (problem_.bBoxMax() - globalPos1) * gravity_ * (density_[nPhaseIdx] - density_[wPhaseIdx]);
                gravityDiff[1] = (problem_.bBoxMax() - globalPos2) * gravity_ * (density_[nPhaseIdx] - density_[wPhaseIdx]);
                gravityDiff[2] = (problem_.bBoxMax() - globalPos4) * gravity_ * (density_[nPhaseIdx] - density_[wPhaseIdx]);

                pc += gravityDiff;

                Dune::FieldVector<Scalar, 3> pcFlux(0);

                Scalar pcPotential12 = 0;
                Scalar pcPotential14 = 0;
                Scalar pcPotential24 = 0;

                DimVector Tu(0);
                Dune::FieldVector<Scalar, 2 * dim - dim + 1> u(0);
                Dune::FieldMatrix<Scalar, dim, 2 * dim - dim + 1> T(0);

                int transmissibilityType = transmissibilityCalculator_.calculateTransmissibility(T, interactionVolume, lambda, 0, 1, 3, 3);

                if (transmissibilityType == TransmissibilityCalculator::rightTriangle)
                {
                    if (innerBoundaryVolumeFaces_[eIdxGlobal1][interactionVolume.getIndexOnElement(0, 0)])
                    {
                        T *= 2;
                    }
                    this->A_[eIdxGlobal1][eIdxGlobal2] += T[1][0];
                    this->A_[eIdxGlobal1][eIdxGlobal4] += T[1][1];
                    this->A_[eIdxGlobal1][eIdxGlobal1] += T[1][2];

                    this->A_[eIdxGlobal2][eIdxGlobal2] -= T[1][0];
                    this->A_[eIdxGlobal2][eIdxGlobal4] -= T[1][1];
                    this->A_[eIdxGlobal2][eIdxGlobal1] -= T[1][2];

                    u[0] = pc[1];
                    u[1] = pc[2];
                    u[2] = pc[0];

                    T.mv(u, Tu);

                    pcFlux[0] = Tu[1];
                    pcPotential12 = Tu[1];
                }
                else if (transmissibilityType == TransmissibilityCalculator::leftTriangle)
                {
                    if (innerBoundaryVolumeFaces_[eIdxGlobal1][interactionVolume.getIndexOnElement(0, 0)])
                    {
                        T *= 2;
                    }
                    this->A_[eIdxGlobal1][eIdxGlobal1] += T[1][0];
                    this->A_[eIdxGlobal1][eIdxGlobal4] += T[1][1];
                    this->A_[eIdxGlobal1][eIdxGlobal2] += T[1][2];

                    this->A_[eIdxGlobal2][eIdxGlobal1] -= T[1][0];
                    this->A_[eIdxGlobal2][eIdxGlobal4] -= T[1][1];
                    this->A_[eIdxGlobal2][eIdxGlobal2] -= T[1][2];

                    u[0] = pc[0];
                    u[1] = pc[2];
                    u[2] = pc[1];

                    T.mv(u, Tu);

                    pcFlux[0] = Tu[1];
                    pcPotential12 = Tu[1];
                }
                else
                {
                    DUNE_THROW(Dune::RangeError, "Could not calculate Tranmissibility!");
                }

                transmissibilityType = transmissibilityCalculator_.calculateLeftHNTransmissibility(T, interactionVolume, lambda, 1, 3, 0);

                if (transmissibilityType == TransmissibilityCalculator::leftTriangle)
                {
                    if (innerBoundaryVolumeFaces_[eIdxGlobal2][interactionVolume.getIndexOnElement(1, 0)])
                    {
                        T *= 2;
                    }
                    this->A_[eIdxGlobal2][eIdxGlobal2] += T[1][0];
                    this->A_[eIdxGlobal2][eIdxGlobal1] += T[1][1];
                    this->A_[eIdxGlobal2][eIdxGlobal4] += T[1][2];

                    this->A_[eIdxGlobal4][eIdxGlobal2] -= T[1][0];
                    this->A_[eIdxGlobal4][eIdxGlobal1] -= T[1][1];
                    this->A_[eIdxGlobal4][eIdxGlobal4] -= T[1][2];

                    u[0] = pc[1];
                    u[1] = pc[0];
                    u[2] = pc[2];

                    T.mv(u, Tu);

                    pcFlux[1] = Tu[1];
                    pcPotential24 = Tu[1];
                }
                else
                {
                    DUNE_THROW(Dune::RangeError, "Could not calculate Tranmissibility!");
                }

                transmissibilityType = transmissibilityCalculator_.calculateRightHNTransmissibility(T, interactionVolume, lambda, 3, 0, 1);

                if (transmissibilityType == TransmissibilityCalculator::rightTriangle)
                {
                    if (innerBoundaryVolumeFaces_[eIdxGlobal1][interactionVolume.getIndexOnElement(0, 1)])
                    {
                        T *= 2;
                    }

                    this->A_[eIdxGlobal4][eIdxGlobal1] += T[1][0];
                    this->A_[eIdxGlobal4][eIdxGlobal2] += T[1][1];
                    this->A_[eIdxGlobal4][eIdxGlobal4] += T[1][2];

                    this->A_[eIdxGlobal1][eIdxGlobal1] -= T[1][0];
                    this->A_[eIdxGlobal1][eIdxGlobal2] -= T[1][1];
                    this->A_[eIdxGlobal1][eIdxGlobal4] -= T[1][2];

                    u[0] = pc[0];
                    u[1] = pc[1];
                    u[2] = pc[2];

                    T.mv(u, Tu);

                    pcFlux[2] = Tu[1];
                    pcPotential14 = -Tu[1];
                }
                else
                {
                    DUNE_THROW(Dune::RangeError, "Could not calculate Tranmissibility!");
                }

                if (pc[0] == 0 && pc[1] == 0 && pc[2] == 0)
                {
                    continue;
                }

                //compute mobilities of face 1
                Dune::FieldVector<Scalar, numPhases> lambda12Upw(0.0);
                lambda12Upw[wPhaseIdx] = (pcPotential12 >= 0) ? lambda1[wPhaseIdx] : lambda2[wPhaseIdx];
                lambda12Upw[nPhaseIdx] = (pcPotential12 >= 0) ? lambda1[nPhaseIdx] : lambda2[nPhaseIdx];

                //compute mobilities of face 4
                Dune::FieldVector<Scalar, numPhases> lambda14Upw(0.0);
                lambda14Upw[wPhaseIdx] = (pcPotential14 >= 0) ? lambda1[wPhaseIdx] : lambda4[wPhaseIdx];
                lambda14Upw[nPhaseIdx] = (pcPotential14 >= 0) ? lambda1[nPhaseIdx] : lambda4[nPhaseIdx];

                //compute mobilities of face 2
                Dune::FieldVector<Scalar, numPhases> lambda24Upw(0.0);
                lambda24Upw[wPhaseIdx] = (pcPotential24 >= 0) ? lambda2[wPhaseIdx] : lambda4[wPhaseIdx];
                lambda24Upw[nPhaseIdx] = (pcPotential24 >= 0) ? lambda2[nPhaseIdx] : lambda4[nPhaseIdx];

                for (int i = 0; i < numPhases; i++)
                {
                    Scalar lambdaT12 = lambda12Upw[wPhaseIdx] + lambda12Upw[nPhaseIdx];
                    Scalar lambdaT14 = lambda14Upw[wPhaseIdx] + lambda14Upw[nPhaseIdx];
                    Scalar lambdaT24 = lambda24Upw[wPhaseIdx] + lambda24Upw[nPhaseIdx];
                    Scalar fracFlow12 = (lambdaT12 > threshold_) ? lambda12Upw[i] / (lambdaT12) : 0.0;
                    Scalar fracFlow14 = (lambdaT14 > threshold_) ? lambda14Upw[i] / (lambdaT14) : 0.0;
                    Scalar fracFlow24 = (lambdaT24 > threshold_) ? lambda24Upw[i] / (lambdaT24) : 0.0;

                    Dune::FieldVector<Scalar, 3> pcFluxReal(pcFlux);

                    pcFluxReal[0] *= fracFlow12;
                    pcFluxReal[1] *= fracFlow24;
                    pcFluxReal[2] *= fracFlow14;

                    switch (pressureType_)
                    {
                    case pw:
                    {
                        if (i == nPhaseIdx)
                        {
                            //add capillary pressure term to right hand side
                            this->f_[eIdxGlobal1] -= (pcFluxReal[0] - pcFluxReal[2]);
                            this->f_[eIdxGlobal2] -= (pcFluxReal[1] - pcFluxReal[0]);
                            this->f_[eIdxGlobal4] -= (pcFluxReal[2] - pcFluxReal[1]);

                        }
                        break;
                    }
                    case pn:
                    {
                        if (i == wPhaseIdx)
                        {
                            //add capillary pressure term to right hand side
                            this->f_[eIdxGlobal1] += (pcFluxReal[0] - pcFluxReal[2]);
                            this->f_[eIdxGlobal2] += (pcFluxReal[1] - pcFluxReal[0]);
                            this->f_[eIdxGlobal4] += (pcFluxReal[2] - pcFluxReal[1]);
                        }
                        break;
                    }
                    }
                }
            }
            else
            {
                DUNE_THROW(Dune::NotImplemented, "Unknown interactionvolume type!");
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
                this->f_[eIdxGlobal] += volume / (4.0) * (source[wPhaseIdx] / density_[wPhaseIdx] + source[nPhaseIdx] / density_[nPhaseIdx]);

                this->f_[eIdxGlobal] += evaluateErrorTerm_(cellData) * volume / (4.0);

                //get mobilities of the phases
                Dune::FieldVector<Scalar, numPhases> lambda(cellData.mobility(wPhaseIdx));
                lambda[nPhaseIdx] = cellData.mobility(nPhaseIdx);

                Scalar pc = cellData.capillaryPressure();

                Scalar gravityDiff = (problem_.bBoxMax() - globalPos) * gravity_ * (density_[nPhaseIdx] - density_[wPhaseIdx]);

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
                            Scalar J = interactionVolume.getNeumannValues(intVolFaceIdx)[wPhaseIdx] / density_[wPhaseIdx];
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
void FvMpfaL2dPressure2pAdaptive<TypeTag>::updateMaterialLaws()
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
