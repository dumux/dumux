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
 * \brief  Finite volume discretization of a two-phase flow pressure equation.
 */
#ifndef DUMUX_FVPRESSURE2P_ADAPTIVE_HH
#define DUMUX_FVPRESSURE2P_ADAPTIVE_HH

#include <dune/common/float_cmp.hh>

// dumux environment
#include <dumux/porousmediumflow/2p/sequential/properties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressure.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/velocity.hh>

namespace Dumux {

/*!
 * \ingroup SequentialTwoPModel
 * \brief Finite volume discretization of a two-phase flow pressure equation of the sequential IMPES model.
 *
 * Details see FVPressure2P
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag> class FVPressure2PAdaptive: public FVPressure2P<TypeTag>
{
    using ParentType = FVPressure2P<TypeTag>;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using PrimaryVariables = typename SolutionTypes::PrimaryVariables;
    using CellData = GetPropType<TypeTag, Properties::CellData>;

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
        sn = Indices::saturationNw
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx, numPhases = getPropValue<TypeTag, Properties::NumPhases>()
    };

    using Intersection = typename GridView::Intersection;

    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GravityVector = Dune::FieldVector<Scalar, dimWorld>;
    using FieldMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

protected:
    //! \cond \private
    using EntryType = typename ParentType::EntryType;
    enum
    {
        rhs = ParentType::rhs, matrix = ParentType::matrix
    };
    //! \endcond

public:
    // Function which calculates the flux entry
    void getFlux(EntryType&, const Intersection&, const CellData&, const bool);

    /*!
     * \brief Initializes the adaptive pressure model
     *
     * \copydetails FVPressure2P::initialize()
     */
    void initialize()
    {
        ParentType::initialize();

        if (!compressibility_)
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
        }

        velocity_.initialize();
    }

    /*!
     * \brief Pressure update
     *
     * \copydetails FVPressure::update()
     *
     * The grid-adaptive implementation also reconstructs the velocity directly after the pressure update.
     * This is necessary to make sure the hanging nodes are treated correctly!
     */
    void update()
    {
        int gridSize = problem_.gridView().size(0);
        // update RHS vector, matrix
        if (problem_.gridAdapt().wasAdapted())
        {
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
                case pGlobal:
                    this->pressure()[i] = cellData.globalPressure();
                    break;
                }
            }

            ParentType::initializeMatrix();
        }


        ParentType::update();

        calculateVelocity();

        return;
    }

    /*!
     * \brief Velocity calculation
     */
    void calculateVelocity()
    {
        velocity_.calculateVelocity();
    }

    /*!
     * \brief Velocity update
     */
    void updateVelocity()
    {
        ParentType::updateVelocity();

        calculateVelocity();
    }

    /*!
     * \brief Adds pressure output to the output file
     *
     *  \copydetails FVPressure2P::addOutputVtkFields(MultiWriter&)
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        ParentType::addOutputVtkFields(writer);

            velocity_.addOutputVtkFields(writer);
    }

    /*!
     * \brief Constructs a FVPressure2PAdaptive object
     *
     * \param problem A problem class object
     */
    FVPressure2PAdaptive(Problem& problem) :
            ParentType(problem), problem_(problem), velocity_(problem), gravity_(problem.gravity())
    {
        if (pressureType_ != pw && pressureType_ != pn && pressureType_ != pGlobal)
        {
            DUNE_THROW(Dune::NotImplemented, "Pressure type not supported!");
        }
        if (pressureType_ == pGlobal && compressibility_)
        {
            DUNE_THROW(Dune::NotImplemented, "Compressibility not supported for global pressure!");
        }
        if (saturationType_ != sw && saturationType_ != sn)
        {
            DUNE_THROW(Dune::NotImplemented, "Saturation type not supported!");
        }

        density_[wPhaseIdx] = 0.;
        density_[nPhaseIdx] = 0.;
        viscosity_[wPhaseIdx] = 0.;
        viscosity_[nPhaseIdx] = 0.;
    }

private:
    Problem& problem_;
    FVVelocity<TypeTag, GetPropType<TypeTag, Properties::Velocity>> velocity_;
    const GravityVector& gravity_; //!< vector including the gravity constant

    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];

    static const bool compressibility_ = getPropValue<TypeTag, Properties::EnableCompressibility>();
    //! gives kind of pressure used (\f$p_w\f$, \f$p_n\f$, \f$p_{global}\f$)
    static const int pressureType_ = getPropValue<TypeTag, Properties::PressureFormulation>();
    //! gives kind of saturation used (\f$S_w\f$, \f$S_n\f$)
    static const int saturationType_ = getPropValue<TypeTag, Properties::SaturationFormulation>();
};

/*!
 * \brief Function which calculates the flux entry.
 *
 * \copydetails FVPressure2P::getFlux(EntryType&,const Intersection&,const CellData&,const bool)
 *
 * Implementation of the getFlux() function for cell-cell interfaces with hanging nodes.
 * In case of hanging nodes the function does not return a vector of entry but directly manipulates the global matrix!
 */
template<class TypeTag>
void FVPressure2PAdaptive<TypeTag>::getFlux(EntryType& entry, const Intersection& intersection
        , const CellData& cellData, const bool first)
{
    auto elementI = intersection.inside();
    auto elementJ = intersection.outside();

    if (elementI.level() == elementJ.level() || dim == 3)
    {
        ParentType::getFlux(entry, intersection, cellData, first);

        // add the entry only once in case the VisitFacesOnlyOnce option is enabled!!!
        if (getPropValue<TypeTag, Properties::VisitFacesOnlyOnce>() && elementI.level() < elementJ.level())
        {
            entry = 0.;
        }
    }
    // hanging node situation: neighbor has higher level
    else if (elementJ.level() == elementI.level() + 1)
    {
        const CellData& cellDataJ = problem_.variables().cellData(problem_.variables().index(elementJ));

        // get global coordinates of cell centers
        const GlobalPosition& globalPosI = elementI.geometry().center();
        const GlobalPosition& globalPosJ = elementJ.geometry().center();

        int globalIdxI = problem_.variables().index(elementI);
        int globalIdxJ = problem_.variables().index(elementJ);

        // get mobilities and fractional flow factors
        Scalar lambdaWI = cellData.mobility(wPhaseIdx);
        Scalar lambdaNwI = cellData.mobility(nPhaseIdx);
        Scalar fractionalWI = cellData.fracFlowFunc(wPhaseIdx);
        Scalar fractionalNwI = cellData.fracFlowFunc(nPhaseIdx);
        Scalar lambdaWJ = cellDataJ.mobility(wPhaseIdx);
        Scalar lambdaNwJ = cellDataJ.mobility(nPhaseIdx);
        Scalar fractionalWJ = cellDataJ.fracFlowFunc(wPhaseIdx);
        Scalar fractionalNwJ = cellDataJ.fracFlowFunc(nPhaseIdx);

        // get capillary pressure
        Scalar pcI = cellData.capillaryPressure();
        Scalar pcJ = cellDataJ.capillaryPressure();

        // get face index
        int isIndexI = intersection.indexInInside();

        // get face normal
        const Dune::FieldVector<Scalar, dim>& unitOuterNormal = intersection.centerUnitOuterNormal();

        // get face area
        Scalar faceArea = intersection.geometry().volume();

        // Count number of hanging nodes
        // not really necessary
        int globalIdxK = 0;
        auto elementK = intersection.outside();
        // We are looking for two things:
        // IsIndexJ, the index of the interface from the neighbor-cell point of view
        // GlobalIdxK, the index of the third cell
        // Intersectioniterator around cell I
        for (const auto& intersectionI : intersections(problem_.gridView(), elementI))
        {
            if (intersectionI.neighbor())
            {
                auto neighbor2 = intersectionI.outside();

                // make sure we do not choose elemntI as third element
                // -> faces with hanging node have more than one intersection but only one face index!
                if (neighbor2 != elementJ && intersectionI.indexInInside() == isIndexI)
                {
                    globalIdxK = problem_.variables().index(neighbor2);
                    elementK = neighbor2;
                    break;
                }
            }
        }

        const CellData& cellDataK = problem_.variables().cellData(globalIdxK);

        // neighbor cell center in global coordinates
        const GlobalPosition& globalPosInterface = intersection.geometry().center();

        Dune::FieldVector < Scalar, dimWorld > distVec = globalPosInterface - globalPosI;
        Scalar lI = distVec * unitOuterNormal;
        distVec = globalPosJ - globalPosInterface;
        Scalar lJ = distVec * unitOuterNormal;
        Scalar l = lI + lJ;

        FieldMatrix permeabilityI(0);
        FieldMatrix permeabilityJ(0);
        FieldMatrix permeabilityK(0);

        problem_.spatialParams().meanK(permeabilityI,
                problem_.spatialParams().intrinsicPermeability(elementI));
        problem_.spatialParams().meanK(permeabilityJ,
                problem_.spatialParams().intrinsicPermeability(elementJ));
        problem_.spatialParams().meanK(permeabilityK,
                problem_.spatialParams().intrinsicPermeability(elementK));

        // Calculate permeablity component normal to interface
        Scalar kI, kJ, kK, kMean, ng;
        Dune::FieldVector < Scalar, dim > permI(0);
        Dune::FieldVector < Scalar, dim > permJ(0);
        Dune::FieldVector < Scalar, dim > permK(0);

        permeabilityI.mv(unitOuterNormal, permI);
        permeabilityJ.mv(unitOuterNormal, permJ);
        permeabilityK.mv(unitOuterNormal, permK);

        // kI,kJ,kK=(n^T)Kn
        kI = unitOuterNormal * permI;
        kJ = unitOuterNormal * permJ;
        kK = unitOuterNormal * permK;

        // See Diplomarbeit Michael Sinsbeck
        kMean = kI * kJ * kK * l / (kJ * kK * lI + kI * (kJ + kK) / 2 * lJ);

        ng = gravity_ * unitOuterNormal;

        Scalar fractionalWK = cellDataK.fracFlowFunc(wPhaseIdx);
        Scalar fractionalNwK = cellDataK.fracFlowFunc(nPhaseIdx);

        Scalar pcK = cellDataK.capillaryPressure();
        Scalar pcJK = (pcJ + pcK) / 2;

        // Potentials from previous time step are not available.
        // Instead, calculate mean density, then find potentials,
        // then upwind density.
        // pressure from previous time step might also be incorrect.

        Scalar rhoMeanWIJ = density_[wPhaseIdx];
        Scalar rhoMeanNwIJ = density_[nPhaseIdx];
        Scalar rhoMeanWIK = density_[wPhaseIdx];
        Scalar rhoMeanNwIK = density_[nPhaseIdx];

        if (compressibility_)
        {
            rhoMeanWIJ = (lJ * cellData.density(wPhaseIdx) + lI * cellDataJ.density(wPhaseIdx)) / l;
            rhoMeanNwIJ = (lJ * cellData.density(nPhaseIdx) + lI * cellDataJ.density(nPhaseIdx)) / l;
            rhoMeanWIK = (lJ * cellData.density(wPhaseIdx) + lI * cellDataK.density(wPhaseIdx)) / l;
            rhoMeanNwIK = (lJ * cellData.density(nPhaseIdx) + lI * cellDataK.density(nPhaseIdx)) / l;
        }

        Scalar densityWIJ = density_[wPhaseIdx];
        Scalar densityNwIJ = density_[nPhaseIdx];
        Scalar densityWIK = density_[wPhaseIdx];
        Scalar densityNwIK = density_[nPhaseIdx];

        Scalar fMeanWIJ = (lJ * fractionalWI + lI * fractionalWJ) / l;
        Scalar fMeanNwIJ = (lJ * fractionalNwI + lI * fractionalNwJ) / l;
        Scalar fMeanWIK = (lJ * fractionalWI + lI * fractionalWK) / l;
        Scalar fMeanNwIK = (lJ * fractionalNwI + lI * fractionalNwK) / l;

        Scalar potentialWIJ = 0;
        Scalar potentialNwIJ = 0;
        Scalar potentialWIK = 0;
        Scalar potentialNwIK = 0;

        // if we are at the very first iteration we can't calculate phase potentials
        if (!first)
        {
            // potentials from previous time step no available.
            if (compressibility_)
            {
                densityWIJ = rhoMeanWIJ;
                densityNwIJ = rhoMeanNwIJ;
                densityWIK = rhoMeanWIK;
                densityNwIK = rhoMeanNwIK;
            }

            switch (pressureType_)
            {
            case pGlobal:
            {
                Scalar pressJK = (cellDataJ.globalPressure() + cellDataK.globalPressure()) / 2;

                potentialWIJ = (cellData.globalPressure() - fMeanNwIJ * pcI - (pressJK - fMeanNwIJ * pcJK)) / l
                        + (densityWIJ - lJ / l * (kI + kK) / kI * (densityWIK - densityWIJ) / 2) * ng;
                potentialNwIJ = (cellData.globalPressure() + fMeanWIJ * pcI - (pressJK + fMeanWIJ * pcJK)) / l
                        + (densityNwIJ - lJ / l * (kI + kK) / kI * (densityNwIK - densityNwIJ) / 2) * ng;
                potentialWIK = (cellData.globalPressure() - fMeanNwIK * pcI - (pressJK - fMeanNwIK * pcJK)) / l
                        + (densityWIK - lJ / l * (kI + kK) / kI * (densityWIJ - densityWIK) / 2) * ng;
                potentialNwIK = (cellData.globalPressure() + fMeanWIK * pcI - (pressJK + fMeanWIK * pcJK)) / l
                        + (densityNwIK - lJ / l * (kI + kK) / kI * (densityNwIJ - densityNwIK) / 2) * ng;
                break;
            }
            default:
            {
                potentialWIJ = (cellData.pressure(wPhaseIdx)
                        - 0.5 * (cellDataJ.pressure(wPhaseIdx) + cellDataK.pressure(wPhaseIdx))) / l
                        + (densityWIJ - lJ / l * (kI + kK) / kI * (densityWIK - densityWIJ) / 2) * ng;
                potentialNwIJ = (cellData.pressure(nPhaseIdx)
                        - (0.5 * (cellDataJ.pressure(nPhaseIdx) + cellDataK.pressure(nPhaseIdx)))) / l
                        + (densityNwIJ - lJ / l * (kI + kK) / kI * (densityNwIK - densityNwIJ) / 2) * ng;
                potentialWIK = (cellData.pressure(wPhaseIdx)
                        - 0.5 * (cellDataJ.pressure(wPhaseIdx) + cellDataK.pressure(wPhaseIdx))) / l
                        + (densityWIK - lJ / l * (kI + kK) / kI * (densityWIJ - densityWIK) / 2) * ng;
                potentialNwIK = (cellData.pressure(nPhaseIdx)
                        - (0.5 * (cellDataJ.pressure(nPhaseIdx) + cellDataK.pressure(nPhaseIdx)))) / l
                        + (densityNwIK - lJ / l * (kI + kK) / kI * (densityNwIJ - densityNwIK) / 2) * ng;
                break;
            }
            }
        }

        // do the upwinding of the mobility depending on the phase potentials
        Scalar lambdaWIJ = (potentialWIJ > 0.) ? lambdaWI : lambdaWJ;
        lambdaWIJ = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialWIJ, 0.0, 1.0e-30)) ? 0.5 * (lambdaWI + lambdaWJ) : lambdaWIJ;
        Scalar lambdaNwIJ = (potentialNwIJ > 0) ? lambdaNwI : lambdaNwJ;
        lambdaNwIJ = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialNwIJ, 0.0, 1.0e-30)) ? 0.5 * (lambdaNwI + lambdaNwJ) : lambdaNwIJ;

        if (compressibility_)
        {
            densityWIJ = (potentialWIJ > 0.) ? cellData.density(wPhaseIdx) : cellDataJ.density(wPhaseIdx);
            densityNwIJ = (potentialNwIJ > 0.) ? cellData.density(nPhaseIdx) : cellDataJ.density(nPhaseIdx);
            densityWIJ = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialWIJ, 0.0, 1.0e-30)) ? rhoMeanWIJ : densityWIJ;
            densityNwIJ = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialNwIJ, 0.0, 1.0e-30)) ? rhoMeanNwIJ : densityNwIJ;
            densityWIK = (potentialWIK > 0.) ? cellData.density(wPhaseIdx) : cellDataK.density(wPhaseIdx);
            densityNwIK = (potentialNwIK > 0.) ? cellData.density(nPhaseIdx) : densityNwIK;
            densityWIK = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialWIK, 0.0, 1.0e-30)) ? rhoMeanWIK : densityWIK;
            densityNwIK = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialNwIK, 0.0, 1.0e-30)) ? rhoMeanNwIK : densityNwIK;
        }

        // update diagonal entry and right hand side
        entry[matrix] = (lambdaWIJ + lambdaNwIJ) * kMean / l * faceArea;
        entry[rhs] = faceArea * lambdaWIJ * kMean * ng
                * ((densityWIJ) - (lJ / l) * (kI + kK) / kI * (densityWIK - densityWIJ) / 2);
        entry[rhs] += faceArea * lambdaNwIJ * kMean * ng
                * ((densityNwIJ) - (lJ / l) * (kI + kK) / kI * (densityNwIK - densityNwIJ) / 2);

        switch (pressureType_)
        {
        case pw:
        {
            entry[rhs] += faceArea * lambdaNwIJ * kMean * (pcJK - pcI) / l;
            break;
        }
        case pn:
        {
            entry[rhs] -= faceArea * lambdaWIJ * kMean * (pcJK - pcI) / l;
            break;
        }
        }

        // write hanging-node-specific stuff directly into matrix and rhs!
        this->f_[globalIdxI] -= entry[rhs];
        this->f_[globalIdxJ] += entry[rhs];

        // set diagonal entry
        this->A_[globalIdxI][globalIdxI] += entry[matrix];
        // set off-diagonal
        this->A_[globalIdxI][globalIdxJ] -= entry[matrix];

        // set entry for cell J
        this->A_[globalIdxJ][globalIdxI] -= entry[matrix];
        this->A_[globalIdxJ][globalIdxJ] += entry[matrix];

        // set entry to zero -> matrix already written!
        entry = 0.;

//        std::cout<<"finished hanging node!\n";
    }
    else
    {
        entry = 0;
    }

    return;
}

} // end namespace Dumux
#endif
