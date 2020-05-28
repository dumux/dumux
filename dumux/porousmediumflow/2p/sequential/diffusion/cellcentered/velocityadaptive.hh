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
 * \brief  Velocity field from a finite volume solution of a pressure equation.
 */
#ifndef DUMUX_FVVELOCITY2P_ADAPTIVE_HH
#define DUMUX_FVVELOCITY2P_ADAPTIVE_HH

#include <dune/common/float_cmp.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/velocity.hh>

namespace Dumux {

/*!
 * \ingroup SequentialTwoPModel
 * \brief Determines the velocity from a finite volume solution of the  pressure equation of a sequential model (IMPES).
 *
 * Details see FVVelocity2P
 */
template<class TypeTag>
class FVVelocity2PAdaptive: public FVVelocity2P<TypeTag>
{
    using ParentType = FVVelocity2P<TypeTag>;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
     using Scalar = GetPropType<TypeTag, Properties::Scalar>;
     using Problem = GetPropType<TypeTag, Properties::Problem>;


     using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

     using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
     using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    using CellData = GetPropType<TypeTag, Properties::CellData>;

    using Intersection = typename GridView::Intersection;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNw,
        pGlobal = Indices::pressureGlobal,
        vw = Indices::velocityW,
        vn = Indices::velocityNw,
        vt = Indices::velocityTotal
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx, numPhases = getPropValue<TypeTag, Properties::NumPhases>()
    };

    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GravityVector = Dune::FieldVector<Scalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

public:
    /*!
     * \brief Constructs a FVVelocity2PAdaptive object
     * \param problem A problem class object
     */
    FVVelocity2PAdaptive(Problem& problem)
    : ParentType(problem), problem_(problem), gravity_(problem.gravity())
    {
        if (getPropValue<TypeTag, Properties::EnableCompressibility>() && velocityType_ == vt)
        {
            DUNE_THROW(Dune::NotImplemented, "Total velocity - global pressure - model cannot be used with compressible fluids!");
        }
        if (velocityType_ != vw && velocityType_ != vn && velocityType_ != vt)
        {
            DUNE_THROW(Dune::NotImplemented, "Velocity type not supported!");
        }

        density_[wPhaseIdx] = 0.;
        density_[nPhaseIdx] = 0.;
        viscosity_[wPhaseIdx] = 0.;
        viscosity_[nPhaseIdx] = 0.;
    }

    //! For initialization
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
    }

    //! Calculates the velocity at a cell-cell interface.
    void calculateVelocity(const Intersection& intersection, CellData& cellData);

    /*!
     * \brief Indicates if velocity is reconstructed in the pressure step or in the transport step
     *
     * Returns false (In the adaptive finite volume scheme the velocity has to be calculated separately
     * to make sure the hanging nodes are treated correctly.)
     */
    bool calculateVelocityInTransport()
    {
            return false;
    }

private:
    Problem& problem_;
    const GravityVector& gravity_; //!< vector including the gravity constant
    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];

    //! Gives kind of velocity used (\f$ 0 = v_w\f$, \f$ 1 = v_n\f$, \f$ 2 = v_t\f$)
    static const int velocityType_ = getPropValue<TypeTag, Properties::VelocityFormulation>();
    static const bool compressibility_ = getPropValue<TypeTag, Properties::EnableCompressibility>();
    //! Gives kind of pressure used (\f$p_w\f$, \f$p_n\f$, \f$p_{global}\f$)
    static const int pressureType_ = getPropValue<TypeTag, Properties::PressureFormulation>();
    //! Gives kind of saturation used (\f$S_w\f$, \f$S_n\f$)
    static const int saturationType_ = getPropValue<TypeTag, Properties::SaturationFormulation>();
};

/*!
 * \brief Calculates the velocity at a cell-cell interface.
 *
 * \copydetails FVVelocity2P::calculateVelocity(const Intersection&,CellData&)
 *
 * Implementation of calculateVelocity() function for cell-cell interfaces with hanging nodes.
 */
template<class TypeTag>
void FVVelocity2PAdaptive<TypeTag>::calculateVelocity(const Intersection& intersection, CellData& cellData)
{
    auto elementI = intersection.inside();
    auto elementJ = intersection.outside();

    if (elementI.level() == elementJ.level())
    {
        ParentType::calculateVelocity(intersection, cellData);
    }
    else if (elementJ.level() == elementI.level() + 1 && dim == 2)
    {
        CellData& cellDataJ = problem_.variables().cellData(problem_.variables().index(elementJ));

        // get global coordinates of cell centers
        const GlobalPosition& globalPosI = elementI.geometry().center();
        const GlobalPosition& globalPosJ = elementJ.geometry().center();

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

        Scalar faceArea = intersection.geometry().volume();
        Scalar faceAreaSum = faceArea;

        // get face normal
        const Dune::FieldVector<Scalar, dim>& unitOuterNormal = intersection.centerUnitOuterNormal();

        // Count number of hanging nodes
        // not really necessary
        int isIndexJ = intersection.indexInOutside();

        int globalIdxK = 0;
        auto elementK = intersection.outside();
        // We are looking for two things:
        // IsIndexJ, the index of the interface from the neighbor-cell point of view
        // GlobalIdxK, the index of the third cell
        // for efficiency this is done in one intersection loop

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
                    faceAreaSum += intersectionI.geometry().volume();

                    break;
                }
            }
        }

        CellData& cellDataK = problem_.variables().cellData(globalIdxK);

        // face global coordinates
        const GlobalPosition& globalPosInterface = intersection.geometry().center();

        Dune::FieldVector < Scalar, dimWorld > distVec = globalPosInterface - globalPosI;
        Scalar lI = distVec * unitOuterNormal;
        distVec = globalPosJ - globalPosInterface;
        Scalar lJ = distVec * unitOuterNormal;
        Scalar l = lI + lJ;

        DimMatrix permeabilityI(0);
        DimMatrix permeabilityJ(0);
        DimMatrix permeabilityK(0);

        problem_.spatialParams().meanK(permeabilityI,
                problem_.spatialParams().intrinsicPermeability(elementI));
        problem_.spatialParams().meanK(permeabilityJ,
                problem_.spatialParams().intrinsicPermeability(elementJ));
        problem_.spatialParams().meanK(permeabilityK,
                problem_.spatialParams().intrinsicPermeability(elementK));

        // Calculate permeablity component normal to interface
        Scalar kI, kJ, kK, ng, kMean; //, gI, gJ, gK;
        Dune::FieldVector < Scalar, dim > permI(0);
        Dune::FieldVector < Scalar, dim > permJ(0);
        Dune::FieldVector < Scalar, dim > permK(0);

        permeabilityI.mv(unitOuterNormal, permI);
        permeabilityJ.mv(unitOuterNormal, permJ);
        permeabilityK.mv(unitOuterNormal, permK);

        // kI,kJ,kK = (n^T)Kn
        kI = unitOuterNormal * permI;
        kJ = unitOuterNormal * permJ;
        kK = unitOuterNormal * permK;
        kMean = kI * kJ * kK * l / (kJ * kK * lI + kI * (kJ + kK) / 2 * lJ);

        ng = gravity_ * unitOuterNormal;

        Scalar fractionalWK = cellDataK.fracFlowFunc(wPhaseIdx);
        Scalar fractionalNwK = cellDataK.fracFlowFunc(nPhaseIdx);

        Scalar pcK = cellDataK.capillaryPressure();
        Scalar pcJK = (pcJ + pcK) / 2;

        // calculate potential gradients
        // reuse potentials from fvpressure2padaptive
        Scalar potentialDiffWIJ = cellData.fluxData().upwindPotential(wPhaseIdx, isIndexI);
        Scalar potentialDiffNwIJ = cellData.fluxData().upwindPotential(nPhaseIdx, isIndexI);
        Scalar potentialDiffWIK = potentialDiffWIJ;
        Scalar potentialDiffNwIK = potentialDiffNwIJ;
        // preliminary potential. The "real" ones are found below

        // Comment: reversed weighting is plausible, too (swap lJ and lI)
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

        Scalar fMeanWIJ = (lJ * fractionalWI + lI * fractionalWJ) / l;
        Scalar fMeanNwIJ = (lJ * fractionalNwI + lI * fractionalNwJ) / l;
        Scalar fMeanWIK = (lJ * fractionalWI + lI * fractionalWK) / l;
        Scalar fMeanNwIK = (lJ * fractionalNwI + lI * fractionalNwK) / l;

        Scalar densityWIJ = density_[wPhaseIdx];
        Scalar densityNwIJ = density_[nPhaseIdx];
        Scalar densityWIK = density_[wPhaseIdx];
        Scalar densityNwIK = density_[nPhaseIdx];

        if (compressibility_)
        {
        // Upwinding for finding the upwinding direction
            densityWIJ = (potentialDiffWIJ > 0.) ? cellData.density(wPhaseIdx) : cellDataJ.density(wPhaseIdx);
            densityNwIJ = (potentialDiffNwIJ > 0.) ? cellData.density(nPhaseIdx) : cellDataJ.density(nPhaseIdx);
            densityWIK = (potentialDiffWIK > 0.) ? cellData.density(wPhaseIdx) : cellDataK.density(wPhaseIdx);
            densityNwIK = (potentialDiffNwIK > 0.) ? cellData.density(nPhaseIdx) : cellDataK.density(nPhaseIdx);

            densityWIJ = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffWIJ, 0.0, 1.0e-30)) ? rhoMeanWIJ : densityWIJ;
            densityNwIJ = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNwIJ, 0.0, 1.0e-30)) ? rhoMeanNwIJ : densityNwIJ;
            densityWIK = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffWIK, 0.0, 1.0e-30)) ? rhoMeanWIK : densityWIK;
            densityNwIK = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNwIK, 0.0, 1.0e-30)) ? rhoMeanNwIK : densityNwIK;
        }

        Scalar fractionalWIJ = (potentialDiffWIJ > 0.) ? fractionalWI : fractionalWJ;
        Scalar fractionalNwIJ = (potentialDiffNwIJ > 0.) ? fractionalNwI : fractionalNwJ;
        Scalar fractionalWIK = (potentialDiffWIK > 0.) ? fractionalWI : fractionalWK;
        Scalar fractionalNwIK = (potentialDiffNwIK > 0.) ? fractionalNwI : fractionalNwK;

        fractionalWIJ = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffWIJ, 0.0, 1.0e-30)) ? fMeanWIJ : fractionalWIJ;
        fractionalNwIJ = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNwIJ, 0.0, 1.0e-30)) ? fMeanNwIJ : fractionalNwIJ;
        fractionalWIK = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffWIK, 0.0, 1.0e-30)) ? fMeanWIK : fractionalWIK;
        fractionalNwIK = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNwIK, 0.0, 1.0e-30)) ? fMeanNwIK : fractionalNwIK;

        switch (pressureType_)
        {
        case pGlobal:
        {
            Scalar pressJK = (cellDataJ.globalPressure() + cellDataK.globalPressure()) / 2;

            potentialDiffWIJ = (cellData.globalPressure() - fMeanNwIJ * pcI - (pressJK - fMeanNwIJ * pcJK)) / l
                    + (densityWIJ - lJ / l * (kI + kK) / kI * (densityWIK - densityWIJ) / 2) * ng;
            potentialDiffNwIJ = (cellData.globalPressure() + fMeanWIJ * pcI - (pressJK + fMeanWIJ * pcJK)) / l
                    + (densityNwIJ - lJ / l * (kI + kK) / kI * (densityNwIK - densityNwIJ) / 2) * ng;
            potentialDiffWIK = (cellData.globalPressure() - fMeanNwIK * pcI - (pressJK - fMeanNwIK * pcJK)) / l
                    + (densityWIK - lJ / l * (kI + kK) / kI * (densityWIJ - densityWIK) / 2) * ng;
            potentialDiffNwIK = (cellData.globalPressure() + fMeanWIK * pcI - (pressJK + fMeanWIK * pcJK)) / l
                    + (densityNwIK - lJ / l * (kI + kK) / kI * (densityNwIJ - densityNwIK) / 2) * ng;
            break;
        }
        default:
        {
            potentialDiffWIJ = (cellData.pressure(wPhaseIdx)
                    - 0.5 * (cellDataJ.pressure(wPhaseIdx) + cellDataK.pressure(wPhaseIdx))) / l
                    + (densityWIJ - lJ / l * (kI + kK) / kI * (densityWIK - densityWIJ) / 2) * ng;
            potentialDiffNwIJ = (cellData.pressure(nPhaseIdx)
                    - (0.5 * (cellDataJ.pressure(nPhaseIdx) + cellDataK.pressure(nPhaseIdx)))) / l
                    + (densityNwIJ - lJ / l * (kI + kK) / kI * (densityNwIK - densityNwIJ) / 2) * ng;
            potentialDiffWIK = (cellData.pressure(wPhaseIdx)
                    - 0.5 * (cellDataJ.pressure(wPhaseIdx) + cellDataK.pressure(wPhaseIdx))) / l
                    + (densityWIK - lJ / l * (kI + kK) / kI * (densityWIJ - densityWIK) / 2) * ng;
            potentialDiffNwIK = (cellData.pressure(nPhaseIdx)
                    - (0.5 * (cellDataJ.pressure(nPhaseIdx) + cellDataK.pressure(nPhaseIdx)))) / l
                    + (densityNwIK - lJ / l * (kI + kK) / kI * (densityNwIJ - densityNwIK) / 2) * ng;
            break;
        }
        }

        // store potentials for further calculations (velocity, saturation, ...)
        // these quantities only have correct sign (needed for upwinding)
        // potentials are defined slightly different for adaptive scheme
        cellData.fluxData().addUpwindPotential(wPhaseIdx, isIndexI, potentialDiffWIJ);
        cellData.fluxData().addUpwindPotential(nPhaseIdx, isIndexI, potentialDiffNwIJ);
        cellDataJ.fluxData().setUpwindPotential(wPhaseIdx, isIndexJ, -potentialDiffWIJ);
        cellDataJ.fluxData().setUpwindPotential(nPhaseIdx, isIndexJ, -potentialDiffNwIJ);

        //do the upwinding of the mobility depending on the phase potentials
        Scalar lambdaWIJ = (potentialDiffWIJ > 0.) ? lambdaWI : lambdaWJ;
        lambdaWIJ = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffWIJ, 0.0, 1.0e-30)) ? 0.5 * (lambdaWI + lambdaWJ) : lambdaWIJ;
        Scalar lambdaNwIJ = (potentialDiffNwIJ > 0.) ? lambdaNwI : lambdaNwJ;
        lambdaNwIJ = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNwIJ, 0.0, 1.0e-30)) ? 0.5 * (lambdaNwI + lambdaNwJ) : lambdaNwIJ;

        if (compressibility_)
        {
            densityWIJ = (potentialDiffWIJ > 0.) ? cellData.density(wPhaseIdx) : cellDataJ.density(wPhaseIdx);
            densityNwIJ = (potentialDiffNwIJ > 0.) ? cellData.density(nPhaseIdx) : cellDataJ.density(nPhaseIdx);
            densityWIJ = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffWIJ, 0.0, 1.0e-30)) ? rhoMeanWIJ : densityWIJ;
            densityNwIJ = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNwIJ, 0.0, 1.0e-30)) ? rhoMeanNwIJ : densityNwIJ;
            densityWIK = (potentialDiffWIK > 0.) ? cellData.density(wPhaseIdx) : cellDataK.density(wPhaseIdx);
            densityNwIK = (potentialDiffNwIK > 0.) ? cellData.density(nPhaseIdx) : cellDataK.density(nPhaseIdx);
            densityWIK = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffWIK, 0.0, 1.0e-30)) ? rhoMeanWIK : densityWIK;
            densityNwIK = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNwIK, 0.0, 1.0e-30)) ? rhoMeanNwIK : densityNwIK;
        }

        // calculate velocities and the gravity term
        Dune::FieldVector < Scalar, dimWorld > velocityW(unitOuterNormal);
        Dune::FieldVector < Scalar, dimWorld > velocityNw(unitOuterNormal);
        Dune::FieldVector < Scalar, dimWorld > gravityTermW(unitOuterNormal);
        Dune::FieldVector < Scalar, dimWorld > gravityTermNw(unitOuterNormal);

        gravityTermW *= lambdaWIJ * kMean * ng;
        gravityTermW *= densityWIJ - (lJ / l) * (kI + kK) / kI * (densityWIK - densityWIJ) / 2;
        gravityTermNw *= lambdaNwIJ * kMean * ng;
        gravityTermNw *= densityNwIJ - (lJ / l) * (kI + kK) / kI * (densityNwIK - densityNwIJ) / 2;

        switch (pressureType_)
        {
        case pw:
        {
            velocityW *= lambdaWIJ * kMean * (cellData.pressure(wPhaseIdx) -
                          (cellDataJ.pressure(wPhaseIdx) + cellDataK.pressure(wPhaseIdx)) / 2.0) / l;
            velocityNw *= lambdaNwIJ * kMean * (cellData.pressure(nPhaseIdx) -
                          (cellDataJ.pressure(nPhaseIdx) + cellDataK.pressure(nPhaseIdx)) / 2.0) / l;

            velocityW += gravityTermW;
            velocityNw += gravityTermNw;
            break;
        }
        case pn:
        {
            velocityNw *= lambdaNwIJ * kMean * (cellData.pressure(nPhaseIdx) -
                           (cellDataJ.pressure(nPhaseIdx) + cellDataK.pressure(nPhaseIdx)) / 2.0) / l;
            velocityW *= lambdaWIJ * kMean * (cellData.pressure(wPhaseIdx) -
                           (cellDataJ.pressure(wPhaseIdx) + cellDataK.pressure(wPhaseIdx)) / 2.0) / l;

            velocityW += gravityTermW;
            velocityNw += gravityTermNw;
            break;
        }
        case pGlobal:
        {
            Scalar pressJK = (cellDataJ.globalPressure() + cellDataK.globalPressure()) / 2;

            velocityW *= lambdaWIJ * kMean * (cellData.globalPressure() - pressJK) / l;
            velocityW += gravityTermW;
            velocityW += gravityTermNw;
            velocityNw = 0;
            break;
        }
        }

        cellDataJ.fluxData().setVelocity(wPhaseIdx, isIndexJ, velocityW);
        cellDataJ.fluxData().setVelocity(nPhaseIdx, isIndexJ, velocityNw);
        cellDataJ.fluxData().setVelocityMarker(isIndexJ);

        //times 0.5 because cell face with hanging node is called twice! Do not set marker because it should be called twice!
        velocityW *= faceArea/faceAreaSum;
        velocityNw *= faceArea/faceAreaSum;
        cellData.fluxData().addVelocity(wPhaseIdx, isIndexI, velocityW);
        cellData.fluxData().addVelocity(nPhaseIdx, isIndexI, velocityNw);
    }
    else if (elementI.level() > elementJ.level() && dim == 3)
    {
        int globalIdxJ = problem_.variables().index(elementJ);

        CellData& cellDataJ = problem_.variables().cellData(globalIdxJ);

        // get global coordinates of cell centers
        const GlobalPosition& globalPosI = elementI.geometry().center();
        const GlobalPosition& globalPosJ = elementJ.geometry().center();

        // get mobilities and fractional flow factors
        Scalar lambdaWI = cellData.mobility(wPhaseIdx);
        Scalar lambdaNwI = cellData.mobility(nPhaseIdx);
        Scalar lambdaWJ = cellDataJ.mobility(wPhaseIdx);
        Scalar lambdaNwJ = cellDataJ.mobility(nPhaseIdx);

        // get capillary pressure
        Scalar pcI = cellData.capillaryPressure();
        Scalar pcJ = cellDataJ.capillaryPressure();

        // get face index
        int isIndexI = intersection.indexInInside();
        int isIndexJ = intersection.indexInOutside();

        // get face normal
        const Dune::FieldVector<Scalar, dim>& unitOuterNormal = intersection.centerUnitOuterNormal();

        // distance vector between barycenters
        GlobalPosition distVec = globalPosJ - globalPosI;

        // compute distance between cell centers
        Scalar dist = distVec.two_norm();

        // compute vectorized permeabilities
        DimMatrix meanPermeability(0);

        problem_.spatialParams().meanK(meanPermeability, problem_.spatialParams().intrinsicPermeability(elementI),
                problem_.spatialParams().intrinsicPermeability(elementJ));

        Dune::FieldVector<Scalar, dim> permeability(0);
        meanPermeability.mv(unitOuterNormal, permeability);

        // calculate potential gradients
        Scalar potentialDiffW = cellData.potential(wPhaseIdx) - cellDataJ.potential(wPhaseIdx);
        Scalar potentialDiffNw = cellData.potential(nPhaseIdx) - cellDataJ.potential(nPhaseIdx);

        if (compressibility_)
        {
            density_[wPhaseIdx] = (potentialDiffW > 0.) ? cellData.density(wPhaseIdx) : cellDataJ.density(wPhaseIdx);
            density_[nPhaseIdx] = (potentialDiffNw > 0.) ? cellData.density(nPhaseIdx) : cellDataJ.density(nPhaseIdx);

            density_[wPhaseIdx] =
                    (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffW, 0.0, 1.0e-30)) ? 0.5 * (cellData.density(wPhaseIdx) + cellDataJ.density(wPhaseIdx)) :
                            density_[wPhaseIdx];
            density_[nPhaseIdx] =
                    (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNw, 0.0, 1.0e-30)) ? 0.5 * (cellData.density(nPhaseIdx) + cellDataJ.density(nPhaseIdx)) :
                            density_[nPhaseIdx];

            potentialDiffW = (cellData.pressure(wPhaseIdx) - cellDataJ.pressure(wPhaseIdx));
            potentialDiffNw = (cellData.pressure(nPhaseIdx) - cellDataJ.pressure(nPhaseIdx));


            potentialDiffW += density_[wPhaseIdx] * (distVec * gravity_); //delta z/delta x in unitOuterNormal[z]
            potentialDiffNw += density_[nPhaseIdx] * (distVec * gravity_);
        }

        // store potentials for further calculations (velocity, saturation, ...)
        cellData.fluxData().setUpwindPotential(wPhaseIdx, isIndexI, potentialDiffW);
        cellData.fluxData().setUpwindPotential(nPhaseIdx, isIndexI, potentialDiffNw);

        cellDataJ.fluxData().setUpwindPotential(wPhaseIdx, isIndexJ, -potentialDiffW);
        cellDataJ.fluxData().setUpwindPotential(nPhaseIdx, isIndexJ, -potentialDiffNw);

        // do the upwinding of the mobility depending on the phase potentials
        Scalar lambdaW = (potentialDiffW > 0.) ? lambdaWI : lambdaWJ;
        lambdaW = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffW, 0.0, 1.0e-30 )) ? 0.5 * (lambdaWI + lambdaWJ) : lambdaW;
        Scalar lambdaNw = (potentialDiffNw > 0.) ? lambdaNwI : lambdaNwJ;
        lambdaNw = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNw, 0.0, 1.0e-30)) ? 0.5 * (lambdaNwI + lambdaNwJ) : lambdaNw;

        if (compressibility_)
        {
            density_[wPhaseIdx] = (potentialDiffW > 0.) ? cellData.density(wPhaseIdx) : cellDataJ.density(wPhaseIdx);
            density_[nPhaseIdx] = (potentialDiffNw > 0.) ? cellData.density(nPhaseIdx) : cellDataJ.density(nPhaseIdx);

            density_[wPhaseIdx] =
                    (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffW, 0.0, 1.0e-30)) ? 0.5 * (cellData.density(wPhaseIdx) + cellDataJ.density(wPhaseIdx)) :
                            density_[wPhaseIdx];
            density_[nPhaseIdx] =
                    (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNw, 0.0, 1.0e-30)) ? 0.5 * (cellData.density(nPhaseIdx) + cellDataJ.density(nPhaseIdx)) :
                            density_[nPhaseIdx];
        }

        Scalar scalarPerm = permeability.two_norm();

        // calculate the gravity term
        Dune::FieldVector<Scalar, dimWorld> velocityW(unitOuterNormal);
        Dune::FieldVector<Scalar, dimWorld> velocityNw(unitOuterNormal);

        // calculate unit distVec
        distVec /= dist;
        Scalar areaScaling = (unitOuterNormal * distVec);
        // this treatment of g allows to account for gravity flux through faces
        // where the face normal has no z component (e.g. parallelepiped grids)
        Scalar gravityTermW = (gravity_ * distVec) * density_[wPhaseIdx] * areaScaling;
        Scalar gravityTermNw = (gravity_ * distVec) * density_[nPhaseIdx] * areaScaling;

        // calculate velocity depending on the pressure used -> use pc = pn - pw
        switch (pressureType_)
        {
        case pw:
        {
            velocityW *= lambdaW * scalarPerm
                    * ((cellData.pressure(wPhaseIdx) - cellDataJ.pressure(wPhaseIdx)) / dist + gravityTermW);
            velocityNw *= lambdaNw * scalarPerm
                    * ((cellData.pressure(wPhaseIdx) - cellDataJ.pressure(wPhaseIdx)) / dist + gravityTermNw)
                    + 0.5 * (lambdaNwI + lambdaNwJ) * scalarPerm * (pcI - pcJ) / dist;
            break;
        }
        case pn:
        {
            velocityW *= lambdaW * scalarPerm
                    * ((cellData.pressure(nPhaseIdx) - cellDataJ.pressure(nPhaseIdx)) / dist + gravityTermW)
                    - 0.5 * (lambdaWI + lambdaWJ) * scalarPerm * (pcI - pcJ) / dist;
            velocityNw *= lambdaNw * scalarPerm
                    * ((cellData.pressure(nPhaseIdx) - cellDataJ.pressure(nPhaseIdx)) / dist + gravityTermNw);
            break;
        }
        case pGlobal:
        {
            velocityW *= (lambdaW + lambdaNw) * scalarPerm * (cellData.globalPressure() - cellDataJ.globalPressure()) / dist
                    + scalarPerm * (lambdaW * gravityTermW + lambdaNw * gravityTermNw);
            velocityNw = 0;
            break;
        }
        }

        // store velocities
        cellData.fluxData().setVelocity(wPhaseIdx, isIndexI, velocityW);
        cellData.fluxData().setVelocity(nPhaseIdx, isIndexI, velocityNw);
        cellData.fluxData().setVelocityMarker(isIndexI);

        using std::pow;
        Scalar weightingFactor = pow(0.5, (dim - 1)*(elementI.level() - elementJ.level()));

        velocityW *= weightingFactor;
        velocityNw *= weightingFactor;

        cellDataJ.fluxData().addVelocity(wPhaseIdx, isIndexJ, velocityW);
        cellDataJ.fluxData().addVelocity(nPhaseIdx, isIndexJ, velocityNw);
        cellDataJ.fluxData().setVelocityMarker(isIndexJ);
    }

    return;
}
} // end namespace Dumux
#endif
