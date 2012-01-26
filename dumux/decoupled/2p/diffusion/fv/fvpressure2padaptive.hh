// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
#ifndef DUMUX_FVPRESSURE2P_ADAPTIVE_HH
#define DUMUX_FVPRESSURE2P_ADAPTIVE_HH

// dumux environment
#include <dumux/decoupled/2p/2pproperties.hh>
#include <dumux/decoupled/2p/diffusion/fv/fvpressure2p.hh>
#include <dumux/decoupled/common/fv/fvvelocity.hh>

/**
 * @file
 * @brief  Finite Volume Discretization of a pressure equation.
 * @author Bernd Flemisch, Jochen Fritz, Markus Wolff
 */

namespace Dumux
{

//! \ingroup FV2p
//! \brief Finite Volume discretization of the pressure equation of the sequential IMPES Model.
/*! Provides a Finite Volume implementation for the evaluation
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
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag> class FVPressure2PAdaptive: public FVPressure2P<TypeTag>
{
    typedef FVPressure2P<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, TwoPIndices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;
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
        Sn = Indices::saturationNW,
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        eqIdxPress = Indices::pressEqIdx,
        eqIdxSat = Indices::satEqIdx
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx, numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    enum
    {
        rhs = ParentType::rhs, matrix = ParentType::matrix
    };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::Intersection Intersection;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;

public:
    void getFlux(Dune::FieldVector<Scalar, 2>&, const Intersection&, const CellData&, const bool);

    //pressure solution routine: update estimate for secants, assemble, solve.
    void update()
    {
        int gridSize = problem_.gridView().size(0);
        // update RHS vector, matrix
        this->A_.setSize(gridSize, gridSize); //
        this->f_.resize(gridSize);
        this->pressure().resize(gridSize);

        ParentType::initializeMatrix();

        ParentType::update();

        velocity_.calculateVelocity();

        return;
    }

    //! \copydoc Dumux::FVPressure1P::addOutputVtkFields(MultiWriter &writer)
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        ParentType::addOutputVtkFields(writer);
        velocity_.addOutputVtkFields(writer);
    }

    //! Constructs a FVPressure2PAdaptive object
    /**
     * \param problem a problem class object
     */
    FVPressure2PAdaptive(Problem& problem) :
            ParentType(problem), problem_(problem), velocity_(problem), gravity_(problem.gravity())
    {
        if (pressureType_ != pw && pressureType_ != pn && pressureType_ != pglobal)
        {
            DUNE_THROW(Dune::NotImplemented, "Pressure type not supported!");
        }
        if (pressureType_ == pglobal && compressibility_)
        {
            DUNE_THROW(Dune::NotImplemented, "Compressibility not supported for global pressure!");
        }
        if (saturationType_ != Sw && saturationType_ != Sn)
        {
            DUNE_THROW(Dune::NotImplemented, "Saturation type not supported!");
        }

        if (!compressibility_)
        {
            const Element& element = *(problem_.gridView().template begin<0>());
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

private:
    Problem& problem_;
    FVVelocity<TypeTag, typename GET_PROP_TYPE(TypeTag, Velocity)> velocity_;
    const GlobalPosition& gravity_; //!< vector including the gravity constant

    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];

    static const bool compressibility_ = GET_PROP_VALUE(TypeTag, EnableCompressibility);
    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PressureFormulation); //!< gives kind of pressure used (\f$p_w\f$, \f$p_n\f$, \f$p_{global}\f$)
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, SaturationFormulation); //!< gives kind of saturation used (\f$S_w\f$, \f$S_n\f$)
};

//!function which assembles the system of equations to be solved. Only works for 2D.
template<class TypeTag>
void FVPressure2PAdaptive<TypeTag>::getFlux(Dune::FieldVector<Scalar, 2>& entries, const Intersection& intersection
        , const CellData& cellDataI, const bool first)
{
    ElementPointer elementI = intersection.inside();
    ElementPointer elementJ = intersection.outside();

    if (elementI->level() == elementJ->level())
    {
        ParentType::getFlux(entries, intersection, cellDataI, first);
    }
    // hanging node situation: neighbor has higher level
    else if (elementJ->level() == elementI->level() + 1)
    {
        const CellData& cellDataJ = problem_.variables().cellData(problem_.variables().index(*elementJ));

        // get global coordinates of cell centers
        const GlobalPosition& globalPosI = elementI->geometry().center();
        const GlobalPosition& globalPosJ = elementJ->geometry().center();

        int globalIdxI = problem_.variables().index(*elementI);
        int globalIdxJ = problem_.variables().index(*elementJ);

        // get mobilities and fractional flow factors
        Scalar lambdaWI = cellDataI.mobility(wPhaseIdx);
        Scalar lambdaNWI = cellDataI.mobility(nPhaseIdx);
        Scalar fractionalWI = cellDataI.fracFlowFunc(wPhaseIdx);
        Scalar fractionalNWI = cellDataI.fracFlowFunc(nPhaseIdx);
        Scalar lambdaWJ = cellDataJ.mobility(wPhaseIdx);
        Scalar lambdaNWJ = cellDataJ.mobility(nPhaseIdx);
        Scalar fractionalWJ = cellDataJ.fracFlowFunc(wPhaseIdx);
        Scalar fractionalNWJ = cellDataJ.fracFlowFunc(nPhaseIdx);

        // get capillary pressure
        Scalar pcI = cellDataI.capillaryPressure();
        Scalar pcJ = cellDataJ.capillaryPressure();

        //get face index
        int isIndexI = intersection.indexInInside();

        //get face normal
        const Dune::FieldVector<Scalar, dim>& unitOuterNormal = intersection.centerUnitOuterNormal();

        // get face area
        Scalar faceArea = intersection.geometry().volume();

        // Count number of hanging nodes
        // not really necessary
        int globalIdxK = 0;
        ElementPointer elementK = intersection.outside();
        // We are looking for two things:
        // IsIndexJ, the index of the interface from the neighbor-cell point of view
        // GlobalIdxK, the index of the third cell
        // Intersectioniterator around cell I
        IntersectionIterator isItEndI = problem_.gridView().iend(*elementI);
        for (IntersectionIterator isItI = problem_.gridView().ibegin(*elementI); isItI != isItEndI; ++isItI)
        {
            if (isItI->neighbor())
            {
                ElementPointer neighborPointer2 = isItI->outside();

                // make sure we do not chose elemntI as third element
                // -> faces with hanging node have more than one intersection but only one face index!
                if (neighborPointer2 != elementJ && isItI->indexInInside() == isIndexI)
                {
                    globalIdxK = problem_.variables().index(*neighborPointer2);
                    elementK = neighborPointer2;
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

        problem_.spatialParameters().meanK(permeabilityI,
                problem_.spatialParameters().intrinsicPermeability(*elementI));
        problem_.spatialParameters().meanK(permeabilityJ,
                problem_.spatialParameters().intrinsicPermeability(*elementJ));
        problem_.spatialParameters().meanK(permeabilityK,
                problem_.spatialParameters().intrinsicPermeability(*elementK));

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
        Scalar fractionalNWK = cellDataK.fracFlowFunc(nPhaseIdx);

        Scalar pcK = cellDataK.capillaryPressure();
        Scalar pcJK = (pcJ + pcK) / 2;

        // Potentials from previous time step are not available.
        // Instead, calculate mean density, then find potentials,
        // then upwind density.
        // pressure from previous time step might also be incorrect.

        Scalar rhoMeanWIJ = density_[wPhaseIdx];
        Scalar rhoMeanNWIJ = density_[nPhaseIdx];
        Scalar rhoMeanWIK = density_[wPhaseIdx];
        Scalar rhoMeanNWIK = density_[nPhaseIdx];

        if (compressibility_)
        {
            rhoMeanWIJ = (lJ * cellDataI.density(wPhaseIdx) + lI * cellDataJ.density(wPhaseIdx)) / l;
            rhoMeanNWIJ = (lJ * cellDataI.density(nPhaseIdx) + lI * cellDataJ.density(nPhaseIdx)) / l;
            rhoMeanWIK = (lJ * cellDataI.density(wPhaseIdx) + lI * cellDataK.density(wPhaseIdx)) / l;
            rhoMeanNWIK = (lJ * cellDataI.density(nPhaseIdx) + lI * cellDataK.density(nPhaseIdx)) / l;
        }

        Scalar densityWIJ = density_[wPhaseIdx];
        Scalar densityNWIJ = density_[nPhaseIdx];
        Scalar densityWIK = density_[wPhaseIdx];
        Scalar densityNWIK = density_[nPhaseIdx];

        Scalar fMeanWIJ = (lJ * fractionalWI + lI * fractionalWJ) / l;
        Scalar fMeanNWIJ = (lJ * fractionalNWI + lI * fractionalNWJ) / l;
        Scalar fMeanWIK = (lJ * fractionalWI + lI * fractionalWK) / l;
        Scalar fMeanNWIK = (lJ * fractionalNWI + lI * fractionalNWK) / l;

        Scalar potentialWIJ = 0;
        Scalar potentialNWIJ = 0;
        Scalar potentialWIK = 0;
        Scalar potentialNWIK = 0;

        //if we are at the very first iteration we can't calculate phase potentials
        if (!first)
        {
            // potentials from previous time step no available.
            if (compressibility_)
            {
                densityWIJ = rhoMeanWIJ;
                densityNWIJ = rhoMeanNWIJ;
                densityWIK = rhoMeanWIK;
                densityNWIK = rhoMeanNWIK;
            }

            switch (pressureType_)
            {
            case pglobal:
            {
                Scalar pressJK = (cellDataJ.globalPressure() + cellDataK.globalPressure()) / 2;

                potentialWIJ = (cellDataI.globalPressure() - fMeanNWIJ * pcI - (pressJK - fMeanNWIJ * pcJK)) / l
                        + (densityWIJ - lJ / l * (kI + kK) / kI * (densityWIK - densityWIJ) / 2) * ng;
                potentialNWIJ = (cellDataI.globalPressure() + fMeanWIJ * pcI - (pressJK + fMeanWIJ * pcJK)) / l
                        + (densityNWIJ - lJ / l * (kI + kK) / kI * (densityNWIK - densityNWIJ) / 2) * ng;
                potentialWIK = (cellDataI.globalPressure() - fMeanNWIK * pcI - (pressJK - fMeanNWIK * pcJK)) / l
                        + (densityWIK - lJ / l * (kI + kK) / kI * (densityWIJ - densityWIK) / 2) * ng;
                potentialNWIK = (cellDataI.globalPressure() + fMeanWIK * pcI - (pressJK + fMeanWIK * pcJK)) / l
                        + (densityNWIK - lJ / l * (kI + kK) / kI * (densityNWIJ - densityNWIK) / 2) * ng;
                break;
            }
            default:
            {
                potentialWIJ = (cellDataI.pressure(wPhaseIdx)
                        - 0.5 * (cellDataJ.pressure(wPhaseIdx) + cellDataK.pressure(wPhaseIdx))) / l
                        + (densityWIJ - lJ / l * (kI + kK) / kI * (densityWIK - densityWIJ) / 2) * ng;
                potentialNWIJ = (cellDataI.pressure(nPhaseIdx)
                        - (0.5 * (cellDataJ.pressure(nPhaseIdx) + cellDataK.pressure(nPhaseIdx)))) / l
                        + (densityNWIJ - lJ / l * (kI + kK) / kI * (densityNWIK - densityNWIJ) / 2) * ng;
                potentialWIK = (cellDataI.pressure(wPhaseIdx)
                        - 0.5 * (cellDataJ.pressure(wPhaseIdx) + cellDataK.pressure(wPhaseIdx))) / l
                        + (densityWIK - lJ / l * (kI + kK) / kI * (densityWIJ - densityWIK) / 2) * ng;
                potentialNWIK = (cellDataI.pressure(nPhaseIdx)
                        - (0.5 * (cellDataJ.pressure(nPhaseIdx) + cellDataK.pressure(nPhaseIdx)))) / l
                        + (densityNWIK - lJ / l * (kI + kK) / kI * (densityNWIJ - densityNWIK) / 2) * ng;
                break;
            }
            }
        }

        //do the upwinding of the mobility depending on the phase potentials
        Scalar lambdaWIJ = (potentialWIJ > 0.) ? lambdaWI : lambdaWJ;
        lambdaWIJ = (potentialWIJ == 0) ? 0.5 * (lambdaWI + lambdaWJ) : lambdaWIJ;
        Scalar lambdaNWIJ = (potentialNWIJ > 0) ? lambdaNWI : lambdaNWJ;
        lambdaNWIJ = (potentialNWIJ == 0) ? 0.5 * (lambdaNWI + lambdaNWJ) : lambdaNWIJ;

        if (compressibility_)
        {
            densityWIJ = (potentialWIJ > 0.) ? cellDataI.density(wPhaseIdx) : cellDataJ.density(wPhaseIdx);
            densityNWIJ = (potentialNWIJ > 0.) ? cellDataI.density(nPhaseIdx) : cellDataJ.density(nPhaseIdx);
            densityWIJ = (potentialWIJ == 0) ? rhoMeanWIJ : densityWIJ;
            densityNWIJ = (potentialNWIJ == 0) ? rhoMeanNWIJ : densityNWIJ;
            densityWIK = (potentialWIK > 0.) ? cellDataI.density(wPhaseIdx) : cellDataK.density(wPhaseIdx);
            densityNWIK = (potentialNWIK > 0.) ? cellDataI.density(nPhaseIdx) : densityNWIK;
            densityWIK = (potentialWIK == 0) ? rhoMeanWIK : densityWIK;
            densityNWIK = (potentialNWIK == 0) ? rhoMeanNWIK : densityNWIK;
        }

        // update diagonal entry and right hand side
        entries[matrix] = (lambdaWIJ + lambdaNWIJ) * kMean / l * faceArea;
        entries[rhs] = faceArea * lambdaWIJ * kMean * ng
                * ((densityWIJ) - (lJ / l) * (kI + kK) / kI * (densityWIK - densityWIJ) / 2);
        entries[rhs] += faceArea * lambdaNWIJ * kMean * ng
                * ((densityNWIJ) - (lJ / l) * (kI + kK) / kI * (densityNWIK - densityNWIJ) / 2);

        switch (pressureType_)
        {
        case pw:
        {
            entries[rhs] += faceArea * lambdaNWIJ * kMean * (pcJK - pcI) / l;
            break;
        }
        case pn:
        {
            entries[rhs] -= faceArea * lambdaWIJ * kMean * (pcJK - pcI) / l;
            break;
        }
        }

        //write hanging-node-specific stuff directly into matrix and rhs!
        this->f_[globalIdxI] -= entries[rhs];
        this->f_[globalIdxJ] += entries[rhs];

        // set diagonal entry
        this->A_[globalIdxI][globalIdxI] += entries[matrix];
        //set off-diagonal
        this->A_[globalIdxI][globalIdxJ] -= entries[matrix];

        // set entries for cell J
        this->A_[globalIdxJ][globalIdxI] -= entries[matrix];
        this->A_[globalIdxJ][globalIdxJ] += entries[matrix];

        //set entries to zero -> matrix already written!
        entries = 0.;

//        std::cout<<"finished hanging node!\n";
    }
    else
    {
        entries = 0;
    }

    return;
}

}
#endif
