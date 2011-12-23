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
#ifndef DUMUX_FVPRESSURE2P_HH
#define DUMUX_FVPRESSURE2P_HH

// dumux environment
#include <dumux/decoupled/2p/2pproperties.hh>

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
template<class TypeTag> class FVPressure2P
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Variables)) Variables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Indices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(WettingPhase)) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NonwettingPhase)) NonwettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(CellData)) CellData;

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
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases))
    };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureCoefficientMatrix)) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PressureRHSVector)) Vector;

    //initializes the matrix to store the system of equations
    void initializeMatrix();

    //function which assembles the system of equations to be solved
    void assemble(bool first);

    //solves the system of equations to get the spatial distribution of the pressure
    void solve();

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
    //! updates and stores constitutive relations
    void updateMaterialLaws();
    //! \copydoc Dumux::FVPressure1P::initialize()
    void initialize(bool solveTwice = true)
    {
        updateMaterialLaws();

        assemble(true);
        solve();
        if (solveTwice)
        {
            Dune::BlockVector < Dune::FieldVector<Scalar, 1>
                    > pressureOld(problem_.variables().primaryVariablesGlobal(eqIdxPress));

            assemble(false);
            solve();

            Dune::BlockVector < Dune::FieldVector<Scalar, 1> > pressureDiff(pressureOld);
            pressureDiff -= problem_.variables().primaryVariablesGlobal(eqIdxPress);
            pressureOld = problem_.variables().primaryVariablesGlobal(eqIdxPress);
            Scalar pressureNorm = pressureDiff.infinity_norm();
            pressureNorm /= pressureOld.infinity_norm();
            int numIter = 0;
            while (pressureNorm > 1e-5 && numIter < 10)
            {
                updateMaterialLaws();
                assemble(false);
                solve();

                pressureDiff = pressureOld;
                pressureDiff -= problem_.variables().primaryVariablesGlobal(eqIdxPress);
                pressureNorm = pressureDiff.infinity_norm();
                pressureOld = problem_.variables().primaryVariablesGlobal(eqIdxPress);
                pressureNorm /= pressureOld.infinity_norm();

                numIter++;
            }
            //            std::cout<<"Pressure defect = "<<pressureNorm<<"; "<<numIter<<" Iterations needed for initial pressure field"<<std::endl;
        }
        return;
    }
    //! \copydoc Dumux::FVPressure1P::pressure()
    void pressure(bool solveTwice = true)
    {
        assemble(false);
        solve();

        return;
    }
    //! updates the pressure field (analog to update function in Dumux::IMPET)
    void update()
    {
        updateMaterialLaws();

        pressure(false);

        return;
    }

    void storePressureSolution()
    {
        int size = problem_.gridView().size(0);
        for (int i = 0; i < size; i++)
        {
            CellData& cellData = problem_.variables().cellData(i);
            updatePressureSolution(i, cellData);
        }
    }

    void storePressureSolution(int globalIdx, CellData& cellData)
    {
        switch (pressureType_)
        {
        case pw:
        {
            Scalar pressW = problem_.variables().primaryVariablesGlobal(eqIdxPress)[globalIdx];
            Scalar pc = cellData.capillaryPressure();

            if (compressibility_)
            {
                cellData.fluidState().setPressure(wPhaseIdx, pressW);
                cellData.fluidState().setPressure(nPhaseIdx, pressW + pc);
            }
            else
            {
                cellData.setPressure(wPhaseIdx, pressW);
                cellData.setPressure(nPhaseIdx, pressW + pc);
            }
            break;
        }
        case pn:
        {
            Scalar pressNW = problem_.variables().primaryVariablesGlobal(eqIdxPress)[globalIdx];
            Scalar pc = cellData.capillaryPressure();
            if (compressibility_)
            {
            cellData.fluidState().setPressure(nPhaseIdx, pressNW);
            cellData.fluidState().setPressure(wPhaseIdx, pressNW - pc);
            }
            else
            {
                cellData.setPressure(nPhaseIdx, pressNW);
                cellData.setPressure(wPhaseIdx, pressNW - pc);
            }
            break;
        }
        case pglobal:
        {
            cellData.globalPressure() = problem_.variables().primaryVariablesGlobal(eqIdxPress)[globalIdx];
            break;
        }
        }
    }

    // serialization methods
    //! \copydoc Dumux::FVPressure1P::serialize(Restarter &res)
    template<class Restarter>
    void serialize(Restarter &res)
    {
        return;
    }

    //! \copydoc Dumux::FVPressure1P::deserialize(Restarter &res)
    template<class Restarter>
    void deserialize(Restarter &res)
    {
        return;
    }

    //! \copydoc Dumux::FVPressure1P::addOutputVtkFields(MultiWriter &writer)
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        int size = problem_.gridView().size(0);
        typename Variables::ScalarSolutionType *pressureW = writer.allocateManagedBuffer(size);
        typename Variables::ScalarSolutionType *pressureN = writer.allocateManagedBuffer(size);
        typename Variables::ScalarSolutionType *pC = writer.allocateManagedBuffer(size);

        for (int i = 0; i < size; i++)
        {
            CellData& cellData = problem_.variables().cellData(i);
            (*pressureW)[i] = cellData.pressure(wPhaseIdx);
            (*pressureN)[i] = cellData.pressure(nPhaseIdx);
            (*pC)[i] = cellData.capillaryPressure();
            if (pressureType_ == pglobal)
            {
                (*pressureW)[i] = cellData.globalPressure();
            }
        }
        if (pressureType_ == pglobal)
        {

            writer.attachCellData(*pressureW, "global pressure");
        }
        else
        {
            writer.attachCellData(*pressureW, "wetting pressure");
            writer.attachCellData(*pressureN, "nonwetting pressure");
        }

        writer.attachCellData(*pC, "capillary pressure");

        if (compressibility_)
        {
            typename Variables::ScalarSolutionType *densityWetting = writer.allocateManagedBuffer(size);
            typename Variables::ScalarSolutionType *densityNonwetting = writer.allocateManagedBuffer(size);
            typename Variables::ScalarSolutionType *viscosityWetting = writer.allocateManagedBuffer(size);
            typename Variables::ScalarSolutionType *viscosityNonwetting = writer.allocateManagedBuffer(size);

            for (int i = 0; i < size; i++)
            {
                CellData& cellData = problem_.variables().cellData(i);
                (*densityWetting)[i] = cellData.density(wPhaseIdx);
                (*densityNonwetting)[i] = cellData.density(nPhaseIdx);
                (*viscosityWetting)[i] = cellData.viscosity(wPhaseIdx);
                (*viscosityNonwetting)[i] = cellData.viscosity(nPhaseIdx);
            }

            writer.attachCellData(*densityWetting, "wetting density");
            writer.attachCellData(*densityNonwetting, "nonwetting density");
            writer.attachCellData(*viscosityWetting, "wetting viscosity");
            writer.attachCellData(*viscosityNonwetting, "nonwetting viscosity");
        }

        return;
    }

    //! Constructs a FVPressure2P object
    /**
     * \param problem a problem class object
     */
    FVPressure2P(Problem& problem) :
            problem_(problem), A_(problem.gridView().size(0), problem.gridView().size(0),
                    (2 * dim + 1) * problem.gridView().size(0), Matrix::random), f_(problem.gridView().size(0)), gravity(
                    problem.gravity())
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

        initializeMatrix();
    }

private:
    Problem& problem_;
    Matrix A_;
    Dune::BlockVector<Dune::FieldVector<Scalar, 1> > f_;

    static const bool compressibility_ = GET_PROP_VALUE(TypeTag, PTAG(EnableCompressibility));
    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PTAG(PressureFormulation)); //!< gives kind of pressure used (\f$p_w\f$, \f$p_n\f$, \f$p_{global}\f$)
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, PTAG(SaturationFormulation)); //!< gives kind of saturation used (\f$S_w\f$, \f$S_n\f$)
protected:
    const GlobalPosition& gravity; //!< vector including the gravity constant
};

//initializes the matrix to store the system of equations
template<class TypeTag>
void FVPressure2P<TypeTag>::initializeMatrix()
{
    // determine matrix row sizes
    ElementIterator eItEnd = problem_.gridView().template end<0>();
    for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        // initialize row size
        int rowSize = 1;

        // run through all intersections with neighbors
        IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
        for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt != isItEnd; ++isIt)
            if (isIt->neighbor())
                rowSize++;
        A_.setrowsize(globalIdxI, rowSize);
    }
    A_.endrowsizes();

    // determine position of matrix entries
    for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        // add diagonal index
        A_.addindex(globalIdxI, globalIdxI);

        // run through all intersections with neighbors
        IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
        for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt != isItEnd; ++isIt)
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer outside = isIt->outside();
                int globalIdxJ = problem_.variables().index(*outside);

                // add off diagonal index
                A_.addindex(globalIdxI, globalIdxJ);
            }
    }
    A_.endindices();

    return;
}

//!function which assembles the system of equations to be solved
template<class TypeTag>
void FVPressure2P<TypeTag>::assemble(bool first)
{
    // initialization: set matrix A_ to zero
    A_ = 0;
    f_ = 0;

    BoundaryTypes bcType;

    ElementIterator eItBegin = problem_.gridView().template begin<0>();
    Scalar densityW = 0;
    Scalar densityNW = 0;
    Scalar viscosityW = 0;
    Scalar viscosityNW = 0;

    if (!compressibility_)
    {
        Scalar temp = problem_.temperature(*eItBegin);
        Scalar pRef = problem_.referencePressure(*eItBegin);
        densityW = WettingPhase::density(temp, pRef);
        densityNW = NonwettingPhase::density(temp, pRef);
        viscosityW = WettingPhase::viscosity(temp, pRef);
        viscosityNW = NonwettingPhase::viscosity(temp, pRef);
    }

    ElementIterator eItEnd = problem_.gridView().template end<0>();
    for (ElementIterator eIt = eItBegin; eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        CellData& cellDataI = problem_.variables().cellData(globalIdxI);

        // get global coordinate of cell center
        const GlobalPosition& globalPos = eIt->geometry().center();

        // cell volume, assume linear map here
        Scalar volume = eIt->geometry().volume();

        // set sources
        PrimaryVariables source(0.0);
        problem().source(source, *eIt);
        if (!compressibility_)
        {
            source[wPhaseIdx] /= densityW;
            source[nPhaseIdx] /= densityNW;
        }
        f_[globalIdxI] = volume * (source[wPhaseIdx] + source[nPhaseIdx]);

        Scalar porosity = problem_.spatialParameters().porosity(*eIt);

        // get mobilities and fractional flow factors
        Scalar lambdaWI = cellDataI.mobility(wPhaseIdx);
        Scalar lambdaNWI = cellDataI.mobility(nPhaseIdx);
        Scalar fractionalWI = cellDataI.fracFlowFunc(wPhaseIdx);
        Scalar fractionalNWI = cellDataI.fracFlowFunc(nPhaseIdx);
        Scalar pcI = cellDataI.capillaryPressure();

        IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
        for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt != isItEnd; ++isIt)
        {
            int isIndex = isIt->indexInInside();

            // get normal vector
            const GlobalPosition& unitOuterNormal = isIt->centerUnitOuterNormal();

            // get face volume
            Scalar faceArea = isIt->geometry().volume();

            // handle interior face
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = problem_.variables().index(*neighborPointer);

                CellData& cellDataJ = problem_.variables().cellData(globalIdxJ);

                // neighbor cell center in global coordinates
                const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().center();

                // distance vector between barycenters
                GlobalPosition distVec = globalPosNeighbor - globalPos;

                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                // compute vectorized permeabilities
                FieldMatrix meanPermeability(0);

                problem_.spatialParameters().meanK(meanPermeability,
                        problem_.spatialParameters().intrinsicPermeability(*eIt),
                        problem_.spatialParameters().intrinsicPermeability(*neighborPointer));

                Dune::FieldVector<Scalar, dim> permeability(0);
                meanPermeability.mv(unitOuterNormal, permeability);

                // get mobilities and fractional flow factors
                Scalar lambdaWJ = cellDataJ.mobility(wPhaseIdx);
                Scalar lambdaNWJ = cellDataJ.mobility(nPhaseIdx);
                Scalar fractionalWJ = cellDataJ.fracFlowFunc(wPhaseIdx);
                Scalar fractionalNWJ = cellDataJ.fracFlowFunc(nPhaseIdx);

                Scalar pcJ = cellDataJ.capillaryPressure();

                Scalar rhoMeanW = 0;
                ;
                Scalar rhoMeanNW = 0;
                if (compressibility_)
                {
                    rhoMeanW = 0.5 * (cellDataI.density(wPhaseIdx) + cellDataJ.density(wPhaseIdx));
                    rhoMeanNW = 0.5 * (cellDataI.density(nPhaseIdx) + cellDataJ.density(nPhaseIdx));
                }
                Scalar fMeanW = 0.5 * (fractionalWI + fractionalWJ);
                Scalar fMeanNW = 0.5 * (fractionalNWI + fractionalNWJ);

                // update diagonal entry
                Scalar entry;

                //calculate potential gradients
                Scalar potentialW = 0;
                Scalar potentialNW = 0;

                //if we are at the very first iteration we can't calculate phase potentials
                if (!first)
                {
                    potentialW = cellDataI.fluxData().potential(wPhaseIdx, isIndex);
                    potentialNW = cellDataI.fluxData().potential(nPhaseIdx, isIndex);

                    if (compressibility_)
                    {
                        densityW =
                                (potentialW > 0.) ? cellDataI.density(wPhaseIdx) : cellDataJ.density(wPhaseIdx);
                        densityNW =
                                (potentialNW > 0.) ? cellDataI.density(nPhaseIdx) : cellDataJ.density(nPhaseIdx);

                        densityW = (potentialW == 0.) ? rhoMeanW : densityW;
                        densityNW = (potentialNW == 0.) ? rhoMeanNW : densityNW;
                    }

                    potentialW = cellDataI.pressure(wPhaseIdx) - cellDataJ.pressure(wPhaseIdx);
                    potentialNW = cellDataI.pressure(nPhaseIdx) - cellDataJ.pressure(nPhaseIdx);

                    if (pressureType_ == pglobal)
                    {
                        potentialW = (cellDataI.globalPressure() - cellDataJ.globalPressure()
                                - fMeanNW * (pcI - pcJ));
                        potentialNW = (cellDataI.globalPressure() - cellDataJ.globalPressure()
                                + fMeanW * (pcI - pcJ));
                    }

                    potentialW += densityW * (distVec * gravity);
                    potentialNW += densityNW * (distVec * gravity);

                    //store potentials for further calculations (velocity, saturation, ...)
                    cellDataI.fluxData().setPotential(wPhaseIdx, isIndex, potentialW);
                    cellDataI.fluxData().setPotential(nPhaseIdx, isIndex, potentialNW);
                }

                //do the upwinding of the mobility depending on the phase potentials
                Scalar lambdaW = (potentialW > 0.) ? lambdaWI : lambdaWJ;
                lambdaW = (potentialW == 0) ? 0.5 * (lambdaWI + lambdaWJ) : lambdaW;
                Scalar lambdaNW = (potentialNW > 0) ? lambdaNWI : lambdaNWJ;
                lambdaNW = (potentialNW == 0) ? 0.5 * (lambdaNWI + lambdaNWJ) : lambdaNW;

                if (compressibility_)
                {
                    densityW = (potentialW > 0.) ? cellDataI.density(wPhaseIdx) : cellDataJ.density(wPhaseIdx);
                    densityNW = (potentialNW > 0.) ? cellDataI.density(nPhaseIdx) : cellDataJ.density(nPhaseIdx);

                    densityW = (potentialW == 0) ? rhoMeanW : densityW;
                    densityNW = (potentialNW == 0) ? rhoMeanNW : densityNW;
                }

                //calculate current matrix entry
                entry = (lambdaW + lambdaNW) * ((permeability * unitOuterNormal) / dist) * faceArea;

                //calculate right hand side
                Scalar rightEntry = (lambdaW * densityW + lambdaNW * densityNW) * (permeability * gravity) * faceArea;

                switch (pressureType_)
                {
                case pw:
                {
                    // calculate capillary pressure gradient
                    Dune::FieldVector<Scalar, dim> pCGradient = unitOuterNormal;
                    pCGradient *= (pcI - pcJ) / dist;

                    //add capillary pressure term to right hand side
                    rightEntry += 0.5 * (lambdaNWI + lambdaNWJ) * (permeability * pCGradient) * faceArea;
                    break;
                }
                case pn:
                {
                    // calculate capillary pressure gradient
                    Dune::FieldVector<Scalar, dim> pCGradient = unitOuterNormal;
                    pCGradient *= (pcI - pcJ) / dist;

                    //add capillary pressure term to right hand side
                    rightEntry -= 0.5 * (lambdaWI + lambdaWJ) * (permeability * pCGradient) * faceArea;
                    break;
                }
                }

                //set right hand side
                f_[globalIdxI] -= rightEntry;

                // set diagonal entry
                A_[globalIdxI][globalIdxI] += entry;

                // set off-diagonal entry
                A_[globalIdxI][globalIdxJ] = -entry;

//                if (eIt->partitionType() == Dune::OverlapEntity)
//                {
//                    A_[globalIdxI][globalIdxI] *= 0.5;
//                    A_[globalIdxI][globalIdxJ] *= 0.5;
//
//                    f_[globalIdxI] *= 0.5;
//                }
            }

            // boundary face

            else
            {
                // center of face in global coordinates
                const GlobalPosition& globalPosFace = isIt->geometry().center();

                problem().boundaryTypes(bcType, *isIt);
                PrimaryVariables boundValues(0.0);

                GlobalPosition distVec(globalPosFace - globalPos);
                Scalar dist = distVec.two_norm();

                if (bcType.isDirichlet(eqIdxPress))
                {
                    problem().dirichlet(boundValues, *isIt);

                    //permeability vector at boundary
                    // compute vectorized permeabilities
                    FieldMatrix meanPermeability(0);

                    problem_.spatialParameters().meanK(meanPermeability,
                            problem_.spatialParameters().intrinsicPermeability(*eIt));

                    Dune::FieldVector<Scalar, dim> permeability(0);
                    meanPermeability.mv(unitOuterNormal, permeability);

                    //determine saturation at the boundary -> if no saturation is known directly at the boundary use the cell saturation
                    Scalar satBound;
                    if (bcType.isDirichlet(eqIdxSat))
                    {
                        satBound = boundValues[saturationIdx];
                    }
                    else
                    {
                        satBound = problem_.variables().primaryVariablesGlobal(eqIdxSat)[globalIdxI];
                    }
                    Scalar temperature = problem_.temperature(*eIt);

                    //get dirichlet pressure boundary condition
                    Scalar pressBound = boundValues[pressureIdx];

                    //calculate consitutive relations depending on the kind of saturation used
                    //determine phase saturations from primary saturation variable
                    Scalar satW = 0;
                    Scalar satNW = 0;
                    switch (saturationType_)
                    {
                    case Sw:
                    {
                        satW = satBound;
                        satNW = 1 - satBound;
                        break;
                    }
                    case Sn:
                    {
                        satW = 1 - satBound;
                        satNW = satBound;
                        break;
                    }
                    }

                    Scalar pcI = cellDataI.capillaryPressure();
                    Scalar pcBound = MaterialLaw::pC(problem_.spatialParameters().materialLawParams(*eIt), satW);

                    //determine phase pressures from primary pressure variable
                    Scalar pressW = 0;
                    Scalar pressNW = 0;
                    switch (pressureType_)
                    {
                    case pw:
                    {
                        pressW = pressBound;
                        pressNW = pressBound + pcBound;
                        break;
                    }
                    case pn:
                    {
                        pressW = pressBound - pcBound;
                        pressNW = pressBound;
                        break;
                    }
                    }

                    Scalar densityWBound = densityW;
                    Scalar densityNWBound = densityNW;
                    Scalar viscosityWBound = viscosityW;
                    Scalar viscosityNWBound = viscosityNW;
                    Scalar rhoMeanW = 0;
                    Scalar rhoMeanNW = 0;

                    if (compressibility_)
                    {
                        FluidState fluidState;
                        fluidState.setSaturation(wPhaseIdx, satW);
                        fluidState.setSaturation(nPhaseIdx, satNW);
                        fluidState.setTemperature(temperature);
                        fluidState.setPressure(wPhaseIdx, pressW);
                        fluidState.setPressure(nPhaseIdx, pressNW);

                        densityWBound = FluidSystem::density(fluidState, wPhaseIdx);
                        densityNWBound = FluidSystem::density(fluidState, nPhaseIdx);
                        viscosityWBound = FluidSystem::viscosity(fluidState, wPhaseIdx) / densityWBound;
                        viscosityNWBound = FluidSystem::viscosity(fluidState, nPhaseIdx) / densityNWBound;

                        rhoMeanW = 0.5 * (cellDataI.density(wPhaseIdx) + densityWBound);
                        rhoMeanNW = 0.5 * (cellDataI.density(nPhaseIdx) + densityNWBound);
                    }

                    Scalar lambdaWBound = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(*eIt), satW)
                            / viscosityWBound;
                    Scalar lambdaNWBound = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(*eIt), satW)
                            / viscosityNWBound;

                    Scalar fractionalWBound = lambdaWBound / (lambdaWBound + lambdaNWBound);
                    Scalar fractionalNWBound = lambdaNWBound / (lambdaWBound + lambdaNWBound);

                    Scalar fMeanW = 0.5 * (fractionalWI + fractionalWBound);
                    Scalar fMeanNW = 0.5 * (fractionalNWI + fractionalNWBound);

                    Scalar potentialW = 0;
                    Scalar potentialNW = 0;

                    if (!first)
                    {
                        potentialW = cellDataI.fluxData().potential(wPhaseIdx, isIndex);
                        potentialNW = cellDataI.fluxData().potential(nPhaseIdx, isIndex);

                        if (compressibility_)
                        {
                            densityW = (potentialW > 0.) ? cellDataI.density(wPhaseIdx) : densityWBound;
                            densityNW = (potentialNW > 0.) ? cellDataI.density(nPhaseIdx) : densityNWBound;

                            densityW = (potentialW == 0.) ? rhoMeanW : densityW;
                            densityNW = (potentialNW == 0.) ? rhoMeanNW : densityNW;
                        }

                        //calculate potential gradient
                        switch (pressureType_)
                        {
                        case pw:
                        {
                            potentialW = (cellDataI.pressure(wPhaseIdx) - pressBound);
                            potentialNW = (cellDataI.pressure(nPhaseIdx) - pressBound - pcBound);
                            break;
                        }
                        case pn:
                        {
                            potentialW = (cellDataI.pressure(wPhaseIdx) - pressBound + pcBound);
                            potentialNW = (cellDataI.pressure(nPhaseIdx) - pressBound);
                            break;
                        }
                        case pglobal:
                        {
                            potentialW = (cellDataI.globalPressure() - pressBound - fMeanNW * (pcI - pcBound));
                            potentialNW = (cellDataI.globalPressure() - pressBound + fMeanW * (pcI - pcBound));
                            break;
                        }
                        }

                        potentialW += densityW * (distVec * gravity);
                        potentialNW += densityNW * (distVec * gravity);

                        //store potential gradients for further calculations
                        cellDataI.fluxData().setPotential(wPhaseIdx, isIndex, potentialW);
                        cellDataI.fluxData().setPotential(nPhaseIdx, isIndex, potentialNW);
                    }

                    //do the upwinding of the mobility depending on the phase potentials
                    Scalar lambdaW = (potentialW > 0.) ? lambdaWI : lambdaWBound;
                    lambdaW = (potentialW == 0) ? 0.5 * (lambdaWI + lambdaWBound) : lambdaW;
                    Scalar lambdaNW = (potentialNW > 0.) ? lambdaNWI : lambdaNWBound;
                    lambdaNW = (potentialNW == 0) ? 0.5 * (lambdaNWI + lambdaNWBound) : lambdaNW;

                    if (compressibility_)
                    {
                        densityW = (potentialW > 0.) ? cellDataI.density(wPhaseIdx) : densityWBound;
                        densityW = (potentialW == 0) ? rhoMeanW : densityW;
                        densityNW = (potentialNW > 0.) ? cellDataI.density(nPhaseIdx) : densityNWBound;
                        densityNW = (potentialNW == 0) ? rhoMeanNW : densityNW;
                    }

                    //calculate current matrix entry
                    Scalar entry = (lambdaW + lambdaNW) * ((permeability * unitOuterNormal) / dist) * faceArea;

                    //calculate right hand side
                    Scalar rightEntry = (lambdaW * densityW + lambdaNW * densityNW) * (permeability * gravity)
                            * faceArea;

                    switch (pressureType_)
                    {
                    case pw:
                    {
                        // calculate capillary pressure gradient
                        Dune::FieldVector<Scalar, dim> pCGradient = unitOuterNormal;
                        pCGradient *= (pcI - pcBound) / dist;

                        //add capillary pressure term to right hand side
                        rightEntry += 0.5 * (lambdaNWI + lambdaNWBound) * (permeability * pCGradient) * faceArea;
                        break;
                    }
                    case pn:
                    {
                        // calculate capillary pressure gradient
                        Dune::FieldVector<Scalar, dim> pCGradient = unitOuterNormal;
                        pCGradient *= (pcI - pcBound) / dist;

                        //add capillary pressure term to right hand side
                        rightEntry -= 0.5 * (lambdaWI + lambdaWBound) * (permeability * pCGradient) * faceArea;
                        break;
                    }
                    }

                    // set diagonal entry and right hand side entry
                    A_[globalIdxI][globalIdxI] += entry;
                    f_[globalIdxI] += entry * pressBound;
                    f_[globalIdxI] -= rightEntry;
                }
                //set neumann boundary condition

                else if (bcType.isNeumann(eqIdxPress))
                {
                    problem().neumann(boundValues, *isIt);

                    if (!compressibility_)
                    {
                        boundValues[wPhaseIdx] /= densityW;
                        boundValues[nPhaseIdx] /= densityNW;
                    }
                    f_[globalIdxI] -= (boundValues[wPhaseIdx] + boundValues[nPhaseIdx]) * faceArea;

                    //Assumes that the phases flow in the same direction at the neumann boundary, which is the direction of the total flux!!!
                    //needed to determine the upwind direction in the saturation equation
                    cellDataI.fluxData().setPotential(wPhaseIdx, isIndex, boundValues[wPhaseIdx]);
                    cellDataI.fluxData().setPotential(nPhaseIdx, isIndex, boundValues[nPhaseIdx]);
                }
                else
                {
                    DUNE_THROW(Dune::NotImplemented, "No valid boundary condition type defined for pressure equation!");
                }
            }
        } // end all intersections

        //volume correction due to density differences
        if (compressibility_)
        {
            switch (saturationType_)
            {
            case Sw:
            {
                f_[globalIdxI] -= cellDataI.volumeCorrection() * porosity * volume * (cellDataI.density(wPhaseIdx) - cellDataI.density(nPhaseIdx));
                break;
            }
            case Sn:
            {
                f_[globalIdxI] -= cellDataI.volumeCorrection() * porosity * volume * (cellDataI.density(nPhaseIdx) - cellDataI.density(wPhaseIdx));
                break;
            }
            }
        }
    } // end grid traversal
    return;
}

//!solves the system of equations to get the spatial distribution of the pressure
template<class TypeTag>
void FVPressure2P<TypeTag>::solve()
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LinearSolver)) Solver;

    int verboseLevelSolver = GET_PARAM(TypeTag, int, LinearSolver, Verbosity);

    if (verboseLevelSolver)
        std::cout << "FVPressure2P: solve for pressure" << std::endl;

//    printmatrix(std::cout, A_, "global stiffness matrix", "row", 11, 3);
//    printvector(std::cout, f_, "right hand side", "row", 200, 1, 3);

    Solver solver(problem_);
    solver.solve(A_, problem_.variables().primaryVariablesGlobal(eqIdxPress), f_);

//                    printvector(std::cout, (problem_.variables().pressure()), "pressure", "row", 200, 1, 3);

    //for parallel use
//    problem_.variables().communicatePressure();

    return;
}
//!constitutive functions are updated once if new saturations are calculated and stored in the variables object
template<class TypeTag>
void FVPressure2P<TypeTag>::updateMaterialLaws()
{
    //for parallel use
//    printvector(std::cout, (problem_.variables().pressure()), "pressure", "row", 200, 1, 3);
//    problem_.variables().communicatePressure();
//    printvector(std::cout, (problem_.variables().pressure()), "pressureComm", "row", 200, 1, 3);

//    printvector(std::cout, (problem_.variables().saturation()), "sat", "row", 200, 1, 3);
//    problem_.variables().communicateTransportedQuantity();
//    printvector(std::cout, (problem_.variables().saturation()), "satComm", "row", 200, 1, 3);

    ElementIterator eItBegin = problem_.gridView().template begin<0>();
    Scalar densityW = 0;
    Scalar densityNW = 0;
    Scalar viscosityW = 0;
    Scalar viscosityNW = 0;

    if (!compressibility_)
    {
        Scalar temp = problem_.temperature(*eItBegin);
        Scalar pRef = problem_.referencePressure(*eItBegin);
        densityW = WettingPhase::density(temp, pRef);
        densityNW = NonwettingPhase::density(temp, pRef);
        viscosityW = WettingPhase::viscosity(temp, pRef);
        viscosityNW = NonwettingPhase::viscosity(temp, pRef);
    }

    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = problem_.gridView().template end<0>();
    for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eItEnd; ++eIt)
    {
        int globalIdx = problem_.variables().index(*eIt);

        CellData& cellData = problem_.variables().cellData(globalIdx);

        Scalar temperature = problem_.temperature(*eIt);

        //determine phase saturations from primary saturation variable
        Scalar satW = 0;
        Scalar satNW = 0;
        switch (saturationType_)
        {
        case Sw:
        {
            satW = problem_.variables().primaryVariablesGlobal(eqIdxSat)[globalIdx];
            satNW = 1 - satW;
            break;
        }
        case Sn:
        {
            satNW = problem_.variables().primaryVariablesGlobal(eqIdxSat)[globalIdx];
            satW = 1 - satNW;
            break;
        }
        }

        Scalar pc = MaterialLaw::pC(problem_.spatialParameters().materialLawParams(*eIt), satW);

        //determine phase pressures from primary pressure variable
        Scalar pressW = 0;
        Scalar pressNW = 0;
        switch (pressureType_)
        {
        case pw:
        {
            pressW = problem_.variables().primaryVariablesGlobal(eqIdxPress)[globalIdx];
            pressNW = pressW + pc;
            break;
        }
        case pn:
        {
            pressNW = problem_.variables().primaryVariablesGlobal(eqIdxPress)[globalIdx];
            pressW = pressNW - pc;
            break;
        }
        }

        if (compressibility_)
        {
            FluidState& fluidState = cellData.fluidState();
            fluidState.setSaturation(wPhaseIdx, satW);
            fluidState.setSaturation(nPhaseIdx, satNW);
            fluidState.setTemperature(temperature);

            fluidState.setPressure(wPhaseIdx, pressW);
            fluidState.setPressure(nPhaseIdx, pressNW);

            densityW = FluidSystem::density(fluidState, wPhaseIdx);
            densityNW = FluidSystem::density(fluidState, nPhaseIdx);

            viscosityW = FluidSystem::viscosity(fluidState, wPhaseIdx);
            viscosityNW = FluidSystem::viscosity(fluidState, nPhaseIdx);

            //store density
            fluidState.setDensity(wPhaseIdx, densityW);
            fluidState.setDensity(nPhaseIdx, densityNW);

            //store viscosity
            fluidState.setViscosity(wPhaseIdx, viscosityW);
            fluidState.setViscosity(nPhaseIdx, viscosityNW);
        }
        else
        {
            cellData.setSaturation(wPhaseIdx, satW);
            cellData.setSaturation(nPhaseIdx, satNW);
            cellData.setCapillaryPressure(pc);

            if (pressureType_ == pglobal)
            {
                cellData.globalPressure() = problem_.variables().primaryVariablesGlobal(eqIdxPress)[globalIdx];
            }
            else
            {
                cellData.setPressure(wPhaseIdx, pressW);
                cellData.setPressure(nPhaseIdx, pressNW);
            }
        }

        // initialize mobilities
        Scalar mobilityW = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(*eIt), satW) / viscosityW;
        Scalar mobilityNW = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(*eIt), satW) / viscosityNW;
        //        std::cout<<"MobilityW: "<<mobilityW <<"\n"
        //                "MobilityNW"<< mobilityNW<<"\n";

        if (compressibility_)
        {
            mobilityW *= densityW;
            mobilityNW *= densityNW;
        }

        // initialize mobilities
        cellData.mobility(wPhaseIdx) = mobilityW;
        cellData.mobility(nPhaseIdx) = mobilityNW;

        //initialize fractional flow functions
        cellData.fracFlowFunc(wPhaseIdx) = mobilityW / (mobilityW + mobilityNW);
        cellData.fracFlowFunc(nPhaseIdx) = mobilityNW / (mobilityW + mobilityNW);
    }
    return;
}

}
#endif
