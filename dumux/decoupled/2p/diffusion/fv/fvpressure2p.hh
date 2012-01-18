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
#include <dumux/decoupled/common/fv/fvpressure.hh>
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
template<class TypeTag> class FVPressure2P: public FVPressure<TypeTag>
{
    typedef FVPressure<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;

    typedef typename SolutionTypes::ScalarSolution ScalarSolutionType;

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
        rhs = ParentType::rhs, matrix = ParentType::matrix,
    };

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;

public:
    void getSource(Dune::FieldVector<Scalar, 2>&, const Element&, const CellData&, const bool);

    void getStorage(Dune::FieldVector<Scalar, 2>&, const Element&, const CellData&, const bool);

    void getFlux(Dune::FieldVector<Scalar, 2>&, const Intersection&, const CellData&, const bool);

    void getFluxOnBoundary(Dune::FieldVector<Scalar, 2>&,
    const Intersection&, const CellData&, const bool);

    //! updates and stores constitutive relations
    void updateMaterialLaws();

    //! \copydoc Dumux::FVPressure1P::initialize()
    void initialize(bool solveTwice = true)
    {
        ParentType::initialize();
        updateMaterialLaws();

        this->assemble(true);
        this->solve();
        if (solveTwice)
        {
            Dune::BlockVector < Dune::FieldVector<Scalar, 1>
                    > pressureOld(this->pressure());

            this->assemble(false);
            this->solve();

            Dune::BlockVector < Dune::FieldVector<Scalar, 1> > pressureDiff(pressureOld);
            pressureDiff -= this->pressure();
            pressureOld = this->pressure();
            Scalar pressureNorm = pressureDiff.infinity_norm();
            pressureNorm /= pressureOld.infinity_norm();
            int numIter = 0;
            while (pressureNorm > 1e-5 && numIter < 10)
            {
                updateMaterialLaws();
                this->assemble(false);
                this->solve();

                pressureDiff = pressureOld;
                pressureDiff -= this->pressure();
                pressureNorm = pressureDiff.infinity_norm();
                pressureOld = this->pressure();
                pressureNorm /= pressureOld.infinity_norm();

                numIter++;
            }
            //            std::cout<<"Pressure defect = "<<pressureNorm<<"; "<<numIter<<" Iterations needed for initial pressure field"<<std::endl;
        }

        storePressureSolution();

        return;
    }

    //! updates the pressure field (analog to update function in Dumux::IMPET)
    void update()
    {
        //error bounds for error term for incompressible models to correct unphysical saturation over/undershoots due to saturation transport
        if (!compressibility_)
        {
            timeStep_ = problem_.timeManager().timeStepSize();
            maxError_ = 0.0;
            int size = problem_.gridView().size(0);
            for (int i = 0; i < size; i++)
            {
                Scalar sat = 0;
                switch (saturationType_)
                {
                case Sw:
                    sat = problem_.variables().cellData(i).saturation(wPhaseIdx);
                    break;
                case Sn:
                    sat = problem_.variables().cellData(i).saturation(nPhaseIdx);
                    break;
                }
                if (sat > 1.0)
                {
                    maxError_ = std::max(maxError_, (sat - 1.0) / timeStep_);
                }
                if (sat < 0.0)
                {
                    maxError_ = std::max(maxError_, (-sat) / timeStep_);
                }
            }
        }

        updateMaterialLaws();

        ParentType::update();

        storePressureSolution();

        return;
    }

    void storePressureSolution()
    {
        int size = problem_.gridView().size(0);
        for (int i = 0; i < size; i++)
        {
            CellData& cellData = problem_.variables().cellData(i);
            storePressureSolution(i, cellData);
        }
    }

    void storePressureSolution(int globalIdx, CellData& cellData)
    {
        switch (pressureType_)
        {
        case pw:
        {
            Scalar pressW = this->pressure()[globalIdx];
            Scalar pc = cellData.capillaryPressure();

                cellData.setPressure(wPhaseIdx, pressW);
                cellData.setPressure(nPhaseIdx, pressW + pc);

            break;
        }
        case pn:
        {
            Scalar pressNW = this->pressure()[globalIdx];
            Scalar pc = cellData.capillaryPressure();

                cellData.setPressure(nPhaseIdx, pressNW);
                cellData.setPressure(wPhaseIdx, pressNW - pc);

            break;
        }
        case pglobal:
        {
            cellData.setGlobalPressure(this->pressure()[globalIdx]);
            break;
        }
        }
    }

    //! \copydoc Dumux::FVPressure1P::addOutputVtkFields(MultiWriter &writer)
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        int size = problem_.gridView().size(0);
        ScalarSolutionType *pressureW = writer.allocateManagedBuffer(size);
        ScalarSolutionType *pressureN = writer.allocateManagedBuffer(size);
        ScalarSolutionType *pC = writer.allocateManagedBuffer(size);

        for (int i = 0; i < size; i++)
        {
            CellData& cellData = problem_.variables().cellData(i);
            storePressureSolution(i, cellData);
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
            ScalarSolutionType *densityWetting = writer.allocateManagedBuffer(size);
            ScalarSolutionType *densityNonwetting = writer.allocateManagedBuffer(size);
            ScalarSolutionType *viscosityWetting = writer.allocateManagedBuffer(size);
            ScalarSolutionType *viscosityNonwetting = writer.allocateManagedBuffer(size);

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
            ParentType(problem), problem_(problem), gravity_(problem.gravity()), maxError_(0.), timeStep_(1.)
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

        ErrorTermFactor_ = GET_PARAM(TypeTag, Scalar, ErrorTermFactor);
        ErrorTermLowerBound_ = GET_PARAM(TypeTag, Scalar, ErrorTermLowerBound);
        ErrorTermUpperBound_ = GET_PARAM(TypeTag, Scalar, ErrorTermUpperBound);

        if (!compressibility_)
        {
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
    }

private:
    Problem& problem_;
    const GlobalPosition& gravity_; //!< vector including the gravity constant

    Scalar maxError_;
    Scalar timeStep_;
    Scalar ErrorTermFactor_; //!< Handling of error term: relaxation factor
    Scalar ErrorTermLowerBound_; //!< Handling of error term: lower bound for error dampening
    Scalar ErrorTermUpperBound_; //!< Handling of error term: upper bound for error dampening

    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];

    static const bool compressibility_ = GET_PROP_VALUE(TypeTag, EnableCompressibility);
    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PressureFormulation); //!< gives kind of pressure used (\f$p_w\f$, \f$p_n\f$, \f$p_{global}\f$)
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, SaturationFormulation); //!< gives kind of saturation used (\f$S_w\f$, \f$S_n\f$)
};

//!function which calculates the source entry
template<class TypeTag>
void FVPressure2P<TypeTag>::getSource(Dune::FieldVector<Scalar, 2>& entries, const Element& element
        , const CellData& cellData, const bool first)
{
    // cell volume, assume linear map here
    Scalar volume = element.geometry().volume();

    // get sources from problem
    PrimaryVariables sourcePhase(0.0);
    problem_.source(sourcePhase, element);

    if (!compressibility_)
    {
        sourcePhase[wPhaseIdx] /= density_[wPhaseIdx];
        sourcePhase[nPhaseIdx] /= density_[nPhaseIdx];
    }

    entries[rhs] = volume * (sourcePhase[wPhaseIdx] + sourcePhase[nPhaseIdx]);

    return;
}

//!function which calculates the storage entry
template<class TypeTag>
void FVPressure2P<TypeTag>::getStorage(Dune::FieldVector<Scalar, 2>& entries, const Element& element
        , const CellData& cellData, const bool first)
{
    //volume correction due to density differences
    if (compressibility_ && !first)
    {
        // cell volume, assume linear map here
        Scalar volume = element.geometry().volume();

        Scalar porosity = problem_.spatialParameters().porosity(element);

        switch (saturationType_)
        {
        case Sw:
        {
            entries[rhs] = -(cellData.volumeCorrection() * porosity * volume
                    * (cellData.density(wPhaseIdx) - cellData.density(nPhaseIdx)));
            break;
        }
        case Sn:
        {
            entries[rhs] = -(cellData.volumeCorrection() * porosity * volume
                    * (cellData.density(nPhaseIdx) - cellData.density(wPhaseIdx)));
            break;
        }
        }
    }
    else if (!compressibility_ && !first)
    {
        //error term for incompressible models to correct unphysical saturation over/undershoots due to saturation transport
        // error reduction routine: volumetric error is damped and inserted to right hand side
        Scalar sat = 0;
        switch (saturationType_)
        {
        case Sw:
            sat = cellData.saturation(wPhaseIdx);
            break;
        case Sn:
            sat = cellData.saturation(nPhaseIdx);
            break;
        }

        Scalar volume = element.geometry().volume();

        Scalar error = (sat > 1.0) ? sat - 1.0 : 0.0;
        if (sat < 0.0) {error =  sat;}
        error /= timeStep_;

        Scalar errorAbs = std::abs(error);

        if ((errorAbs*timeStep_ > 1e-6) && (errorAbs > ErrorTermLowerBound_ * maxError_) && (!problem_.timeManager().willBeFinished()))
        {
            entries[rhs] = ErrorTermFactor_ * error * volume;
        }
    }

    return;
}

//!function which calculates internal flux entries
template<class TypeTag>
void FVPressure2P<TypeTag>::getFlux(Dune::FieldVector<Scalar, 2>& entries, const Intersection& intersection
        , const CellData& cellDataI, const bool first)
{
    ElementPointer elementI = intersection.inside();
    ElementPointer elementJ = intersection.outside();

    const CellData& cellDataJ = problem_.variables().cellData(problem_.variables().index(*elementJ));

    // get global coordinates of cell centers
    const GlobalPosition& globalPosI = elementI->geometry().center();
    const GlobalPosition& globalPosJ = elementJ->geometry().center();

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

    // distance vector between barycenters
    GlobalPosition distVec = globalPosJ - globalPosI;

    // compute distance between cell centers
    Scalar dist = distVec.two_norm();

    // compute vectorized permeabilities
    FieldMatrix meanPermeability(0);

    problem_.spatialParameters().meanK(meanPermeability, problem_.spatialParameters().intrinsicPermeability(*elementI),
            problem_.spatialParameters().intrinsicPermeability(*elementJ));

    Dune::FieldVector<Scalar, dim> permeability(0);
    meanPermeability.mv(unitOuterNormal, permeability);

    Scalar rhoMeanW = 0;
    Scalar rhoMeanNW = 0;
    if (compressibility_)
    {
        rhoMeanW = 0.5 * (cellDataI.density(wPhaseIdx) + cellDataJ.density(wPhaseIdx));
        rhoMeanNW = 0.5 * (cellDataI.density(nPhaseIdx) + cellDataJ.density(nPhaseIdx));
    }
    Scalar fMeanW = 0.5 * (fractionalWI + fractionalWJ);
    Scalar fMeanNW = 0.5 * (fractionalNWI + fractionalNWJ);

    //calculate potential gradients
    Scalar potentialW = 0;
    Scalar potentialNW = 0;

    //if we are at the very first iteration we can't calculate phase potentials
    if (!first)
    {
        potentialW = cellDataI.fluxData().potential(wPhaseIdx, isIndexI);
        potentialNW = cellDataI.fluxData().potential(nPhaseIdx, isIndexI);

        if (compressibility_)
        {
            density_[wPhaseIdx] = (potentialW > 0.) ? cellDataI.density(wPhaseIdx) : cellDataJ.density(wPhaseIdx);
            density_[nPhaseIdx] = (potentialNW > 0.) ? cellDataI.density(nPhaseIdx) : cellDataJ.density(nPhaseIdx);

            density_[wPhaseIdx] = (potentialW == 0.) ? rhoMeanW : density_[wPhaseIdx];
            density_[nPhaseIdx] = (potentialNW == 0.) ? rhoMeanNW : density_[nPhaseIdx];
        }

        potentialW = cellDataI.pressure(wPhaseIdx) - cellDataJ.pressure(wPhaseIdx);
        potentialNW = cellDataI.pressure(nPhaseIdx) - cellDataJ.pressure(nPhaseIdx);

        if (pressureType_ == pglobal)
        {
            potentialW = (cellDataI.globalPressure() - cellDataJ.globalPressure() - fMeanNW * (pcI - pcJ));
            potentialNW = (cellDataI.globalPressure() - cellDataJ.globalPressure() + fMeanW * (pcI - pcJ));
        }

        potentialW += density_[wPhaseIdx] * (distVec * gravity_);
        potentialNW += density_[nPhaseIdx] * (distVec * gravity_);
    }

    //do the upwinding of the mobility depending on the phase potentials
    Scalar lambdaW = (potentialW > 0.) ? lambdaWI : lambdaWJ;
    lambdaW = (potentialW == 0) ? 0.5 * (lambdaWI + lambdaWJ) : lambdaW;
    Scalar lambdaNW = (potentialNW > 0) ? lambdaNWI : lambdaNWJ;
    lambdaNW = (potentialNW == 0) ? 0.5 * (lambdaNWI + lambdaNWJ) : lambdaNW;

    if (compressibility_)
    {
        density_[wPhaseIdx] = (potentialW > 0.) ? cellDataI.density(wPhaseIdx) : cellDataJ.density(wPhaseIdx);
        density_[nPhaseIdx] = (potentialNW > 0.) ? cellDataI.density(nPhaseIdx) : cellDataJ.density(nPhaseIdx);

        density_[wPhaseIdx] = (potentialW == 0) ? rhoMeanW : density_[wPhaseIdx];
        density_[nPhaseIdx] = (potentialNW == 0) ? rhoMeanNW : density_[nPhaseIdx];
    }

    //calculate current matrix entry
    entries[matrix] = (lambdaW + lambdaNW) * ((permeability * unitOuterNormal) / dist) * faceArea;

    //calculate right hand side
    entries[rhs] = (lambdaW * density_[wPhaseIdx] + lambdaNW * density_[nPhaseIdx]) * (permeability * gravity_) * faceArea;

    switch (pressureType_)
    {
    case pw:
    {
        // calculate capillary pressure gradient
        Dune::FieldVector<Scalar, dim> pCGradient = unitOuterNormal;
        pCGradient *= (pcI - pcJ) / dist;

        //add capillary pressure term to right hand side
        entries[rhs] += 0.5 * (lambdaNWI + lambdaNWJ) * (permeability * pCGradient) * faceArea;
        break;
    }
    case pn:
    {
        // calculate capillary pressure gradient
        Dune::FieldVector<Scalar, dim> pCGradient = unitOuterNormal;
        pCGradient *= (pcI - pcJ) / dist;

        //add capillary pressure term to right hand side
        entries[rhs] -= 0.5 * (lambdaWI + lambdaWJ) * (permeability * pCGradient) * faceArea;
        break;
    }
    }

    return;
}

//!function which calculates internal flux entries
template<class TypeTag>
void FVPressure2P<TypeTag>::getFluxOnBoundary(Dune::FieldVector<Scalar, 2>& entries,
const Intersection& intersection, const CellData& cellData, const bool first)
{
    ElementPointer element = intersection.inside();

    // get global coordinates of cell centers
    const GlobalPosition& globalPosI = element->geometry().center();

    // center of face in global coordinates
    const GlobalPosition& globalPosJ = intersection.geometry().center();

    // get mobilities and fractional flow factors
    Scalar lambdaWI = cellData.mobility(wPhaseIdx);
    Scalar lambdaNWI = cellData.mobility(nPhaseIdx);
    Scalar fractionalWI = cellData.fracFlowFunc(wPhaseIdx);
    Scalar fractionalNWI = cellData.fracFlowFunc(nPhaseIdx);

    // get capillary pressure
    Scalar pcI = cellData.capillaryPressure();

    //get face index
    int isIndexI = intersection.indexInInside();

    //get face normal
    const Dune::FieldVector<Scalar, dim>& unitOuterNormal = intersection.centerUnitOuterNormal();

    // get face area
    Scalar faceArea = intersection.geometry().volume();

    // distance vector between barycenters
    GlobalPosition distVec = globalPosJ - globalPosI;

    // compute distance between cell centers
    Scalar dist = distVec.two_norm();

    BoundaryTypes bcType;
    problem_.boundaryTypes(bcType, intersection);
    PrimaryVariables boundValues(0.0);

    if (bcType.isDirichlet(eqIdxPress))
    {
        problem_.dirichlet(boundValues, intersection);

        //permeability vector at boundary
        // compute vectorized permeabilities
        FieldMatrix meanPermeability(0);

        problem_.spatialParameters().meanK(meanPermeability,
                problem_.spatialParameters().intrinsicPermeability(*element));

        Dune::FieldVector<Scalar, dim> permeability(0);
        meanPermeability.mv(unitOuterNormal, permeability);

        //determine saturation at the boundary -> if no saturation is known directly at the boundary use the cell saturation
        Scalar satW = 0;
        Scalar satNW = 0;
        if (bcType.isDirichlet(eqIdxSat))
        {
            switch (saturationType_)
            {
            case Sw:
            {
                satW = boundValues[saturationIdx];
                satNW = 1 - boundValues[saturationIdx];
                break;
            }
            case Sn:
            {
                satW = 1 - boundValues[saturationIdx];
                satNW = boundValues[saturationIdx];
                break;
            }
            }
        }
        else
        {
            satW = cellData.saturation(wPhaseIdx);
            satNW = cellData.saturation(nPhaseIdx);
        }
        Scalar temperature = problem_.temperature(*element);

        //get dirichlet pressure boundary condition
        Scalar pressBound = boundValues[pressureIdx];

        //calculate consitutive relations depending on the kind of saturation used
        Scalar pcBound = MaterialLaw::pC(problem_.spatialParameters().materialLawParams(*element), satW);

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

        Scalar densityWBound = density_[wPhaseIdx];
        Scalar densityNWBound = density_[nPhaseIdx];
        Scalar viscosityWBound = viscosity_[wPhaseIdx];
        Scalar viscosityNWBound = viscosity_[nPhaseIdx];
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

            rhoMeanW = 0.5 * (cellData.density(wPhaseIdx) + densityWBound);
            rhoMeanNW = 0.5 * (cellData.density(nPhaseIdx) + densityNWBound);
        }

        Scalar lambdaWBound = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(*element), satW)
                / viscosityWBound;
        Scalar lambdaNWBound = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(*element), satW)
                / viscosityNWBound;

        Scalar fractionalWBound = lambdaWBound / (lambdaWBound + lambdaNWBound);
        Scalar fractionalNWBound = lambdaNWBound / (lambdaWBound + lambdaNWBound);

        Scalar fMeanW = 0.5 * (fractionalWI + fractionalWBound);
        Scalar fMeanNW = 0.5 * (fractionalNWI + fractionalNWBound);

        Scalar potentialW = 0;
        Scalar potentialNW = 0;

        if (!first)
        {
            potentialW = cellData.fluxData().potential(wPhaseIdx, isIndexI);
            potentialNW = cellData.fluxData().potential(nPhaseIdx, isIndexI);

            if (compressibility_)
            {
                density_[wPhaseIdx] = (potentialW > 0.) ? cellData.density(wPhaseIdx) : densityWBound;
                density_[nPhaseIdx] = (potentialNW > 0.) ? cellData.density(nPhaseIdx) : densityNWBound;

                density_[wPhaseIdx] = (potentialW == 0.) ? rhoMeanW : density_[wPhaseIdx];
                density_[nPhaseIdx] = (potentialNW == 0.) ? rhoMeanNW : density_[nPhaseIdx];
            }

            //calculate potential gradient
            switch (pressureType_)
            {
            case pw:
            {
                potentialW = (cellData.pressure(wPhaseIdx) - pressBound);
                potentialNW = (cellData.pressure(nPhaseIdx) - pressBound - pcBound);
                break;
            }
            case pn:
            {
                potentialW = (cellData.pressure(wPhaseIdx) - pressBound + pcBound);
                potentialNW = (cellData.pressure(nPhaseIdx) - pressBound);
                break;
            }
            case pglobal:
            {
                potentialW = (cellData.globalPressure() - pressBound - fMeanNW * (pcI - pcBound));
                potentialNW = (cellData.globalPressure() - pressBound + fMeanW * (pcI - pcBound));
                break;
            }
            }

            potentialW += density_[wPhaseIdx] * (distVec * gravity_);
            potentialNW += density_[nPhaseIdx] * (distVec * gravity_);
        }

        //do the upwinding of the mobility depending on the phase potentials
        Scalar lambdaW = (potentialW > 0.) ? lambdaWI : lambdaWBound;
        lambdaW = (potentialW == 0) ? 0.5 * (lambdaWI + lambdaWBound) : lambdaW;
        Scalar lambdaNW = (potentialNW > 0.) ? lambdaNWI : lambdaNWBound;
        lambdaNW = (potentialNW == 0) ? 0.5 * (lambdaNWI + lambdaNWBound) : lambdaNW;

        if (compressibility_)
        {
            density_[wPhaseIdx] = (potentialW > 0.) ? cellData.density(wPhaseIdx) : densityWBound;
            density_[wPhaseIdx] = (potentialW == 0) ? rhoMeanW : density_[wPhaseIdx];
            density_[nPhaseIdx] = (potentialNW > 0.) ? cellData.density(nPhaseIdx) : densityNWBound;
            density_[nPhaseIdx] = (potentialNW == 0) ? rhoMeanNW : density_[nPhaseIdx];
        }

        //calculate current matrix entry
        entries[matrix] = (lambdaW + lambdaNW) * ((permeability * unitOuterNormal) / dist) * faceArea;
        entries[rhs] = entries[matrix] * pressBound;

        //calculate right hand side
        entries[rhs] -= (lambdaW * density_[wPhaseIdx] + lambdaNW * density_[nPhaseIdx]) * (permeability * gravity_)
                * faceArea;

        switch (pressureType_)
        {
        case pw:
        {
            // calculate capillary pressure gradient
            Dune::FieldVector<Scalar, dim> pCGradient = unitOuterNormal;
            pCGradient *= (pcI - pcBound) / dist;

            //add capillary pressure term to right hand side
            entries[rhs] -= 0.5 * (lambdaNWI + lambdaNWBound) * (permeability * pCGradient) * faceArea;
            break;
        }
        case pn:
        {
            // calculate capillary pressure gradient
            Dune::FieldVector<Scalar, dim> pCGradient = unitOuterNormal;
            pCGradient *= (pcI - pcBound) / dist;

            //add capillary pressure term to right hand side
            entries[rhs] += 0.5 * (lambdaWI + lambdaWBound) * (permeability * pCGradient) * faceArea;
            break;
        }
        }
    }
    //set neumann boundary condition
    else if (bcType.isNeumann(eqIdxPress))
    {
        problem_.neumann(boundValues, intersection);

        if (!compressibility_)
        {
            boundValues[wPhaseIdx] /= density_[wPhaseIdx];
            boundValues[nPhaseIdx] /= density_[nPhaseIdx];
        }
        entries[rhs] = -(boundValues[wPhaseIdx] + boundValues[nPhaseIdx]) * faceArea;
    }
    else
    {
        DUNE_THROW(Dune::NotImplemented, "No valid boundary condition type defined for pressure equation!");
    }

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


    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = problem_.gridView().template end<0>();
    for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eItEnd; ++eIt)
    {
        int globalIdx = problem_.variables().index(*eIt);

        CellData& cellData = problem_.variables().cellData(globalIdx);

        Scalar temperature = problem_.temperature(*eIt);

        //determine phase saturations from primary saturation variable

        Scalar satW = cellData.saturation(wPhaseIdx);
        Scalar satNW = cellData.saturation(nPhaseIdx);

        Scalar pc = MaterialLaw::pC(problem_.spatialParameters().materialLawParams(*eIt), satW);

        //determine phase pressures from primary pressure variable
        Scalar pressW = 0;
        Scalar pressNW = 0;
        switch (pressureType_)
        {
        case pw:
        {
            pressW = cellData.pressure(wPhaseIdx);
            pressNW = pressW + pc;
            break;
        }
        case pn:
        {
            pressNW = cellData.pressure(nPhaseIdx);
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

            density_[wPhaseIdx] = FluidSystem::density(fluidState, wPhaseIdx);
            density_[nPhaseIdx] = FluidSystem::density(fluidState, nPhaseIdx);

            viscosity_[wPhaseIdx]= FluidSystem::viscosity(fluidState, wPhaseIdx);
            viscosity_[nPhaseIdx] = FluidSystem::viscosity(fluidState, nPhaseIdx);

            //store density
            fluidState.setDensity(wPhaseIdx, density_[wPhaseIdx]);
            fluidState.setDensity(nPhaseIdx, density_[nPhaseIdx]);

            //store viscosity
            fluidState.setViscosity(wPhaseIdx, viscosity_[wPhaseIdx]);
            fluidState.setViscosity(nPhaseIdx, viscosity_[nPhaseIdx]);
        }
        else
        {
            cellData.setSaturation(wPhaseIdx, satW);
            cellData.setSaturation(nPhaseIdx, satNW);
            cellData.setCapillaryPressure(pc);

            if (pressureType_ != pglobal)
            {
                cellData.setPressure(wPhaseIdx, pressW);
                cellData.setPressure(nPhaseIdx, pressNW);
            }
        }

        // initialize mobilities
        Scalar mobilityW = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(*eIt), satW) / viscosity_[wPhaseIdx];
        Scalar mobilityNW = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(*eIt), satW) / viscosity_[nPhaseIdx];
//                std::cout<<"MobilityW: "<<mobilityW <<"\n"
//                        "MobilityNW"<< mobilityNW<<"\n";

        if (compressibility_)
        {
            mobilityW *= density_[wPhaseIdx];
            mobilityNW *= density_[nPhaseIdx];
        }

        // initialize mobilities
        cellData.setMobility(wPhaseIdx, mobilityW);
        cellData.setMobility(nPhaseIdx, mobilityNW);

        //initialize fractional flow functions
        cellData.setFracFlowFunc(wPhaseIdx, mobilityW / (mobilityW + mobilityNW));
        cellData.setFracFlowFunc(nPhaseIdx, mobilityNW / (mobilityW + mobilityNW));
    }
    return;
}

}
#endif
