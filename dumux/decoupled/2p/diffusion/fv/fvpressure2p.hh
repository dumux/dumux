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
#ifndef DUMUX_FVPRESSURE2P_HH
#define DUMUX_FVPRESSURE2P_HH

// dumux environment
#include <dumux/decoupled/common/fv/fvpressure.hh>
#include <dumux/decoupled/2p/diffusion/diffusionproperties2p.hh>

/**
 * @file
 * @brief  Finite Volume Discretization of a two-phase flow pressure equation.
 * @author Markus Wolff, Bernd Flemisch, Jochen Fritz,
 */

namespace Dumux
{
//! \ingroup FVPressure2p
/*!  \brief Finite Volume discretization of a two-phase flow pressure equation of the sequential IMPES Model.
 *
 * This class provides a finite volume (FV) implementation for solving equations of the form
 * \f[
 * \phi \left( \rho_w  \frac{\partial S_w}{\partial t} + \rho_n \frac{\partial S_n}{\partial t}\right) + \text{div}\, \boldsymbol{v}_{total} = q.
 * \f]
 * The definition of the total velocity \f$\boldsymbol{v}_{total}\f$ depends on the choice of the primary pressure variable.
 * Further, fluids can be assumed to be compressible or incompressible (Property: <tt>EnableCompressibility</tt>).
 * In the incompressible case a wetting (\f$ w \f$) phase pressure as primary variable leads to
 *
 * \f[
 * - \text{div}\,  \left[\lambda \boldsymbol K \left(\text{grad}\, p_w + f_n \text{grad}\, p_c + \sum f_\alpha \rho_\alpha g  \text{grad}\, z\right)\right] = q,
 * \f]
 *
 * a non-wetting (\f$ n \f$) phase pressure yields
 * \f[
 *  - \text{div}\,  \left[\lambda \boldsymbol K  \left(\text{grad}\, p_n - f_w \text{grad}\, p_c + \sum f_\alpha \rho_\alpha g  \text{grad}\, z\right)\right] = q,
 *  \f]
 * and a global pressure leads to
 * \f[
 * - \text{div}\, \left[\lambda \boldsymbol K \left(\text{grad}\, p_{global} + \sum f_\alpha \rho_\alpha g  \text{grad}\, z\right)\right] = q.
 * \f]
 * Here, \f$ p_\alpha \f$ is a phase pressure, \f$ p_ {global} \f$ the global pressure of a classical fractional flow formulation
 * (see e.g. ﻿P. Binning and M. A. Celia, “Practical implementation of the fractional flow approach to multi-phase flow simulation,” Advances in water resources, vol. 22, no. 5, pp. 461-478, 1999.),
 * \f$ p_c = p_n - p_w \f$ is the capillary pressure, \f$ \boldsymbol K \f$ the absolute permeability, \f$ \lambda = \lambda_w +  \lambda_n \f$ the total mobility depending on the
 * saturation (\f$ \lambda_\alpha = k_{r_\alpha} / \mu_\alpha \f$),\f$ f_\alpha = \lambda_\alpha / \lambda \f$ the fractional flow function of a phase,
 * \f$ \rho_\alpha \f$ a phase density, \f$ g \f$ the gravity constant and \f$ q \f$ the source term.
 *
 * For all cases, \f$ p = p_D \f$ on \f$ \Gamma_{Neumann} \f$, and \f$ \boldsymbol{v}_{total}  = q_N \f$
 * on \f$ \Gamma_{Dirichlet} \f$.
 *
 * The slightly compressible case is only implemented for phase pressures! In this case for a wetting (\f$ w \f$) phase pressure as primary variable the equations are formulated as
 * \f[
 * \phi \left( \rho_w  \frac{\partial S_w}{\partial t} + \rho_n \frac{\partial S_n}{\partial t}\right) - \text{div}\,  \left[\lambda \boldsymbol{K} \left(\text{grad}\, p_w + f_n \text{grad}\, p_c + \sum f_\alpha \rho_\alpha g  \text{grad}\, z\right)\right] = q,
 * \f]
 * and for a non-wetting (\f$ n \f$) phase pressure as
 *  \f[
 *  \phi \left( \rho_w  \frac{\partial S_w}{\partial t} + \rho_n \frac{\partial S_n}{\partial t}\right) - \text{div}\,  \left[\lambda \boldsymbol{K}  \left(\text{grad}\, p_n - f_w \text{grad}\, p_c + \sum f_\alpha \rho_\alpha g  \text{grad}\, z\right)\right] = q,
 *  \f]
 * In this slightly compressible case the following definitions are valid:  \f$ \lambda = \rho_w \lambda_w + \rho_n \lambda_n \f$, \f$ f_\alpha = (\rho_\alpha \lambda_\alpha) / \lambda \f$
 * This model assumes that temporal changes in density are very small and thus terms of temporal derivatives are negligible in the pressure equation.
 * Depending on the formulation the terms including time derivatives of saturations are simplified by inserting  \f$ S_w + S_n = 1 \f$.
 *
 *  In the IMPES models the default setting is:
 *
 *      - formulation: \f$ p_w-S_w \f$ (Property: <tt>Formulation</tt> defined as <tt>DecoupledTwoPCommonIndices::pwSw</tt>)
 *      - compressibility: disabled (Property: <tt>EnableCompressibility</tt> set to <tt>false</tt>)
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

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;

protected:
    //! \cond \private
    enum
    {
        rhs = ParentType::rhs, matrix = ParentType::matrix
    };
    //! \endcond

public:

    //! \cond \private

    // Function which calculates the source entry
    void getSource(Dune::FieldVector<Scalar, 2>& entry, const Element& element, const CellData& cellData, const bool first);

    // Function which calculates the storage entry
    void getStorage(Dune::FieldVector<Scalar, 2>& entry, const Element& element, const CellData& cellData, const bool first);

    // Function which calculates the flux entry
    void getFlux(Dune::FieldVector<Scalar, 2>& entry, const Intersection& intersection, const CellData& cellData, const bool first);

    // Function which calculates the boundary flux entry
    void getFluxOnBoundary(Dune::FieldVector<Scalar, 2>& entry,
    const Intersection& intersection, const CellData& cellData, const bool first);

    // updates and stores constitutive relations
    void updateMaterialLaws();
    //! \endcond

    /*! \brief Initializes the pressure model
     *
     * \copydoc ParentType::initialize()
     *
     * \param solveTwice indicates if more than one iteration is allowed to get an initial pressure solution
     */
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

    /*! \brief Pressure update
     *
     * \copydoc ParentType::update()
     *
     */
    void update()
    {
        timeStep_ = problem_.timeManager().timeStepSize();
        //error bounds for error term for incompressible models to correct unphysical saturation over/undershoots due to saturation transport
        if (!compressibility_)
        {
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

        ParentType::update();

        storePressureSolution();

        return;
    }

    /*! \brief Globally stores the pressure solution
     *
     */
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
     * Calculates secondary pressure variables and stores pressures.
     *
     * \param globalIdx Global cell index
     * \param cellData A CellData object
     */
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

    /*! \brief Adds pressure output to the output file
     *
     * Adds the phase pressures or a global pressure (depending on the formulation) as well as the capillary pressure to the output.
     * In the compressible case additionally density and viscosity are added.
     *
     * \tparam MultiWriter Class defining the output writer
     * \param writer The output writer (usually a <tt>VTKMultiWriter</tt> object)
     *
     */
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

/*! \brief Function which calculates the source entry
 *
 * \copydoc FVPressure::getSource(Dune::FieldVector<Scalar,2>&,const Element&,const CellData&,const bool)
 *
 * Source of each fluid phase has to be added as mass flux (\f$\text{kg}/(\text{m}^3 \text{s}\f$).
 */
template<class TypeTag>
void FVPressure2P<TypeTag>::getSource(Dune::FieldVector<Scalar, 2>& entry, const Element& element
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

    entry[rhs] = volume * (sourcePhase[wPhaseIdx] + sourcePhase[nPhaseIdx]);

    return;
}

/** \brief Function which calculates the storage entry
 *
 * \copydoc FVPressure::getStorage(Dune::FieldVector<Scalar,2>&,const Element&,const CellData&,const bool)
 *
 * If compressibility is enabled this functions calculates the term
 * \f[
 *      \phi \sum_\alpha \rho_\alpha \frac{\partial S_\alpha}{\partial t} V
 * \f]
 *
 * In the incompressible case an volume correction term is calculated which corrects for unphysical saturation overshoots/undershoots.
 * These can occur if the estimated time step for the explicit transport was too large. Correction by an artificial source term allows to correct
 * this errors due to wrong time-stepping without losing mass conservation. The error term looks as follows:
 * \f[
 *  q_{error} = \begin{cases}
 *          S < 0 & a_{error} \frac{S}{\Delta t} V \\
 *          S > 1 & a_{error} \frac{(S - 1)}{\Delta t} V \\
 *          0 \le S \le 1 & 0
 *      \end{cases}
 *  \f]
 *  where \f$a_{error}\f$ is a weighting factor (default: \f$a_{error} = 0.5\f$)
 */
template<class TypeTag>
void FVPressure2P<TypeTag>::getStorage(Dune::FieldVector<Scalar, 2>& entry, const Element& element
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
            entry[rhs] = -(cellData.volumeCorrection()/timeStep_ * porosity * volume
                    * (cellData.density(wPhaseIdx) - cellData.density(nPhaseIdx)));
            break;
        }
        case Sn:
        {
            entry[rhs] = -(cellData.volumeCorrection()/timeStep_ * porosity * volume
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
            entry[rhs] = ErrorTermFactor_ * error * volume;
        }
    }

    return;
}

/*! \brief Function which calculates the flux entry
 *
 * \copydoc FVPressure::getFlux(Dune::FieldVector<Scalar,2>&,const Intersection&,const CellData&,const bool)
 *
 */
template<class TypeTag>
void FVPressure2P<TypeTag>::getFlux(Dune::FieldVector<Scalar, 2>& entry, const Intersection& intersection
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
    entry[matrix] = (lambdaW + lambdaNW) * ((permeability * unitOuterNormal) / dist) * faceArea;

    //calculate right hand side
    entry[rhs] = (lambdaW * density_[wPhaseIdx] + lambdaNW * density_[nPhaseIdx]) * (permeability * gravity_) * faceArea;

    switch (pressureType_)
    {
    case pw:
    {
        // calculate capillary pressure gradient
        Dune::FieldVector<Scalar, dim> pCGradient = unitOuterNormal;
        pCGradient *= (pcI - pcJ) / dist;

        //add capillary pressure term to right hand side
        entry[rhs] += 0.5 * (lambdaNWI + lambdaNWJ) * (permeability * pCGradient) * faceArea;
        break;
    }
    case pn:
    {
        // calculate capillary pressure gradient
        Dune::FieldVector<Scalar, dim> pCGradient = unitOuterNormal;
        pCGradient *= (pcI - pcJ) / dist;

        //add capillary pressure term to right hand side
        entry[rhs] -= 0.5 * (lambdaWI + lambdaWJ) * (permeability * pCGradient) * faceArea;
        break;
    }
    }

    return;
}

/*! \brief Function which calculates the flux entry at a boundary
 *
 * \copydoc FVPressure::getFluxOnBoundary(Dune::FieldVector<Scalar,2>&,const Intersection&,const CellData&,const bool)
 *
 * Dirichlet boundary condition is a pressure depending on the formulation (\f$p_w\f$ (default), \f$p_n\f$, \f$p_{global}\f$),
 * Neumann boundary condition are the phase mass fluxes (\f$q_w\f$ and \f$q_n\f$, [\f$\text{kg}/(\text{m}^2 \text{s}\f$])
 */
template<class TypeTag>
void FVPressure2P<TypeTag>::getFluxOnBoundary(Dune::FieldVector<Scalar, 2>& entry,
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
        entry[matrix] = (lambdaW + lambdaNW) * ((permeability * unitOuterNormal) / dist) * faceArea;
        entry[rhs] = entry[matrix] * pressBound;

        //calculate right hand side
        entry[rhs] -= (lambdaW * density_[wPhaseIdx] + lambdaNW * density_[nPhaseIdx]) * (permeability * gravity_)
                * faceArea;

        switch (pressureType_)
        {
        case pw:
        {
            // calculate capillary pressure gradient
            Dune::FieldVector<Scalar, dim> pCGradient = unitOuterNormal;
            pCGradient *= (pcI - pcBound) / dist;

            //add capillary pressure term to right hand side
            entry[rhs] -= 0.5 * (lambdaNWI + lambdaNWBound) * (permeability * pCGradient) * faceArea;
            break;
        }
        case pn:
        {
            // calculate capillary pressure gradient
            Dune::FieldVector<Scalar, dim> pCGradient = unitOuterNormal;
            pCGradient *= (pcI - pcBound) / dist;

            //add capillary pressure term to right hand side
            entry[rhs] += 0.5 * (lambdaWI + lambdaWBound) * (permeability * pCGradient) * faceArea;
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
        entry[rhs] = -(boundValues[wPhaseIdx] + boundValues[nPhaseIdx]) * faceArea;
    }
    else
    {
        DUNE_THROW(Dune::NotImplemented, "No valid boundary condition type defined for pressure equation!");
    }

    return;
}

/*! \brief Updates constitutive relations and stores them in the variable class
 *
 * Stores mobility, fractional flow function and capillary pressure for all grid cells. In the compressible case additionally the densities and viscosities are stored.
 */
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
