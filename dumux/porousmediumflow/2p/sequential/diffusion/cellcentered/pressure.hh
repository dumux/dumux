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
 * \brief  Finite Volume discretization of a two-phase flow pressure equation.
 */
#ifndef DUMUX_FVPRESSURE2P_HH
#define DUMUX_FVPRESSURE2P_HH

#include <dune/common/float_cmp.hh>

// dumux environment
#include <dumux/porousmediumflow/sequential/cellcentered/pressure.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/properties.hh>

#include <dumux/common/deprecated.hh>

namespace Dumux {
/*!
 * \ingroup SequentialTwoPModel
 * \brief Finite Volume discretization of a two-phase flow pressure equation of the sequential IMPES model.
 *
 * This model implements two-phase flow of two immiscible fluids \f$\alpha \in \{ w, n \}\f$ using
 * a standard multi-phase Darcy approach as the equation for the conservation of momentum, i.e.
 \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \textbf{K}
 \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} {\textbf g} \right).
 \f]
 * By inserting this into the equation for the conservation of the phase mass, one gets
  \f[
 \phi \frac{\partial \varrho_\alpha S_\alpha}{\partial t}
 -
 \text{div} \left\{
 \varrho_\alpha \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K} \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} \mathbf{g} \right)
 \right\} = q_\alpha \;.
 \f]

 * In the incompressible case the phase densities \f$ \varrho_\alpha \f$ can be eliminated from the equations.
 * The two equations are then added up, and because \f$S_w + S_n = 1\f$, the first term cancels out.
 * This leads to the so-called pressure equation. For the wetting (\f$w\f$)
 * phase pressure as primary variable this yields
 \f[
 - \text{div}\,  \left[\lambda \boldsymbol K \left(\textbf{grad}\, p_w + f_n \textbf{grad}\,
   p_c - \sum f_\alpha \varrho_\alpha {\textbf g}\right)\right] = q,
 \f]
 *
 * for the nonwetting (\f$ n \f$) phase pressure as primary variable it yields
 \f[
  - \text{div}\,  \left[\lambda \boldsymbol K  \left(\textbf{grad}\, p_n - f_w \textbf{grad}\,
    p_c - \sum f_\alpha \varrho_\alpha {\textbf g}\right)\right] = q,
  \f]
 * and for the global pressure as primary variable it leads to
 \f[
 - \text{div}\, \left[\lambda \boldsymbol K \left(\textbf{grad}\,
   p_{global} - \sum f_\alpha \varrho_\alpha {\textbf g}\right)\right] = q.
 \f]
 * Here, \f$ p_\alpha \f$ is a phase pressure, \f$ p_ {global} \f$ the global pressure of a classical fractional
 * flow formulation (see e.g. Binning and Celia (1999) \cite Binning1999),
 * \f$ p_c = p_n - p_w \f$ is the capillary pressure, \f$ \boldsymbol K \f$ the absolute permeability tensor,
 * \f$ \lambda = \lambda_w +  \lambda_n \f$ the total mobility depending on the
 * saturation (\f$ \lambda_\alpha = k_{r_\alpha} / \mu_\alpha \f$),
 * \f$ f_\alpha = \lambda_\alpha / \lambda \f$ the fractional flow function of a phase,
 * \f$ \varrho_\alpha \f$ a phase density, \f$ {\textbf g} \f$ the gravitational acceleration vector and
 * \f$ q = q_w + q_n \f$ the total source term.
 * Depending on the primary variable chosen, one of the pressure equations above is solved sequentially
 * together with the conservation of the phase mass for one phase.
 *
 * For all cases, \f$ p = p_D \f$ on \f$ \Gamma_{Dirichlet} \f$, and \f$ \boldsymbol v_{total} \cdot  \boldsymbol n  = q_N \f$
 * on \f$ \Gamma_{Neumann} \f$.
 *
 * The slightly compressible case is only implemented for phase pressures! In this case for a wetting
 * \f$(w) \f$ phase pressure as primary variable the pressure equation is formulated as
 * \f[
 * \phi \left( \varrho_w  \frac{\partial S_w}{\partial t} + \varrho_n \frac{\partial S_n}{\partial t}\right) - \text{div}\,
 * \left[\lambda \boldsymbol{K} \left(\textbf{grad}\, p_w + f_n \, \textbf{grad}\, p_c - \sum f_\alpha \varrho_\alpha
 * {\textbf g}\right)\right] = q,
 * \f]
 * and for a nonwetting (\f$ n \f$) phase pressure as primary variable as
 *  \f[
 *  \phi \left( \varrho_w  \frac{\partial S_w}{\partial t} + \varrho_n \frac{\partial S_n}{\partial t}\right) - \text{div}\,
 * \left[\lambda \boldsymbol{K}  \left(\textbf{grad}\, p_n - f_w \textbf{grad}\, p_c - \sum f_\alpha \varrho_\alpha
 * {\textbf g}\right)\right] = q.
 *  \f]
 * In this slightly compressible case the following definitions are valid:
 * \f$ \lambda = \varrho_w \lambda_w + \varrho_n \lambda_n \f$, \f$ f_\alpha = (\varrho_\alpha \lambda_\alpha) / \lambda \f$
 * This model assumes that temporal changes in density are very small and thus terms of temporal derivatives are negligible in the pressure equation.
 * Depending on the formulation the terms including time derivatives of saturations are simplified by inserting  \f$ S_w + S_n = 1 \f$.
 *
 *  In the IMPES models the default setting is:
 *
 *  - formulation: \f$ p_w-S_w \f$ (Property: \a Formulation defined as \a SequentialTwoPCommonIndices::pwsw)
 *
 *  - compressibility: disabled (Property: \a EnableCompressibility set to \a false)
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag> class FVPressure2P: public FVPressure<TypeTag>
{
    using ParentType = FVPressure<TypeTag>;

    //the model implementation
    using Implementation = GetPropType<TypeTag, Properties::PressureModel>;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using PrimaryVariables = typename SolutionTypes::PrimaryVariables;
    using CellData = GetPropType<TypeTag, Properties::CellData>;
    using PressureSolutionVector = GetPropType<TypeTag, Properties::PressureSolutionVector>;

    using ScalarSolutionType = typename SolutionTypes::ScalarSolution;

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
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        eqIdxPress = Indices::pressureEqIdx,
        eqIdxSat = Indices::satEqIdx
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx, numPhases = getPropValue<TypeTag, Properties::NumPhases>()
    };

    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GravityVector = Dune::FieldVector<Scalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

protected:
    //! \cond \private
    using EntryType = typename ParentType::EntryType;
    enum
    {
        rhs = ParentType::rhs, matrix = ParentType::matrix
    };
    //! \endcond

public:
    //! Function which calculates the source entry
    void getSource(EntryType& entry, const Element& element, const CellData& cellData, const bool first);

    //! Function which calculates the storage entry
    void getStorage(EntryType& entry, const Element& element, const CellData& cellData, const bool first);

    //! Function which calculates the flux entry
    void getFlux(EntryType& entry, const Intersection& intersection, const CellData& cellData, const bool first);

    //! Function which calculates the boundary flux entry
    void getFluxOnBoundary(EntryType& entry,
    const Intersection& intersection, const CellData& cellData, const bool first);

    //! Updates and stores constitutive relations
    void updateMaterialLaws();

    /*!
     * \brief Initializes the pressure model
     *
     * \copydetails ParentType::initialize()
     *
     * \param solveTwice indicates if more than one iteration is allowed to get an initial pressure solution
     */
    void initialize(bool solveTwice = true)
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

        storePressureSolution();
        updateMaterialLaws();

        this->assemble(true);
        this->solve();
        if (solveTwice)
        {
            PressureSolutionVector pressureOld(this->pressure());

            this->assemble(false);
            this->solve();

            PressureSolutionVector pressureDiff(pressureOld);
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
            // std::cout<<"Pressure defect = "<<pressureNorm<<"; "<<
            //        numIter<<" Iterations needed for initial pressure field"<<std::endl;
        }

        storePressureSolution();
    }

    /*!
     * \brief Pressure update
     *
     * \copydetails FVPressure::update()
     */
    void update()
    {
        timeStep_ = problem_.timeManager().timeStepSize();
        //error bounds for error term for incompressible models
        //to correct unphysical saturation over/undershoots due to saturation transport
        if (!compressibility_)
        {
            maxError_ = 0.0;
            int size = problem_.gridView().size(0);
            for (int i = 0; i < size; i++)
            {
                Scalar sat = 0;
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
                    using std::max;
                    maxError_ = max(maxError_, (sat - 1.0) / timeStep_);
                }
                if (sat < 0.0)
                {
                    using std::max;
                    maxError_ = max(maxError_, (-sat) / timeStep_);
                }
            }
        }

        ParentType::update();

        storePressureSolution();
    }

    /*!
     * \brief Velocity update
     *
     * Resets the velocities in the cellData.
     */
    void updateVelocity()
    {
        updateMaterialLaws();

        //reset velocities
        int size = problem_.gridView().size(0);
        for (int i = 0; i < size; i++)
        {
            CellData& cellData = problem_.variables().cellData(i);
            cellData.fluxData().resetVelocity();
        }
    }

    //! \brief Globally stores the pressure solution
    void storePressureSolution()
    {
        // iterate through leaf grid
        for (const auto& element : elements(problem_.gridView()))
        {
            asImp_().storePressureSolution(element);
        }
    }

    /*!
     * \brief Stores the pressure solution of a cell
     *
     * Calculates secondary pressure variables and stores pressures.
     *
     * \param element Grid element
     */
    void storePressureSolution(const Element& element)
    {
        int eIdxGlobal = problem_.variables().index(element);
        CellData& cellData = problem_.variables().cellData(eIdxGlobal);

        if (compressibility_)
        {
            density_[wPhaseIdx] = cellData.density(wPhaseIdx);
            density_[nPhaseIdx] = cellData.density(nPhaseIdx);
        }

        switch (pressureType_)
        {
        case pw:
        {
            Scalar pressW = this->pressure()[eIdxGlobal];
            Scalar pc = cellData.capillaryPressure();

            cellData.setPressure(wPhaseIdx, pressW);
            cellData.setPressure(nPhaseIdx, pressW + pc);

            Scalar gravityDiff = (problem_.bBoxMax() - element.geometry().center()) * gravity_;
            Scalar potW = pressW + gravityDiff * density_[wPhaseIdx];
            Scalar potNw = pressW + pc + gravityDiff * density_[nPhaseIdx];

            cellData.setPotential(wPhaseIdx, potW);
            cellData.setPotential(nPhaseIdx, potNw);

            break;
        }
        case pn:
        {
            Scalar pressNw = this->pressure()[eIdxGlobal];
            Scalar pc = cellData.capillaryPressure();

            cellData.setPressure(nPhaseIdx, pressNw);
            cellData.setPressure(wPhaseIdx, pressNw - pc);

            Scalar gravityDiff = (problem_.bBoxMax() - element.geometry().center()) * gravity_;
            Scalar potW = pressNw - pc + gravityDiff * density_[wPhaseIdx];
            Scalar potNw = pressNw + gravityDiff * density_[nPhaseIdx];

            cellData.setPotential(wPhaseIdx, potW);
            cellData.setPotential(nPhaseIdx, potNw);

            break;
        }
        case pGlobal:
        {
            Scalar press = this->pressure()[eIdxGlobal];
            cellData.setGlobalPressure(press);

            Scalar pc = cellData.capillaryPressure();
            Scalar gravityDiff = (problem_.bBoxMax() - element.geometry().center()) * gravity_;

            //This is only an estimation!!! -> only used for uwind directions and time-step estimation!!!
            Scalar potW = press - cellData.fracFlowFunc(nPhaseIdx) * pc + gravityDiff * density_[wPhaseIdx];
            Scalar potNw = press - cellData.fracFlowFunc(wPhaseIdx) * pc + gravityDiff * density_[nPhaseIdx];

            cellData.setPotential(wPhaseIdx, potW);
            cellData.setPotential(nPhaseIdx, potNw);

            break;
        }
        }
        cellData.fluxData().resetVelocity();
    }

    /*!
     * \brief Adds pressure output to the output file
     *
     * Adds the phase pressures or a global pressure (depending on the formulation) as well as the capillary pressure to the output.
     * In the compressible case additionally density and viscosity are added.
     *
     * \tparam MultiWriter Class defining the output writer
     * \param writer The output writer (usually a <tt>VTKMultiWriter</tt> object)
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        int size = problem_.gridView().size(0);
        ScalarSolutionType *pressure = writer.allocateManagedBuffer(size);
        ScalarSolutionType *pressureSecond = 0;
        ScalarSolutionType *pc = 0;
        ScalarSolutionType *potentialW = 0;
        ScalarSolutionType *potentialNw = 0;
        ScalarSolutionType *mobilityW = 0;
        ScalarSolutionType *mobilityNW = 0;

        if (vtkOutputLevel_ > 0)
        {
            pressureSecond = writer.allocateManagedBuffer(size);
            pc = writer.allocateManagedBuffer(size);
            mobilityW = writer.allocateManagedBuffer(size);
            mobilityNW = writer.allocateManagedBuffer(size);
        }
        if (vtkOutputLevel_ > 1)
        {
            potentialW = writer.allocateManagedBuffer(size);
            potentialNw = writer.allocateManagedBuffer(size);
        }

        for (int i = 0; i < size; i++)
        {
            CellData& cellData = problem_.variables().cellData(i);

            if (pressureType_ == pw)
            {
                (*pressure)[i] = cellData.pressure(wPhaseIdx);
                if (vtkOutputLevel_ > 0)
                {
                    (*pressureSecond)[i] = cellData.pressure(nPhaseIdx);
                }
            }
            else if (pressureType_ == pn)
            {
                (*pressure)[i] = cellData.pressure(nPhaseIdx);
                if (vtkOutputLevel_ > 0)
                {
                    (*pressureSecond)[i] = cellData.pressure(wPhaseIdx);
                }
            }
            else if (pressureType_ == pGlobal)
            {
                (*pressure)[i] = cellData.globalPressure();
            }
            if (vtkOutputLevel_ > 0)
            {
                (*pc)[i] = cellData.capillaryPressure();
                (*mobilityW)[i]  = cellData.mobility(wPhaseIdx);
                (*mobilityNW)[i]  = cellData.mobility(nPhaseIdx);
                if (compressibility_)
                {
                    (*mobilityW)[i] = (*mobilityW)[i]/cellData.density(wPhaseIdx);
                    (*mobilityNW)[i] = (*mobilityNW)[i]/cellData.density(nPhaseIdx);
                }
            }
            if (vtkOutputLevel_ > 1)
            {
                (*potentialW)[i]  = cellData.potential(wPhaseIdx);
                (*potentialNw)[i]  = cellData.potential(nPhaseIdx);
            }
        }

        if (pressureType_ == pw)
        {
            writer.attachCellData(*pressure, "wetting pressure");
            if (vtkOutputLevel_ > 0)
            {
                writer.attachCellData(*pressureSecond, "nonwetting pressure");
            }
        }
        else if (pressureType_ == pn)
        {
            writer.attachCellData(*pressure, "nonwetting pressure");
            if (vtkOutputLevel_ > 0)
            {
                writer.attachCellData(*pressureSecond, "wetting pressure");
            }
        }
        if (pressureType_ == pGlobal)
        {

            writer.attachCellData(*pressure, "global pressure");
        }

        if (vtkOutputLevel_ > 0)
        {
            writer.attachCellData(*pc, "capillary pressure");
            writer.attachCellData(*mobilityW, "wetting mobility");
            writer.attachCellData(*mobilityNW, "nonwetting mobility");
        }
        if (vtkOutputLevel_ > 1)
        {
            writer.attachCellData(*potentialW, "wetting potential");
            writer.attachCellData(*potentialNw, "nonwetting potential");
        }

        if (compressibility_)
        {
            if (vtkOutputLevel_ > 0)
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
        }
    }

    /*!
     * \brief Constructs a FVPressure2P object
     *
     * \param problem A problem class object
     */
    FVPressure2P(Problem& problem) :
            ParentType(problem), problem_(problem), gravity_(problem.gravity()), maxError_(0.), timeStep_(1.)
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
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc IMPETProblem::asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    Problem& problem_;
    const GravityVector& gravity_; //!< vector including the gravity constant

    Scalar maxError_;
    Scalar timeStep_;
    Scalar ErrorTermFactor_; //!< Handling of error term: relaxation factor
    Scalar ErrorTermLowerBound_; //!< Handling of error term: lower bound for error dampening
    Scalar ErrorTermUpperBound_; //!< Handling of error term: upper bound for error dampening

    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];

    int vtkOutputLevel_;

    static const bool compressibility_ = getPropValue<TypeTag, Properties::EnableCompressibility>();
    //! Gives kind of pressure used (\f$p_w\f$, \f$p_n\f$, \f$p_{global}\f$)
    static const int pressureType_ = getPropValue<TypeTag, Properties::PressureFormulation>();
    //! Gives kind of saturation used (\f$S_w\f$, \f$S_n\f$)
    static const int saturationType_ = getPropValue<TypeTag, Properties::SaturationFormulation>();
};

/*!
 * \brief Function which calculates the source entry
 *
 * \copydetails FVPressure::getSource(EntryType&,const Element&,const CellData&,const bool)
 *
 * Source of each fluid phase has to be added as mass flux (\f$\text{kg}/(\text{m}^3 \text{s}\f$).
 */
template<class TypeTag>
void FVPressure2P<TypeTag>::getSource(EntryType& entry, const Element& element
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
}

/*!
 * \brief Function which calculates the storage entry
 *
 * \copydetails FVPressure::getStorage(EntryType&,const Element&,const CellData&,const bool)
 *
 * If compressibility is enabled this functions calculates the term
 * \f[
 *      \phi \sum_\alpha \varrho_\alpha \frac{\partial S_\alpha}{\partial t} V
 * \f]
 *
 * In the incompressible case an volume correction term is calculated which corrects
 * for non-physical saturation overshoots/undershoots.
 * These can occur if the estimated time step for the explicit transport was too large.
 * Correction by an artificial source term allows to correct this errors due to wrong time-stepping
 * without losing mass conservation. The error term looks as follows:
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
void FVPressure2P<TypeTag>::getStorage(EntryType& entry, const Element& element
        , const CellData& cellData, const bool first)
{
    //volume correction due to density differences
    if (compressibility_ && !first)
    {
        // cell volume, assume linear map here
        Scalar volume = element.geometry().volume();

        Scalar porosity = problem_.spatialParams().porosity(element);

        switch (saturationType_)
        {
        case sw:
        {
            entry[rhs] = -(cellData.volumeCorrection()/timeStep_ * porosity * volume
                    * (cellData.density(wPhaseIdx) - cellData.density(nPhaseIdx)));
            break;
        }
        case sn:
        {
            entry[rhs] = -(cellData.volumeCorrection()/timeStep_ * porosity * volume
                    * (cellData.density(nPhaseIdx) - cellData.density(wPhaseIdx)));
            break;
        }
        }
    }
    else if (!compressibility_ && !first)
    {
        // error term for incompressible models to correct non-physical saturation over/undershoots due to saturation transport
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

        Scalar volume = element.geometry().volume();

        Scalar error = (sat > 1.0) ? sat - 1.0 : 0.0;
        if (sat < 0.0) {error =  sat;}
        error /= timeStep_;

        using std::abs;
        Scalar errorAbs = abs(error);

        if ((errorAbs*timeStep_ > 1e-6) && (errorAbs > ErrorTermLowerBound_ * maxError_) && (!problem_.timeManager().willBeFinished()))
        {
            entry[rhs] = ErrorTermFactor_ * error * volume;
        }
    }
}

/*!
 * \brief Function which calculates the flux entry
 *
 * \copydetails FVPressure::getFlux(EntryType&,const Intersection&,const CellData&,const bool)
 */
template<class TypeTag>
void FVPressure2P<TypeTag>::getFlux(EntryType& entry, const Intersection& intersection
        , const CellData& cellData, const bool first)
{
    auto elementI = intersection.inside();
    auto elementJ = intersection.outside();

    const CellData& cellDataJ = problem_.variables().cellData(problem_.variables().index(elementJ));

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

    // get face normal
    const Dune::FieldVector<Scalar, dim>& unitOuterNormal = intersection.centerUnitOuterNormal();

    // get face area
    Scalar faceArea = intersection.geometry().volume();

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

    Scalar rhoMeanW = 0;
    Scalar rhoMeanNw = 0;
    if (compressibility_)
    {
        rhoMeanW = 0.5 * (cellData.density(wPhaseIdx) + cellDataJ.density(wPhaseIdx));
        rhoMeanNw = 0.5 * (cellData.density(nPhaseIdx) + cellDataJ.density(nPhaseIdx));
    }

    // calculate potential gradients
    Scalar potentialDiffW = 0;
    Scalar potentialDiffNw = 0;

    // if we are at the very first iteration we can't calculate phase potentials
    if (!first)
    {
        potentialDiffW = cellData.potential(wPhaseIdx) - cellDataJ.potential(wPhaseIdx);
        potentialDiffNw = cellData.potential(nPhaseIdx) - cellDataJ.potential(nPhaseIdx);

        if (compressibility_)
        {
            density_[wPhaseIdx] = (potentialDiffW > 0.) ? cellData.density(wPhaseIdx) : cellDataJ.density(wPhaseIdx);
            density_[nPhaseIdx] = (potentialDiffNw > 0.) ? cellData.density(nPhaseIdx) : cellDataJ.density(nPhaseIdx);

            density_[wPhaseIdx] = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffW, 0.0, 1.0e-30)) ? rhoMeanW : density_[wPhaseIdx];
            density_[nPhaseIdx] = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNw, 0.0, 1.0e-30)) ? rhoMeanNw : density_[nPhaseIdx];

            potentialDiffW = cellData.pressure(wPhaseIdx) - cellDataJ.pressure(wPhaseIdx);
            potentialDiffNw = cellData.pressure(nPhaseIdx) - cellDataJ.pressure(nPhaseIdx);

            potentialDiffW += density_[wPhaseIdx] * (distVec * gravity_);
            potentialDiffNw += density_[nPhaseIdx] * (distVec * gravity_);
        }
    }

    // do the upwinding of the mobility depending on the phase potentials
    Scalar lambdaW = (potentialDiffW > 0.) ? lambdaWI : lambdaWJ;
    lambdaW = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffW, 0.0, 1.0e-30)) ? 0.5 * (lambdaWI + lambdaWJ) : lambdaW;
    Scalar lambdaNw = (potentialDiffNw > 0) ? lambdaNwI : lambdaNwJ;
    lambdaNw = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNw, 0.0, 1.0e-30)) ? 0.5 * (lambdaNwI + lambdaNwJ) : lambdaNw;

    if (compressibility_)
    {
        density_[wPhaseIdx] = (potentialDiffW > 0.) ? cellData.density(wPhaseIdx) : cellDataJ.density(wPhaseIdx);
        density_[nPhaseIdx] = (potentialDiffNw > 0.) ? cellData.density(nPhaseIdx) : cellDataJ.density(nPhaseIdx);

        density_[wPhaseIdx] = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffW, 0.0, 1.0e-30)) ? rhoMeanW : density_[wPhaseIdx];
        density_[nPhaseIdx] = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNw, 0.0, 1.0e-30)) ? rhoMeanNw : density_[nPhaseIdx];
    }

    Scalar scalarPerm = permeability.two_norm();
    // calculate current matrix entry
    entry[matrix] = (lambdaW + lambdaNw) * scalarPerm / dist * faceArea;

    // calculate right hand side
    // calculate unit distVec
    distVec /= dist;
    Scalar areaScaling = (unitOuterNormal * distVec);
    // this treatment of g allows to account for gravity flux through faces where the face normal has no z component (e.g. parallelepiped grids)
    entry[rhs] = (lambdaW * density_[wPhaseIdx] + lambdaNw * density_[nPhaseIdx]) * scalarPerm * (gravity_ * distVec) * faceArea * areaScaling;

    if (pressureType_ == pw)
    {
        // add capillary pressure term to right hand side
        entry[rhs] += 0.5 * (lambdaNwI + lambdaNwJ) * scalarPerm * (pcI - pcJ) / dist * faceArea;
    }
    else if (pressureType_ == pn)
    {
        // add capillary pressure term to right hand side
        entry[rhs] -= 0.5 * (lambdaWI + lambdaWJ) * scalarPerm * (pcI - pcJ) / dist * faceArea;
    }
}

/*!
 * \brief Function which calculates the flux entry at a boundary
 *
 * \copydetails FVPressure::getFluxOnBoundary(EntryType&,const Intersection&,const CellData&,const bool)
 *
 * Dirichlet boundary condition for pressure depends on the formulation (\f$p_w\f$ (default), \f$p_n\f$, \f$p_{global}\f$),
 * Neumann boundary condition are the phase mass fluxes (\f$q_w\f$ and \f$q_n\f$, [\f$\text{kg}/(\text{m}^2 \text{s}\f$])
 */
template<class TypeTag>
void FVPressure2P<TypeTag>::getFluxOnBoundary(EntryType& entry,
const Intersection& intersection, const CellData& cellData, const bool first)
{
    auto element = intersection.inside();

    // get global coordinates of cell centers
    const GlobalPosition& globalPosI = element.geometry().center();

    // center of face in global coordinates
    const GlobalPosition& globalPosJ = intersection.geometry().center();

    // get mobilities and fractional flow factors
    Scalar lambdaWI = cellData.mobility(wPhaseIdx);
    Scalar lambdaNwI = cellData.mobility(nPhaseIdx);
    Scalar fractionalWI = cellData.fracFlowFunc(wPhaseIdx);
    Scalar fractionalNwI = cellData.fracFlowFunc(nPhaseIdx);

    // get capillary pressure
    Scalar pcI = cellData.capillaryPressure();

    // get face index
    int isIndexI = intersection.indexInInside();

    // get face normal
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
        DimMatrix meanPermeability(0);

        problem_.spatialParams().meanK(meanPermeability,
                problem_.spatialParams().intrinsicPermeability(element));

        Dune::FieldVector<Scalar, dim> permeability(0);
        meanPermeability.mv(unitOuterNormal, permeability);

        // determine saturation at the boundary -> if no saturation is known directly at the boundary use the cell saturation
        Scalar satW = 0;
        Scalar satNw = 0;
        if (bcType.isDirichlet(eqIdxSat))
        {
            switch (saturationType_)
            {
            case sw:
            {
                satW = boundValues[saturationIdx];
                satNw = 1 - boundValues[saturationIdx];
                break;
            }
            case sn:
            {
                satW = 1 - boundValues[saturationIdx];
                satNw = boundValues[saturationIdx];
                break;
            }
            }
        }
        else
        {
            satW = cellData.saturation(wPhaseIdx);
            satNw = cellData.saturation(nPhaseIdx);
        }

        const Scalar temperature = problem_.temperature(element);

        // get Dirichlet pressure boundary condition
        const Scalar pressBound = boundValues[pressureIdx];

        //calculate constitutive relations depending on the kind of saturation used

        // old material law interface is deprecated: Replace this by
        // const auto& fluidMatrixInteraction = spatialParams.fluidMatrixInteractionAtPos(element.geometry().center());
        // after the release of 3.3, when the deprecated interface is no longer supported
        const auto fluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem_.spatialParams(), element);

        const Scalar pcBound = fluidMatrixInteraction.pc(satW);

        // determine phase pressures from primary pressure variable
        Scalar pressW = 0;
        Scalar pressNw = 0;
        if (pressureType_ == pw)
        {
            pressW = pressBound;
            pressNw = pressBound + pcBound;
        }
        else if (pressureType_ ==  pn)
        {
            pressW = pressBound - pcBound;
            pressNw = pressBound;
        }

        Scalar densityWBound = density_[wPhaseIdx];
        Scalar densityNwBound = density_[nPhaseIdx];
        Scalar viscosityWBound = viscosity_[wPhaseIdx];
        Scalar viscosityNwBound = viscosity_[nPhaseIdx];
        Scalar rhoMeanW = 0;
        Scalar rhoMeanNw = 0;

        if (compressibility_)
        {
            FluidState fluidState;
            fluidState.setSaturation(wPhaseIdx, satW);
            fluidState.setSaturation(nPhaseIdx, satNw);
            fluidState.setTemperature(temperature);
            fluidState.setPressure(wPhaseIdx, pressW);
            fluidState.setPressure(nPhaseIdx, pressNw);

            densityWBound = FluidSystem::density(fluidState, wPhaseIdx);
            densityNwBound = FluidSystem::density(fluidState, nPhaseIdx);
            viscosityWBound = FluidSystem::viscosity(fluidState, wPhaseIdx) / densityWBound;
            viscosityNwBound = FluidSystem::viscosity(fluidState, nPhaseIdx) / densityNwBound;

            rhoMeanW = 0.5 * (cellData.density(wPhaseIdx) + densityWBound);
            rhoMeanNw = 0.5 * (cellData.density(nPhaseIdx) + densityNwBound);
        }

        Scalar lambdaWBound = fluidMatrixInteraction.krw(satW)
                / viscosityWBound;
        Scalar lambdaNwBound = fluidMatrixInteraction.krn(satW)
                / viscosityNwBound;

        Scalar fractionalWBound = lambdaWBound / (lambdaWBound + lambdaNwBound);
        Scalar fractionalNwBound = lambdaNwBound / (lambdaWBound + lambdaNwBound);

        Scalar fMeanW = 0.5 * (fractionalWI + fractionalWBound);
        Scalar fMeanNw = 0.5 * (fractionalNwI + fractionalNwBound);

        Scalar potentialDiffW = 0;
        Scalar potentialDiffNw = 0;

        if (!first)
        {
            potentialDiffW = cellData.fluxData().upwindPotential(wPhaseIdx, isIndexI);
            potentialDiffNw = cellData.fluxData().upwindPotential(nPhaseIdx, isIndexI);

            if (compressibility_)
            {
                density_[wPhaseIdx] = (potentialDiffW > 0.) ? cellData.density(wPhaseIdx) : densityWBound;
                density_[nPhaseIdx] = (potentialDiffNw > 0.) ? cellData.density(nPhaseIdx) : densityNwBound;

                density_[wPhaseIdx] = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffW, 0.0, 1.0e-30)) ? rhoMeanW : density_[wPhaseIdx];
                density_[nPhaseIdx] = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNw, 0.0, 1.0e-30)) ? rhoMeanNw : density_[nPhaseIdx];
            }

            // calculate potential gradient
            switch (pressureType_)
            {
            case pw:
            {
                potentialDiffW = (cellData.pressure(wPhaseIdx) - pressBound);
                potentialDiffNw = (cellData.pressure(nPhaseIdx) - pressBound - pcBound);
                break;
            }
            case pn:
            {
                potentialDiffW = (cellData.pressure(wPhaseIdx) - pressBound + pcBound);
                potentialDiffNw = (cellData.pressure(nPhaseIdx) - pressBound);
                break;
            }
            case pGlobal:
            {
                potentialDiffW = (cellData.globalPressure() - pressBound - fMeanNw * (pcI - pcBound));
                potentialDiffNw = (cellData.globalPressure() - pressBound + fMeanW * (pcI - pcBound));
                break;
            }
            }

            potentialDiffW += density_[wPhaseIdx] * (distVec * gravity_);
            potentialDiffNw += density_[nPhaseIdx] * (distVec * gravity_);
        }

        // do the upwinding of the mobility depending on the phase potentials
        Scalar lambdaW = (potentialDiffW > 0.) ? lambdaWI : lambdaWBound;
        lambdaW = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffW, 0.0, 1.0e-30)) ? 0.5 * (lambdaWI + lambdaWBound) : lambdaW;
        Scalar lambdaNw = (potentialDiffNw > 0.) ? lambdaNwI : lambdaNwBound;
        lambdaNw = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNw, 0.0, 1.0e-30)) ? 0.5 * (lambdaNwI + lambdaNwBound) : lambdaNw;

        if (compressibility_)
        {
            density_[wPhaseIdx] = (potentialDiffW > 0.) ? cellData.density(wPhaseIdx) : densityWBound;
            density_[wPhaseIdx] = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffW, 0.0, 1.0e-30)) ? rhoMeanW : density_[wPhaseIdx];
            density_[nPhaseIdx] = (potentialDiffNw > 0.) ? cellData.density(nPhaseIdx) : densityNwBound;
            density_[nPhaseIdx] = (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potentialDiffNw, 0.0, 1.0e-30)) ? rhoMeanNw : density_[nPhaseIdx];
        }

        Scalar scalarPerm = permeability.two_norm();
        // calculate current matrix entry
        entry[matrix] = (lambdaW + lambdaNw) * scalarPerm / dist * faceArea;
        entry[rhs] = entry[matrix] * pressBound;

        // calculate right hand side
        // calculate unit distVec
        distVec /= dist;
        Scalar areaScaling = (unitOuterNormal * distVec);
        // this treatment of g allows to account for gravity flux through faces where the face normal has no z component (e.g. parallelepiped grids)
        entry[rhs] -= (lambdaW * density_[wPhaseIdx] + lambdaNw * density_[nPhaseIdx]) * scalarPerm * (gravity_ * distVec)
                * faceArea * areaScaling;

        if (pressureType_ == pw)
        {
            //add capillary pressure term to right hand side
            entry[rhs] -= 0.5 * (lambdaNwI + lambdaNwBound) * scalarPerm * (pcI - pcBound) / dist * faceArea;
        }
        else if (pressureType_ ==  pn)
        {
            //add capillary pressure term to right hand side
            entry[rhs] += 0.5 * (lambdaWI + lambdaWBound) * scalarPerm * (pcI - pcBound) / dist * faceArea;
        }
    }
    // set Neumann boundary condition
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
}

/*!
 * \brief Updates constitutive relations and stores them in the variable class
 *
 * Stores mobility, fractional flow function and capillary pressure for all grid cells.
 * In the compressible case additionally the densities and viscosities are stored.
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
void FVPressure2P<TypeTag>::updateMaterialLaws()
{
    for (const auto& element : elements(problem_.gridView()))
    {
        int eIdxGlobal = problem_.variables().index(element);

        CellData& cellData = problem_.variables().cellData(eIdxGlobal);

        const Scalar temperature = problem_.temperature(element);

        // determine phase saturations from primary saturation variable
        const Scalar satW = cellData.saturation(wPhaseIdx);

        // old material law interface is deprecated: Replace this by
        // const auto& fluidMatrixInteraction = spatialParams.fluidMatrixInteractionAtPos(element.geometry().center());
        // after the release of 3.3, when the deprecated interface is no longer supported
        const auto fluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem_.spatialParams(), element);

        const Scalar pc = fluidMatrixInteraction.pc(satW);

        // determine phase pressures from primary pressure variable
        Scalar pressW = 0;
        Scalar pressNw = 0;
        if (pressureType_ == pw)
        {
            pressW = cellData.pressure(wPhaseIdx);
            pressNw = pressW + pc;
        }
        else if (pressureType_ == pn)
        {
            pressNw = cellData.pressure(nPhaseIdx);
            pressW = pressNw - pc;
        }

        if (compressibility_)
        {
            FluidState& fluidState = cellData.fluidState();
            fluidState.setTemperature(temperature);

            fluidState.setPressure(wPhaseIdx, pressW);
            fluidState.setPressure(nPhaseIdx, pressNw);

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

            if (pressureType_ != pGlobal)
            {
                cellData.setPressure(wPhaseIdx, pressW);
                cellData.setPressure(nPhaseIdx, pressNw);
            }
        }

        // initialize mobilities
        Scalar mobilityW = fluidMatrixInteraction.krw(satW) / viscosity_[wPhaseIdx];
        Scalar mobilityNw = fluidMatrixInteraction.krn(satW) / viscosity_[nPhaseIdx];

        if (compressibility_)
        {
            mobilityW *= density_[wPhaseIdx];
            mobilityNw *= density_[nPhaseIdx];
        }

        // initialize mobilities
        cellData.setMobility(wPhaseIdx, mobilityW);
        cellData.setMobility(nPhaseIdx, mobilityNw);

        // initialize fractional flow functions
        cellData.setFracFlowFunc(wPhaseIdx, mobilityW / (mobilityW + mobilityNw));
        cellData.setFracFlowFunc(nPhaseIdx, mobilityNw / (mobilityW + mobilityNw));

        const Scalar gravityDiff = (problem_.bBoxMax() - element.geometry().center()) * gravity_;

        Scalar potW = pressW + gravityDiff * density_[wPhaseIdx];
        Scalar potNw = pressNw + gravityDiff * density_[nPhaseIdx];

        if (pressureType_ == pGlobal)
        {
            potW = cellData.globalPressure() - cellData.fracFlowFunc(nPhaseIdx) * pc + gravityDiff * density_[wPhaseIdx];
            potNw = cellData.globalPressure() - cellData.fracFlowFunc(wPhaseIdx) * pc + gravityDiff * density_[nPhaseIdx];
        }

        cellData.setPotential(wPhaseIdx, potW);
        cellData.setPotential(nPhaseIdx, potNw);
    }
}

} // end namespace Dumux
#endif
