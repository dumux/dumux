// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 *
 * \brief A test problem for the one-phase pore network model.
 */
#ifndef DUMUX_PNM2P_PROBLEM_HH
#define DUMUX_PNM2P_PROBLEM_HH

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porenetwork/2p/model.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/porenetwork/common/utilities.hh>

namespace Dumux
{
template <class TypeTag>
class DrainageProblem;

namespace Properties
{
// Create new type tags
namespace TTag {
struct DrainageProblem { using InheritsFrom = std::tuple<PNMTwoP>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DrainageProblem> { using type = Dumux::DrainageProblem<TypeTag>; };

template<class TypeTag>
struct Formulation<TypeTag, TTag::DrainageProblem>
{ static constexpr auto value = TwoPFormulation::p0s1; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DrainageProblem>
 {
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = Dumux::FluidSystems::H2OAir<Scalar, Dumux::Components::SimpleH2O<Scalar>>;
 };

 // Set the grid type
 template<class TypeTag>
 struct Grid<TypeTag, TTag::DrainageProblem> { using type = Dune::FoamGrid<1, 3>; };

}

template <class TypeTag>
class DrainageProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Labels = GetPropType<TypeTag, Properties::Labels>;
    enum {
        pwIdx = Indices::pressureIdx,
        snIdx = Indices::saturationIdx,
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<GridView::dimension>::Entity;

public:
    template<class SpatialParams>
    DrainageProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<SpatialParams> spatialParams)
    : ParentType(gridGeometry, spatialParams)
    {
        vtpOutputFrequency_ = getParam<int>("Problem.VtpOutputFrequency");
        initialPc_ = getParam<Scalar>("Problem.InitialPc");
        finalPc_ = getParam<Scalar>("Problem.FinalPc");
        numSteps_ = getParam<int>("Problem.NumSteps");
        swShiftThreshold_ = getParam<Scalar>("Problem.RelShiftThreshold", 1e-12);
        writeOnlyEqPoints_ = getParam<bool>("Problem.WriteOnlyEquilibriumPoints", false);

        pcEpisopde_.resize(numSteps_ + 1);
        for (int i = 0 ; i < pcEpisopde_.size(); i++)
              pcEpisopde_[i] = initialPc_ + i*(finalPc_ - initialPc_)/numSteps_;

        std::cout << "The following global PCs are applied: " << std::endl;
        for (auto x: pcEpisopde_)
        {
            std::cout << x << std::endl;
        }

        if (!writeOnlyEqPoints_)
        {
            logfile_.open("logfile_" + this->name() + ".txt"); //for the logfile
            logfile_ <<"Logfile for: " + this->name()  << std::endl;
            logfile_ << std::left << std::setw(20) << std::setfill(' ') << "Time"
                     << std::left << std::setw(20) << std::setfill(' ') << "globalPc"
                     << std::left << std::setw(20) << std::setfill(' ') << "swAveraged"
                     << std::left << std::setw(20) << std::setfill(' ') << "pwAveraged"
                     << std::left << std::setw(20) << std::setfill(' ') << "pnAveraged"
                     << std::left << std::setw(20) << std::setfill(' ') << "pcAveraged"
                     << std::left << std::setw(20) << std::setfill(' ') << "numThroatsInvaded"
                     << std::endl;
        }

        logfileEqPoints_.open("eqPoints_" + this->name() + ".txt");
        logfileEqPoints_ << std::left << std::setw(20) << std::setfill(' ') << "Time"
                         << std::left << std::setw(20) << std::setfill(' ') << "globalPc"
                         << std::left << std::setw(20) << std::setfill(' ') << "swAveraged"
                         << std::left << std::setw(20) << std::setfill(' ') << "pwAveraged"
                         << std::left << std::setw(20) << std::setfill(' ') << "pnAveraged"
                         << std::left << std::setw(20) << std::setfill(' ') << "pcAveraged"
                         << std::left << std::setw(20) << std::setfill(' ') << "numThroatsInvaded"
                         << std::endl;
        step_ = 0;
    }

    /*!
     * \name Simulation steering
     */
    // \{

    /*!
     * \name Problem parameters
     */
    // \{

    bool shouldWriteOutput(const int timeStepIndex, const GridVariables& gridVariables) const
    {
        const  bool wasInEquilibrium = inEquilibrium_;
        inEquilibrium_ = false;
        if (vtpOutputFrequency_ < 0)
            return true;

        if (vtpOutputFrequency_ == 0)
            return (timeStepIndex == 0 || gridVariables.gridFluxVarsCache().invasionState().hasChanged() || wasInEquilibrium);
        else
            return (timeStepIndex % vtpOutputFrequency_ == 0 || gridVariables.gridFluxVarsCache().invasionState().hasChanged() || wasInEquilibrium);
    }


    /*!
     * \brief Called at the end of each time step
     */
    template<class AveragedValues>
    void postTimeStep(const Scalar time, const AveragedValues& avgValues, std::size_t numThroatsInvaded, const Scalar dt)
    {
        const Scalar avgSw = avgValues["avgSat"];

        if (!writeOnlyEqPoints_)
        {
            logfile_ << std::fixed << std::left << std::setw(20) << std::setfill(' ') << time
                                   << std::left << std::setw(20) << std::setfill(' ') <<  pcEpisopde_[step_]
                                   << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgSat"]
                                   << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgPw"]
                                   << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgPn"]
                                   << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgPn"] - avgValues["avgPw"]
                                   << std::left << std::setw(20) << std::setfill(' ') << numThroatsInvaded
                                   << std::endl;
        }

        // store the three most recent averaged saturations
        std::rotate(swAvg_.rbegin(), swAvg_.rbegin()+1, swAvg_.rend());
        swAvg_[0]= avgSw;

        // Check for steady state and end episode
        dSwDt_ = std::abs(swAvg_[0]-swAvg_[1])/dt;
        const Scalar pc = pcEpisopde_[step_];
        std::cout << "global pC applied: " << pc << " / " << finalPc_ << " (step " << step_ << " of " << numSteps_ << ")" << std::endl;
        std::cout << "swAverage: " << swAvg_[0] << " (relative shift: " << dSwDt_ << "). " << std::endl;
        std::cout << numThroatsInvaded << " of " << this->gridGeometry().gridView().size(0) << " throats invaded." << std::endl;
        if(dSwDt_ < swShiftThreshold_)
        {
            std::cout << "Equlibrium point reached!" << std::endl;

            logfileEqPoints_ << std::fixed << std::left << std::setw(20) << std::setfill(' ') << time
                                           << std::left << std::setw(20) << std::setfill(' ') <<  pcEpisopde_[step_]
                                           << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgSat"]
                                           << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgPw"]
                                           << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgPn"]
                                           << std::left << std::setw(20) << std::setfill(' ') << avgValues["avgPn"] - avgValues["avgPw"]
                                           << std::left << std::setw(20) << std::setfill(' ') << numThroatsInvaded
                                           << std::endl;
            inEquilibrium_ = true;
            ++step_;
        }
    }

    bool equilibriumPointReached() const
    { return dSwDt_ < swShiftThreshold_; }

     /*!
     * \name Boundary conditions
     */
    // \{
    //! Specifies which kind of boundary condition should be used for
    //! which equation for a finite volume on the boundary.
    BoundaryTypes boundaryTypes(const Element &element, const SubControlVolume &scv) const
    {
        BoundaryTypes bcTypes;
        if (isInletPore_(scv) || isOutletPore_(scv))
           bcTypes.setAllDirichlet();

        else // neuman for the remaining boundaries
           bcTypes.setAllNeumann();

        return bcTypes;
    }

    /*!
    * \brief Evaluate the boundary conditions for a dirichlet
    *        control volume.
    *
    * \param values The dirichlet values for the primary variables
    * \param vertex The vertex (pore body) for which the condition is evaluated
    *
    * For this method, the \a values parameter stores primary variables.
    */
    PrimaryVariables dirichlet(const Element &element,
                               const SubControlVolume &scv) const
    {
        PrimaryVariables values(0.0);
        values[pwIdx] = 1e5;
        values[snIdx] = 0.0;

        const auto& fluidMatrixInteraction = this->spatialParams().fluidMatrixInteraction(element, scv, 0);
        if (isInletPore_(scv))
            values[snIdx] = 1.0 - fluidMatrixInteraction.sw(pcEpisopde_[step_]);

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{
    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub control volume
     *
     * For this method, the return parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
    PrimaryVariables source(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume &scv) const
    {
        return PrimaryVariables(0.0);
    }
    // \}

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initial(const Vertex& vertex) const
    {
        PrimaryVariables values(0.0);
        values[pwIdx] = 1e5;
        return values;
    }

    /*!
     * \brief Evaluate the initial invasion state of a pore throat
     *
     * Returns true for a invaded throat and false elsewise.
     */
    bool initialInvasionState(const Element &element) const
    { return false; }

    // \}

    /*!
     * \brief Calls gnuplot and plots the pc-S curve
     */
    void plotPcS()
    {
        std::FILE * pipe_;
        pipe_ = popen("gnuplot -persist", "w");
        std::string command = "set xrange [0:1] \n";
        command += "set xlabel 'S_w' \n";
        command += "set ylabel 'p_c' \n";
        std::string filename = "'logfile_"+ this->name() +".txt'";
        command += "plot " + filename + " using 4:3 with lines title 'p_c(S_w)'";
        fputs((command + "\n").c_str(), pipe_);
        pclose(pipe_);
    }

    bool simulationFinished() const
    { return (step_ > numSteps_) && (dSwDt_ < swShiftThreshold_) ; }

    /*!
     * \brief Imposes a global capillary pressure for the invasion mechanism. Only throats
     *        with an entry pressure smaller than this can be invaded.
     */
    Scalar globalCapillaryPressure() const
    { return pcEpisopde_[step_]; }

    Scalar dSwDt() const
    { return dSwDt_; }

private:

    bool isInletPore_(const SubControlVolume& scv) const
    {
        return isInletPore_(scv.dofIndex());
    }

    bool isInletPore_(const std::size_t dofIdxGlobal) const
    {
        return this->gridGeometry().poreLabel(dofIdxGlobal) == Labels::inlet;
    }

    bool isOutletPore_(const SubControlVolume& scv) const
    {
        return this->gridGeometry().poreLabel(scv.dofIndex()) == Labels::outlet;
    }

    int vtpOutputFrequency_;
    Scalar initialPc_;
    Scalar finalPc_;
    int numSteps_;
    std::vector<Scalar> pcEpisopde_;
    std::array<Scalar, 2> swAvg_ = {{1.0, 1.0}};
    std::ofstream logfile_;
    std::ofstream logfileEqPoints_;
    Scalar dSwDt_ = 0.0;
    mutable bool inEquilibrium_ = false;
    Scalar swShiftThreshold_;
    bool writeOnlyEqPoints_;

    int step_;
};
} //end namespace

#endif
