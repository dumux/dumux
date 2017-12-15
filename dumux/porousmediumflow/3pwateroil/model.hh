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
 * \brief Adaption of the fully implicit scheme to the three-phase three-component flow model.
 *
 * The model is designed for simulating three fluid phases with water, gas, and
 * a liquid contaminant (NAPL - non-aqueous phase liquid)
 */
#ifndef DUMUX_3PWATEROIL_MODEL_HH
#define DUMUX_3PWATEROIL_MODEL_HH

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup ThreePWaterOilModel
 * \brief Adaption of the fully implicit scheme to the three-phase two-component flow model.
 *
 * This model implements three-phase two-component flow of three fluid phases
 * \f$\alpha \in \{ water, gas, NAPL \}\f$ each composed of up to two components
 * \f$\kappa \in \{ water, contaminant \}\f$. The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum:
 * \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one transport equation for each component is obtained as
 * \f{eqnarray*}
 && \phi \frac{\partial (\sum_\alpha \varrho_\alpha X_\alpha^\kappa
 S_\alpha )}{\partial t}
 - \sum\limits_\alpha \text{div} \left\{ \frac{k_{r\alpha}}{\mu_\alpha}
 \varrho_\alpha x_\alpha^\kappa \mbox{\bf K}
 (\textbf{grad}\, p_\alpha - \varrho_\alpha \mbox{\bf g}) \right\}
 \nonumber \\
 \nonumber \\
 && - \sum\limits_\alpha \text{div} \left\{ D_\text{pm}^\kappa \varrho_\alpha \frac{M^\kappa}{M_\alpha}
 \textbf{grad} x^\kappa_{\alpha} \right\}
 - q^\kappa = 0 \qquad \forall \kappa , \; \forall \alpha
 \f}
 *
 * Note that these balance equations are molar.
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 *
 * The model uses commonly applied auxiliary conditions like
 * \f$S_w + S_n + S_g = 1\f$ for the saturations and
 * \f$x^w_\alpha + x^c_\alpha = 1\f$ for the mole fractions.
 * Furthermore, the phase pressures are related to each other via
 * capillary pressures between the fluid phases, which are functions of
 * the saturation, e.g. according to the approach of Parker et al.
 *
 * The used primary variables are dependent on the locally present fluid phases
 * An adaptive primary variable switch is included. The phase state is stored for all nodes
 * of the system. Different cases can be distinguished:
 * <ul>
 *  <li> ... to be completed ... </li>
 * </ul>
 */
template<class TypeTag>
class ThreePWaterOilModel: public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseModel) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        pressureIdx = Indices::pressureIdx,
        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,


        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,

        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx,

        threePhases = Indices::threePhases,
        wPhaseOnly  = Indices::wPhaseOnly,
        gnPhaseOnly = Indices::gnPhaseOnly,
        wnPhaseOnly = Indices::wnPhaseOnly,
        gPhaseOnly  = Indices::gPhaseOnly,
        wgPhaseOnly = Indices::wgPhaseOnly

    };


    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<0>::Entity Element;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);


    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    /*!
     * \brief Initialize the static data with the initial solution.
     *
     * \param problem The problem to be solved
     */
    void init(Problem &problem)
    {
        ParentType::init(problem);

        staticDat_.resize(this->numDofs());

        setSwitched_(false);

        if (isBox)
        {
            for (const auto& vertex : vertices(this->gridView_()))
            {
                int globalIdx = this->dofMapper().index(vertex);
                const GlobalPosition &globalPos = vertex.geometry().corner(0);

                // initialize phase presence
                staticDat_[globalIdx].phasePresence
                    = this->problem_().initialPhasePresence(vertex, globalIdx,
                                                        globalPos);
                staticDat_[globalIdx].wasSwitched = false;

                staticDat_[globalIdx].oldPhasePresence
                    = staticDat_[globalIdx].phasePresence;
            }
        }
        else
        {
            for (const auto& element : elements(this->gridView_()))
            {
                int globalIdx = this->dofMapper().index(element);
                const GlobalPosition &globalPos = element.geometry().center();

                // initialize phase presence
                staticDat_[globalIdx].phasePresence
                    = this->problem_().initialPhasePresence(*this->gridView_().template begin<dim> (),
                                                            globalIdx, globalPos);
                staticDat_[globalIdx].wasSwitched = false;

                staticDat_[globalIdx].oldPhasePresence
                    = staticDat_[globalIdx].phasePresence;
            }
        }
    }

    /*!
     * \brief Compute the total storage inside one phase of all
     *        conservation quantities.
     *
     * \param storage Contains the storage of each component for one phase
     * \param phaseIdx The phase index
     */
    void globalPhaseStorage(PrimaryVariables &storage, const int phaseIdx)
    {
        storage = 0;

        for (const auto& element : elements(this->gridView_(), Dune::Partitions::interior))
        {
            this->localResidual().evalPhaseStorage(element, phaseIdx);

            for (unsigned int i = 0; i < this->localResidual().storageTerm().size(); ++i)
                storage += this->localResidual().storageTerm()[i];
        }
        if (this->gridView_().comm().size() > 1)
            storage = this->gridView_().comm().sum(storage);
    }


    /*!
     * \brief Called by the update() method if applying the newton
     *         method was unsuccessful.
     */
    void updateFailed()
    {
        ParentType::updateFailed();

        setSwitched_(false);
        resetPhasePresence_();
    };

    /*!
     * \brief Called by the problem if a time integration was
     *        successful, post processing of the solution is done and the
     *        result has been written to disk.
     *
     * This should prepare the model for the next time integration.
     */
    void advanceTimeLevel()
    {
        ParentType::advanceTimeLevel();

        // update the phase state
        updateOldPhasePresence_();
        setSwitched_(false);
    }

    /*!
     * \brief Return true if the primary variables were switched for
     *        at least one vertex after the last timestep.
     */
    bool switched() const
    {
        return switchFlag_;
    }

    /*!
     * \brief Returns the phase presence of the current or the old solution of a degree of freedom.
     *
     * \param globalIdx The global index of the degree of freedom
     * \param oldSol Evaluate function with solution of current or previous time step
     */
    int phasePresence(int globalIdx, bool oldSol) const
    {
        return
            oldSol
            ? staticDat_[globalIdx].oldPhasePresence
            : staticDat_[globalIdx].phasePresence;
    }

    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     *
     * \param sol The solution vector
     * \param writer The writer for multi-file VTK datasets
     */
    template<class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol,
                            MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<double, dim> > VectorField;

        // get the number of degrees of freedom
        unsigned numDofs = this->numDofs();

        // create the required scalar fields
        ScalarField *saturation[numPhases];
        ScalarField *pressure[numPhases];
        ScalarField *density[numPhases];

        ScalarField *mobility[numPhases];
        ScalarField *viscosity[numPhases];


        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            saturation[phaseIdx] = writer.allocateManagedBuffer(numDofs);
            pressure[phaseIdx] = writer.allocateManagedBuffer(numDofs);
            density[phaseIdx] = writer.allocateManagedBuffer(numDofs);

            mobility[phaseIdx] = writer.allocateManagedBuffer(numDofs);
            viscosity[phaseIdx] = writer.allocateManagedBuffer(numDofs);
        }

        ScalarField *phasePresence = writer.allocateManagedBuffer (numDofs);
        ScalarField *moleFraction[numPhases][numComponents];
        for (int i = 0; i < numPhases; ++i)
            for (int j = 0; j < numComponents; ++j)
                moleFraction[i][j] = writer.allocateManagedBuffer (numDofs);
        ScalarField *temperature = writer.allocateManagedBuffer (numDofs);
        ScalarField *poro = writer.allocateManagedBuffer(numDofs);
        ScalarField *perm = writer.allocateManagedBuffer(numDofs);
        VectorField *velocityN = writer.template allocateManagedBuffer<double, dim>(numDofs);
        VectorField *velocityW = writer.template allocateManagedBuffer<double, dim>(numDofs);
        VectorField *velocityG = writer.template allocateManagedBuffer<double, dim>(numDofs);
        ImplicitVelocityOutput<TypeTag> velocityOutput(this->problem_());

        if (velocityOutput.enableOutput()) // check if velocity output is demanded
        {
            // initialize velocity fields
            for (unsigned int i = 0; i < numDofs; ++i)
            {
                (*velocityN)[i] = Scalar(0);
                (*velocityW)[i] = Scalar(0);
                (*velocityG)[i] = Scalar(0);
            }
        }

        unsigned numElements = this->gridView_().size(0);
        ScalarField *rank = writer.allocateManagedBuffer (numElements);

        for (const auto& element : elements(this->gridView_(), Dune::Partitions::interior))
        {
            int idx = this->problem_().elementMapper().index(element);
            (*rank)[idx] = this->gridView_().comm().rank();

            FVElementGeometry fvGeometry;
            fvGeometry.update(this->gridView_(), element);


            ElementVolumeVariables elemVolVars;
            elemVolVars.update(this->problem_(),
                               element,
                               fvGeometry,
                               false /* oldSol? */);

            for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
            {
                int globalIdx = this->dofMapper().subIndex(element, scvIdx, dofCodim);

                for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                    (*saturation[phaseIdx])[globalIdx] = elemVolVars[scvIdx].fluidState().saturation(phaseIdx);
                    (*pressure[phaseIdx])[globalIdx] = elemVolVars[scvIdx].fluidState().pressure(phaseIdx);
                    (*density[phaseIdx])[globalIdx] = elemVolVars[scvIdx].fluidState().density(phaseIdx);

                    (*mobility[phaseIdx])[globalIdx] = elemVolVars[scvIdx].mobility(phaseIdx);
                    (*viscosity[phaseIdx])[globalIdx] = elemVolVars[scvIdx].fluidState().viscosity(phaseIdx);
                }

                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                        (*moleFraction[phaseIdx][compIdx])[globalIdx] =
                            elemVolVars[scvIdx].fluidState().moleFraction(phaseIdx,
                                                              compIdx);

                        Valgrind::CheckDefined((*moleFraction[phaseIdx][compIdx])[globalIdx]);
                    }
                }

                (*poro)[globalIdx] = elemVolVars[scvIdx].porosity();
                (*perm)[globalIdx] = elemVolVars[scvIdx].permeability();
                (*temperature)[globalIdx] = elemVolVars[scvIdx].temperature();
                (*phasePresence)[globalIdx] = staticDat_[globalIdx].phasePresence;
            }

            // velocity output
            velocityOutput.calculateVelocity(*velocityW, elemVolVars, fvGeometry, element, wPhaseIdx);
            velocityOutput.calculateVelocity(*velocityN, elemVolVars, fvGeometry, element, nPhaseIdx);
            velocityOutput.calculateVelocity(*velocityN, elemVolVars, fvGeometry, element, gPhaseIdx);

        }

        writer.attachDofData(*saturation[wPhaseIdx], "sw", isBox);
        writer.attachDofData(*saturation[nPhaseIdx], "sn", isBox);
        writer.attachDofData(*saturation[gPhaseIdx], "sg", isBox);
        writer.attachDofData(*pressure[wPhaseIdx], "pw", isBox);
        writer.attachDofData(*pressure[nPhaseIdx], "pn", isBox);
        writer.attachDofData(*pressure[gPhaseIdx], "pg", isBox);
        writer.attachDofData(*density[wPhaseIdx], "rhow", isBox);
        writer.attachDofData(*density[nPhaseIdx], "rhon", isBox);
        writer.attachDofData(*density[gPhaseIdx], "rhog", isBox);

        writer.attachDofData(*mobility[wPhaseIdx], "MobW", isBox);
        writer.attachDofData(*mobility[nPhaseIdx], "MobN", isBox);
        writer.attachDofData(*mobility[gPhaseIdx], "MobG", isBox);

        writer.attachDofData(*viscosity[wPhaseIdx], "ViscosW", isBox);
        writer.attachDofData(*viscosity[nPhaseIdx], "ViscosN", isBox);
        writer.attachDofData(*viscosity[gPhaseIdx], "ViscosG", isBox);
        for (int i = 0; i < numPhases; ++i)
        {
            for (int j = 0; j < numComponents; ++j)
            {
                std::ostringstream oss;
                oss << "x^"
                    << FluidSystem::phaseName(i)
                    << "_"
                    << FluidSystem::componentName(j);
                writer.attachDofData(*moleFraction[i][j], oss.str().c_str(), isBox);
            }
        }
        writer.attachDofData(*poro, "porosity", isBox);
        writer.attachDofData(*perm, "permeability", isBox);
        writer.attachDofData(*temperature, "temperature", isBox);
        writer.attachDofData(*phasePresence, "phase presence", isBox);

        if (velocityOutput.enableOutput()) // check if velocity output is demanded
        {
            writer.attachDofData(*velocityW,  "velocityW", isBox, dim);
            writer.attachDofData(*velocityN,  "velocityN", isBox, dim);
            writer.attachDofData(*velocityG,  "velocityG", isBox, dim);
        }

        writer.attachCellData(*rank, "process rank");
    }

    /*!
     * \brief Write the current solution to a restart file.
     *
     * \param outStream The output stream of one entity for the restart file
     * \param entity The entity, either a vertex or an element
     */
    template<class Entity>
    void serializeEntity(std::ostream &outStream, const Entity &entity)
    {
        // write primary variables
        ParentType::serializeEntity(outStream, entity);

        int globalIdx = this->dofMapper().index(entity);
        if (!outStream.good())
            DUNE_THROW(Dune::IOError, "Could not serialize entity " << globalIdx);

        outStream << staticDat_[globalIdx].phasePresence << " ";
    }

    /*!
     * \brief Reads the current solution from a restart file.
     *
     * \param inStream The input stream of one entity from the restart file
     * \param entity The entity, either a vertex or an element
     */
    template<class Entity>
    void deserializeEntity(std::istream &inStream, const Entity &entity)
    {
        // read primary variables
        ParentType::deserializeEntity(inStream, entity);

        // read phase presence
        int globalIdx = this->dofMapper().index(entity);
        if (!inStream.good())
            DUNE_THROW(Dune::IOError,
                       "Could not deserialize entity " << globalIdx);

        inStream >> staticDat_[globalIdx].phasePresence;
        staticDat_[globalIdx].oldPhasePresence
            = staticDat_[globalIdx].phasePresence;

    }

    /*!
     * \brief Update the static data of all vertices in the grid.
     *
     * \param curGlobalSol The current global solution
     * \param oldGlobalSol The previous global solution
     */
    void updateStaticData(SolutionVector &curGlobalSol,
                          const SolutionVector &oldGlobalSol)
    {
        bool wasSwitched = false;

        bool useSimpleModel = GET_PROP_VALUE(TypeTag, UseSimpleModel);

        for (unsigned i = 0; i < staticDat_.size(); ++i)
            staticDat_[i].visited = false;

        FVElementGeometry fvGeometry;
        static VolumeVariables volVars;
        for (const auto& element : elements(this->gridView_()))
        {
            fvGeometry.update(this->gridView_(), element);
            for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
            {
                int globalIdx = this->dofMapper().subIndex(element, scvIdx, dofCodim);

                if (staticDat_[globalIdx].visited)
                    continue;

                staticDat_[globalIdx].visited = true;
                volVars.update(curGlobalSol[globalIdx],
                               this->problem_(),
                               element,
                               fvGeometry,
                               scvIdx,
                               false);
                const GlobalPosition &globalPos = fvGeometry.subContVol[scvIdx].global;

                if(useSimpleModel)
                {
                    if (primaryVarSwitchSimple_(curGlobalSol,
                                          volVars,
                                          globalIdx,
                                          globalPos))
                    {
                        this->jacobianAssembler().markDofRed(globalIdx);
                        wasSwitched = true;
                    }
                }
                else
                {
                    if (primaryVarSwitch_(curGlobalSol,
                                          volVars,
                                          globalIdx,
                                          globalPos))
                    {
                        this->jacobianAssembler().markDofRed(globalIdx);
                        wasSwitched = true;
                    }
                }
            }
        }

        // make sure that if there was a variable switch in an
        // other partition we will also set the switch flag
        // for our partition.
        if (this->gridView_().comm().size() > 1)
            wasSwitched = this->gridView_().comm().max(wasSwitched);

        setSwitched_(wasSwitched);
    }

protected:
    /*!
     * \brief Data which is attached to each vertex and is not only
     *        stored locally.
     */
    struct StaticVars
    {
        int phasePresence;
        bool wasSwitched;

        int oldPhasePresence;
        bool visited;
    };

    /*!
     * \brief Reset the current phase presence of all vertices to the old one.
     *
     * This is done after an update failed.
     */
    void resetPhasePresence_()
    {
        for (unsigned int i = 0; i < this->numDofs(); ++i)
        {
            staticDat_[i].phasePresence
                = staticDat_[i].oldPhasePresence;
            staticDat_[i].wasSwitched = false;
        }
    }

    /*!
     * \brief Set the old phase of all verts state to the current one.
     */
    void updateOldPhasePresence_()
    {
        for (unsigned int i = 0; i < this->numDofs(); ++i)
        {
            staticDat_[i].oldPhasePresence
                = staticDat_[i].phasePresence;
            staticDat_[i].wasSwitched = false;
        }
    }

    /*!
     * \brief Set whether there was a primary variable switch after in
     *        the last timestep.
     */
    void setSwitched_(bool yesno)
    {
        switchFlag_ = yesno;
    }

    //  perform variable switch at a vertex; Returns true if a
    //  variable switch was performed.
    bool primaryVarSwitch_(SolutionVector &globalSol,
                           const VolumeVariables &volVars,
                           int globalIdx,
                           const GlobalPosition &globalPos)
    {
        // evaluate primary variable switch
        bool wouldSwitch = false;
        int phasePresence = staticDat_[globalIdx].phasePresence;
        int newPhasePresence = phasePresence;

        // check if a primary var switch is necessary
        if (phasePresence == threePhases)
        {
            Scalar Smin = 0;
            if (staticDat_[globalIdx].wasSwitched)
                Smin = -0.01;

            if (volVars.saturation(gPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // gas phase disappears
                std::cout << "Gas phase disappears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", sg: "
                          << volVars.saturation(gPhaseIdx) << std::endl;
                newPhasePresence = wnPhaseOnly;

                globalSol[globalIdx][pressureIdx] = volVars.fluidState().pressure(wPhaseIdx);
                globalSol[globalIdx][switch1Idx] = volVars.fluidState().temperature();
            }
            else if (volVars.saturation(wPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // water phase disappears
                std::cout << "Water phase disappears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", sw: "
                          << volVars.saturation(wPhaseIdx) << std::endl;
                newPhasePresence = gnPhaseOnly;

                globalSol[globalIdx][switch1Idx] = volVars.fluidState().saturation(nPhaseIdx);
                globalSol[globalIdx][switch2Idx] = volVars.fluidState().moleFraction(nPhaseIdx, wCompIdx);
            }
            else if (volVars.saturation(nPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // NAPL phase disappears
                std::cout << "NAPL phase disappears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", sn: "
                          << volVars.saturation(nPhaseIdx) << std::endl;
                newPhasePresence = wgPhaseOnly;

                globalSol[globalIdx][switch2Idx] = volVars.fluidState().moleFraction(wPhaseIdx, nCompIdx);
            }
        }
        else if (phasePresence == wPhaseOnly)
        {
            bool gasFlag = 0;
            bool nonwettingFlag = 0;
            // calculate fractions of the partial pressures in the
            // hypothetical gas phase
            Scalar xwg = volVars.fluidState().moleFraction(gPhaseIdx, wCompIdx);
            Scalar xng = volVars.fluidState().moleFraction(gPhaseIdx, nCompIdx);

            Scalar xgMax = 1.0;
            if (xwg + xng > xgMax)
                wouldSwitch = true;
            if (staticDat_[globalIdx].wasSwitched)
                xgMax *= 1.02;

            // if the sum of the mole fractions would be larger than
            // 100%, gas phase appears
            if (xwg + xng > xgMax)
            {
                // gas phase appears
                std::cout << "gas phase appears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", xwg + xng: "
                          << xwg + xng << std::endl;
                gasFlag = 1;
            }

            // calculate fractions in the hypothetical NAPL phase
            Scalar xnn = volVars.fluidState().moleFraction(nPhaseIdx, nCompIdx);
            /* take care:
               for xnn in case wPhaseOnly we compute xnn=henry_mesitylene*x1w,
               where a hypothetical gas pressure is assumed for the Henry
               x0n is set to NULL  (all NAPL phase is dirty)
               x2n is set to NULL  (all NAPL phase is dirty)
            */

            Scalar xnMax = 1.0;
            if (xnn > xnMax)
                wouldSwitch = true;
            if (staticDat_[globalIdx].wasSwitched)
                xnMax *= 1.02;

            // if the sum of the hypothetical mole fractions would be larger than
            // 100%, NAPL phase appears
            if (xnn > xnMax)
            {
                // NAPL phase appears
                std::cout << "NAPL phase appears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", xnn: "
                          << xnn << std::endl;
                nonwettingFlag = 1;
            }

            if ((gasFlag == 1) && (nonwettingFlag == 0))
            {
                newPhasePresence = wgPhaseOnly;
                globalSol[globalIdx][switch1Idx] = 0.9999;
                globalSol[globalIdx][pressureIdx] = volVars.fluidState().pressure(gPhaseIdx);
            }
            else if ((gasFlag == 1) && (nonwettingFlag == 1))
            {
                newPhasePresence = threePhases;
                globalSol[globalIdx][pressureIdx] = volVars.fluidState().pressure(gPhaseIdx);
                globalSol[globalIdx][switch1Idx] = 0.999;
            }
            else if ((gasFlag == 0) && (nonwettingFlag == 1))
            {
                newPhasePresence = wnPhaseOnly;
                globalSol[globalIdx][switch2Idx] = 0.0001;
            }
        }
        else if (phasePresence == gnPhaseOnly)
        {
            bool nonwettingFlag = 0;
            bool wettingFlag = 0;

            Scalar Smin = 0.0;
            if (staticDat_[globalIdx].wasSwitched)
                Smin = -0.01;

            if (volVars.saturation(nPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // NAPL phase disappears
                std::cout << "NAPL phase disappears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", sn: "
                          << volVars.saturation(nPhaseIdx) << std::endl;
                nonwettingFlag = 1;
            }


            // calculate fractions of the hypothetical water phase
            Scalar xww = volVars.fluidState().moleFraction(wPhaseIdx, wCompIdx);
            /*
              take care:, xww, if no water is present, then take xww=xwg*pg/pwsat .
              If this is larger than 1, then water appears
            */
            Scalar xwMax = 1.0;
            if (xww > xwMax)
                wouldSwitch = true;
            if (staticDat_[globalIdx].wasSwitched)
                xwMax *= 1.02;

            // if the sum of the mole fractions would be larger than
            // 100%, water phase appears
            if (xww > xwMax)
            {
                // water phase appears
                std::cout << "water phase appears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", xww=xwg*pg/pwsat : "
                          << xww << std::endl;
                wettingFlag = 1;
            }

            if ((wettingFlag == 1) && (nonwettingFlag == 0))
            {
                newPhasePresence = threePhases;
                globalSol[globalIdx][switch1Idx] = 0.0001;
                globalSol[globalIdx][switch2Idx] = volVars.saturation(nPhaseIdx);
            }
            else if ((wettingFlag == 1) && (nonwettingFlag == 1))
            {
                newPhasePresence = wgPhaseOnly;
                globalSol[globalIdx][switch1Idx] = 0.0001;
                globalSol[globalIdx][switch2Idx] = volVars.fluidState().moleFraction(wPhaseIdx, nCompIdx);
            }
            else if ((wettingFlag == 0) && (nonwettingFlag == 1))
            {
                newPhasePresence = gPhaseOnly;
                globalSol[globalIdx][switch1Idx] = volVars.fluidState().temperature();
                globalSol[globalIdx][switch2Idx] = volVars.fluidState().moleFraction(gPhaseIdx, nCompIdx);
            }
        }
        else if (phasePresence == wnPhaseOnly)
        {
            bool nonwettingFlag = 0;
            bool gasFlag = 0;

            Scalar Smin = 0.0;
            if (staticDat_[globalIdx].wasSwitched)
                Smin = -0.01;

            if (volVars.saturation(nPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // NAPL phase disappears
                std::cout << "NAPL phase disappears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", sn: "
                          << volVars.saturation(nPhaseIdx) << std::endl;
                nonwettingFlag = 1;
            }

            // calculate fractions of the partial pressures in the
            // hypothetical gas phase
            Scalar xwg = volVars.fluidState().moleFraction(gPhaseIdx, wCompIdx);
            Scalar xng = volVars.fluidState().moleFraction(gPhaseIdx, nCompIdx);

            Scalar xgMax = 1.0;
            if (xwg + xng > xgMax)
                wouldSwitch = true;
            if (staticDat_[globalIdx].wasSwitched)
                xgMax *= 1.02;

            // if the sum of the mole fractions would be larger than
            // 100%, gas phase appears
            if (xwg + xng > xgMax)
            {
                // gas phase appears
                std::cout << "gas phase appears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", xwg + xng: "
                          << xwg + xng << std::endl;
                gasFlag = 1;
            }

            if ((gasFlag == 1) && (nonwettingFlag == 0))
            {
                newPhasePresence = threePhases;
                globalSol[globalIdx][pressureIdx] = volVars.pressure(gPhaseIdx);
                globalSol[globalIdx][switch1Idx] = volVars.saturation(wPhaseIdx);
            }
            else if ((gasFlag == 1) && (nonwettingFlag == 1))
            {
                newPhasePresence = wgPhaseOnly;
                globalSol[globalIdx][pressureIdx] = volVars.pressure(gPhaseIdx);
                globalSol[globalIdx][switch1Idx] = volVars.temperature();
                globalSol[globalIdx][switch2Idx] = volVars.fluidState().moleFraction(wPhaseIdx, nCompIdx);
            }
            else if ((gasFlag == 0) && (nonwettingFlag == 1))
            {
                newPhasePresence = wPhaseOnly;
                globalSol[globalIdx][switch2Idx] = volVars.fluidState().moleFraction(wPhaseIdx, nCompIdx);
            }
        }
        else if (phasePresence == gPhaseOnly)
        {
            bool nonwettingFlag = 0;
            bool wettingFlag = 0;

            // calculate fractions in the hypothetical NAPL phase
            Scalar xnn = volVars.fluidState().moleFraction(nPhaseIdx, nCompIdx);
            /*
              take care:, xnn, if no NAPL phase is there, take xnn=xng*pg/pcsat
              if this is larger than 1, then NAPL appears
            */

            Scalar xnMax = 1.0;
            if (xnn > xnMax)
                wouldSwitch = true;
            if (staticDat_[globalIdx].wasSwitched)
                xnMax *= 1.02;

            // if the sum of the hypothetical mole fraction would be larger than
            // 100%, NAPL phase appears
            if (xnn > xnMax)
            {
                // NAPL phase appears
                std::cout << "NAPL phase appears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", xnn: "
                          << xnn << std::endl;
                nonwettingFlag = 1;
            }
            // calculate fractions of the hypothetical water phase
            Scalar xww = volVars.fluidState().moleFraction(wPhaseIdx, wCompIdx);
            /*
              take care:, xww, if no water is present, then take xww=xwg*pg/pwsat .
              If this is larger than 1, then water appears
            */
            Scalar xwMax = 1.0;
            if (xww > xwMax)
                wouldSwitch = true;
            if (staticDat_[globalIdx].wasSwitched)
                xwMax *= 1.02;

            // if the sum of the mole fractions would be larger than
            // 100%, gas phase appears
            if (xww > xwMax)
            {
                // water phase appears
                std::cout << "water phase appears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", xww=xwg*pg/pwsat : "
                          << xww << std::endl;
                wettingFlag = 1;
            }
            if ((wettingFlag == 1) && (nonwettingFlag == 0))
            {
                newPhasePresence = wgPhaseOnly;
                globalSol[globalIdx][switch1Idx] = 0.0001;
                globalSol[globalIdx][switch2Idx] = volVars.fluidState().moleFraction(wPhaseIdx, nCompIdx);
            }
            else if ((wettingFlag == 1) && (nonwettingFlag == 1))
            {
                newPhasePresence = threePhases;
                globalSol[globalIdx][switch1Idx] = 0.0001;
                globalSol[globalIdx][switch2Idx] = 0.0001;
            }
            else if ((wettingFlag == 0) && (nonwettingFlag == 1))
            {
                newPhasePresence = gnPhaseOnly;
                globalSol[globalIdx][switch2Idx]
                    = volVars.fluidState().moleFraction(wPhaseIdx, nCompIdx);
            }
        }
        else if (phasePresence == wgPhaseOnly)
        {
            bool nonwettingFlag = 0;
            bool gasFlag = 0;
            bool wettingFlag = 0;

            // get the fractions in the hypothetical NAPL phase
            Scalar xnn = volVars.fluidState().moleFraction(nPhaseIdx, nCompIdx);

            // take care: if the NAPL phase is not present, take
            // xnn=xng*pg/pcsat if this is larger than 1, then NAPL
            // appears
            Scalar xnMax = 1.0;
            if (xnn > xnMax)
                wouldSwitch = true;
            if (staticDat_[globalIdx].wasSwitched)
                xnMax *= 1.02;

            // if the sum of the hypothetical mole fraction would be larger than
            // 100%, NAPL phase appears
            if (xnn > xnMax)
            {
                // NAPL phase appears
                std::cout << "NAPL phase appears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", xnn: "
                          << xnn << std::endl;
                nonwettingFlag = 1;
            }

            Scalar Smin = -1.e-6;
            if (staticDat_[globalIdx].wasSwitched)
                Smin = -0.01;

            if (volVars.saturation(gPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // gas phase disappears
                std::cout << "Gas phase disappears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", sg: "
                          << volVars.saturation(gPhaseIdx) << std::endl;
                gasFlag = 1;
            }

            Smin = 0.0;
            if (staticDat_[globalIdx].wasSwitched)
                Smin = -0.01;

            if (volVars.saturation(wPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // gas phase disappears
                std::cout << "Water phase disappears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", sw: "
                          << volVars.saturation(wPhaseIdx) << std::endl;
                wettingFlag = 1;
            }

            if ((gasFlag == 0) && (nonwettingFlag == 1) && (wettingFlag == 1))
            {
                newPhasePresence = gnPhaseOnly;
                globalSol[globalIdx][switch1Idx] = 0.0001;
                globalSol[globalIdx][switch2Idx] = volVars.fluidState().moleFraction(nPhaseIdx, wCompIdx);
;
            }
            else if ((gasFlag == 0) && (nonwettingFlag == 1) && (wettingFlag == 0))
            {
                newPhasePresence = threePhases;
                globalSol[globalIdx][switch2Idx] = 0.0001;
            }
            else if ((gasFlag == 1) && (nonwettingFlag == 0) && (wettingFlag == 0))
            {
                newPhasePresence = wPhaseOnly;
                globalSol[globalIdx][pressureIdx] = volVars.fluidState().pressure(wPhaseIdx);
                globalSol[globalIdx][switch1Idx] = volVars.fluidState().temperature();
                globalSol[globalIdx][switch2Idx] = volVars.fluidState().moleFraction(wPhaseIdx, nCompIdx);
            }
            else if ((gasFlag == 0) && (nonwettingFlag == 0) && (wettingFlag == 1))
            {
                newPhasePresence = gPhaseOnly;
                globalSol[globalIdx][switch1Idx]
                    = volVars.fluidState().temperature();
                globalSol[globalIdx][switch2Idx]
                    = volVars.fluidState().moleFraction(gPhaseIdx, nCompIdx);
            }
        }

        staticDat_[globalIdx].phasePresence = newPhasePresence;
        staticDat_[globalIdx].wasSwitched = wouldSwitch;
        return phasePresence != newPhasePresence;
    }

    //  perform variable switch at a vertex; Returns true if a
    //  variable switch was performed.
    // Note that this one is the simplified version with only two phase states
    bool primaryVarSwitchSimple_(SolutionVector &globalSol,
                           const VolumeVariables &volVars,
                           int globalIdx,
                           const GlobalPosition &globalPos)
    {
        // evaluate primary variable switch
        bool wouldSwitch = false;
        int phasePresence = staticDat_[globalIdx].phasePresence;
        int newPhasePresence = phasePresence;

        // check if a primary var switch is necessary
        if (phasePresence == threePhases)
        {
            Scalar Smin = 0;
            if (staticDat_[globalIdx].wasSwitched)
                Smin = -0.01;

            if (volVars.saturation(gPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // gas phase disappears
                std::cout << "Gas phase disappears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", sg: "
                          << volVars.saturation(gPhaseIdx) << std::endl;
                newPhasePresence = wnPhaseOnly;

                globalSol[globalIdx][pressureIdx] = volVars.fluidState().pressure(wPhaseIdx);
                globalSol[globalIdx][switch1Idx] = volVars.fluidState().temperature();
            }
        }
        else if (phasePresence == wnPhaseOnly)
        {
            bool gasFlag = 0;

            // calculate fractions of the partial pressures in the
            // hypothetical gas phase
            Scalar xwg = volVars.fluidState().moleFraction(gPhaseIdx, wCompIdx);
            Scalar xng = volVars.fluidState().moleFraction(gPhaseIdx, nCompIdx);

            Scalar xgMax = 1.0;
            if (xwg + xng > xgMax)
                wouldSwitch = true;
            if (staticDat_[globalIdx].wasSwitched)
                xgMax *= 1.02;

            // if the sum of the mole fractions would be larger than
            // 100%, gas phase appears
            if (xwg + xng > xgMax)
            {
                // gas phase appears
                std::cout << "gas phase appears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", xwg + xng: "
                          << xwg + xng << std::endl;
                gasFlag = 1;
            }

            if (gasFlag == 1)
            {
                newPhasePresence = threePhases;
                globalSol[globalIdx][pressureIdx] = volVars.pressure(gPhaseIdx);
                globalSol[globalIdx][switch1Idx] = volVars.saturation(wPhaseIdx);
            }
        }

        staticDat_[globalIdx].phasePresence = newPhasePresence;
        staticDat_[globalIdx].wasSwitched = wouldSwitch;
        return phasePresence != newPhasePresence;
    }

    // parameters given in constructor
    std::vector<StaticVars> staticDat_;
    bool switchFlag_;
};

}

#include "propertydefaults.hh"

#endif
