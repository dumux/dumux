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
 * \brief Adaption of the fully implicit scheme to the
 *        two-phase two-component fully implicit model.
 */
#ifndef DUMUX_2P2C_MODEL_HH
#define DUMUX_2P2C_MODEL_HH

#include "properties.hh"
#include "indices.hh"
#include <dumux/porousmediumflow/implicit/velocityoutput.hh>

namespace Dumux
{
/*!
 * \ingroup TwoPTwoCModel
 * \brief Adaption of the fully implicit scheme to the
 *        two-phase two-component fully implicit model.
 *
 * This model implements two-phase two-component flow of two compressible and
 * partially miscible fluids \f$\alpha \in \{ w, n \}\f$ composed of the two components
 * \f$\kappa \in \{ w, a \}\f$. The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum:
 * \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one gets one transport equation for each component
 * \f{eqnarray*}
 && \phi \frac{\partial (\sum_\alpha \varrho_\alpha \frac{M^\kappa}{M_\alpha} x_\alpha^\kappa S_\alpha )}
 {\partial t}
 - \sum_\alpha  \text{div} \left\{ \varrho_\alpha \frac{M^\kappa}{M_\alpha} x_\alpha^\kappa
 \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 (\textbf{grad}\, p_\alpha - \varrho_{\alpha}  \mbox{\bf g}) \right\}
 \nonumber \\ \nonumber \\
 &-& \sum_\alpha \text{div} \left\{ D_{\alpha,\text{pm}}^\kappa \varrho_{\alpha} \frac{M^\kappa}{M_\alpha}
 \textbf{grad} x^\kappa_{\alpha} \right\}
 - \sum_\alpha q_\alpha^\kappa = 0 \qquad \kappa \in \{w, a\} \, ,
 \alpha \in \{w, g\}
 \f}
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$ and \f$x^\kappa_w + x^\kappa_n = 1\f$, the number of
 * unknowns can be reduced to two.
 * The used primary variables are, like in the two-phase model, either \f$p_w\f$ and \f$S_n\f$
 * or \f$p_n\f$ and \f$S_w\f$. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * TwoPTwoCIndices::pWsN or TwoPTwoCIndices::pNsW. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 * Moreover, the second primary variable depends on the phase state, since a
 * primary variable switch is included. The phase state is stored for all nodes
 * of the system.
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. useMoles is set to true by default.
 * Following cases can be distinguished:
 * <ul>
 *  <li> Both phases are present: The saturation is used (either \f$S_n\f$ or \f$S_w\f$, dependent on the chosen <tt>Formulation</tt>),
 *      as long as \f$ 0 < S_\alpha < 1\f$</li>.
 *  <li> Only wetting phase is present: The mole fraction of, e.g., air in the wetting phase \f$x^a_w\f$ is used,
 *      as long as the maximum mole fraction is not exceeded \f$(x^a_w<x^a_{w,max})\f$</li>
 *  <li> Only non-wetting phase is present: The mole fraction of, e.g., water in the non-wetting phase, \f$x^w_n\f$, is used,
 *      as long as the maximum mole fraction is not exceeded \f$(x^w_n<x^w_{n,max})\f$</li>
 * </ul>
 */

template<class TypeTag>
class TwoPTwoCModel: public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseModel) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    enum {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        switchIdx = Indices::switchIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx,

        wPhaseOnly = Indices::wPhaseOnly,
        nPhaseOnly = Indices::nPhaseOnly,
        bothPhases = Indices::bothPhases,

        pwsn = TwoPTwoCFormulation::pwsn,
        pnsw = TwoPTwoCFormulation::pnsw,
        formulation = GET_PROP_VALUE(TypeTag, Formulation)
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
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

        unsigned numDofs = this->numDofs();

        staticDat_.resize(numDofs);

        setSwitched_(false);

        // check, if velocity output can be used (works only for cubes so far)
        for (const auto& element : Dune::elements(this->gridView_()))
        {
            if (!isBox) // i.e. cell-centered discretization
            {
                int eIdxGlobal = this->dofMapper().index(element);
                const GlobalPosition &globalPos = element.geometry().center();

                // initialize phase presence
                staticDat_[eIdxGlobal].phasePresence
                    = this->problem_().initialPhasePresence(*(this->gridView_().template begin<dim>()),
                                                            eIdxGlobal, globalPos);
                staticDat_[eIdxGlobal].wasSwitched = false;

                staticDat_[eIdxGlobal].oldPhasePresence
                    = staticDat_[eIdxGlobal].phasePresence;
            }
        }

        if (isBox) // i.e. vertex-centered discretization
        {
            for (const auto& vertex : Dune::vertices(this->gridView_()))
            {
                int vIdxGlobal = this->dofMapper().index(vertex);
                const GlobalPosition &globalPos = vertex.geometry().corner(0);

                // initialize phase presence
                staticDat_[vIdxGlobal].phasePresence
                    = this->problem_().initialPhasePresence(vertex, vIdxGlobal,
                                                            globalPos);
                staticDat_[vIdxGlobal].wasSwitched = false;

                staticDat_[vIdxGlobal].oldPhasePresence
                    = staticDat_[vIdxGlobal].phasePresence;
            }
        }
    }

    /*!
     * \brief Compute the total storage of all conservation quantities in one phase
     *
     * \param storage Contains the storage of each component in one phase
     * \param phaseIdx The phase index
     */
    void globalPhaseStorage(PrimaryVariables &storage, const int phaseIdx)
    {
        storage = 0;

        for (const auto& element : Dune::elements(this->gridView_())) {
            if(element.partitionType() == Dune::InteriorEntity)
            {


                this->localResidual().evalPhaseStorage(element, phaseIdx);

                for (unsigned int i = 0; i < this->localResidual().storageTerm().size(); ++i)
                    storage += this->localResidual().storageTerm()[i];
            }
        }
        if (this->gridView_().comm().size() > 1)
            storage = this->gridView_().comm().sum(storage);
    }

    /*!
     * \brief Called by the update() method if applying the Newton
     *        method was unsuccessful.
     */
    void updateFailed()
    {
        ParentType::updateFailed();

        setSwitched_(false);
        resetPhasePresence_();
    }

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
     * \brief Returns true if the primary variables were switched for
     *        at least one vertex after the last timestep.
     */
    bool switched() const
    {
        return switchFlag_;
    }

    /*!
     * \brief Returns the phase presence of the current or the old solution of a degree of freedom.
     *
     * \param dofIdxGlobal The global index of the degree of freedom
     * \param oldSol Based on oldSol current or previous time step is used
     */
    int phasePresence(int dofIdxGlobal, bool oldSol) const
    {
        return oldSol ? staticDat_[dofIdxGlobal].oldPhasePresence
            : staticDat_[dofIdxGlobal].phasePresence;
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
        typedef Dune::BlockVector<Dune::FieldVector<double, dimWorld> > VectorField;

        // get the number of degrees of freedom
        unsigned numDofs = this->numDofs();

        // create the required scalar fields
        ScalarField *sN    = writer.allocateManagedBuffer(numDofs);
        ScalarField *sW    = writer.allocateManagedBuffer(numDofs);
        ScalarField *pn    = writer.allocateManagedBuffer(numDofs);
        ScalarField *pw    = writer.allocateManagedBuffer(numDofs);
        ScalarField *pc    = writer.allocateManagedBuffer(numDofs);
        ScalarField *rhoW  = writer.allocateManagedBuffer(numDofs);
        ScalarField *rhoN  = writer.allocateManagedBuffer(numDofs);
        ScalarField *mobW  = writer.allocateManagedBuffer(numDofs);
        ScalarField *mobN = writer.allocateManagedBuffer(numDofs);
        ScalarField *phasePresence = writer.allocateManagedBuffer(numDofs);
        ScalarField *massFrac[numPhases][numComponents];
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                massFrac[phaseIdx][compIdx] = writer.allocateManagedBuffer(numDofs);
        ScalarField *moleFrac[numPhases][numComponents];
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                    for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                        moleFrac[phaseIdx][compIdx] = writer.allocateManagedBuffer(numDofs);
        ScalarField *temperature = writer.allocateManagedBuffer(numDofs);
        ScalarField *poro = writer.allocateManagedBuffer(numDofs);
        VectorField *velocityN = writer.template allocateManagedBuffer<double, dimWorld>(numDofs);
        VectorField *velocityW = writer.template allocateManagedBuffer<double, dimWorld>(numDofs);
        ImplicitVelocityOutput<TypeTag> velocityOutput(this->problem_());

        if (velocityOutput.enableOutput()) // check if velocity output is demanded
        {
            // initialize velocity fields
            for (unsigned int i = 0; i < numDofs; ++i)
            {
                (*velocityN)[i] = Scalar(0);
                (*velocityW)[i] = Scalar(0);
            }
        }

        unsigned numElements = this->gridView_().size(0);
        ScalarField *rank = writer.allocateManagedBuffer(numElements);

        for (const auto& element : Dune::elements(this->gridView_()))
        {
            if(element.partitionType() == Dune::InteriorEntity)
            {
                int eIdx = this->elementMapper().index(element);
                (*rank)[eIdx] = this->gridView_().comm().rank();

                FVElementGeometry fvGeometry;
                fvGeometry.update(this->gridView_(), element);

                ElementVolumeVariables elemVolVars;
                elemVolVars.update(this->problem_(),
                                   element,
                                   fvGeometry,
                                   false /* oldSol? */);

                for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                {
                    int dofIdxGlobal = this->dofMapper().subIndex(element, scvIdx, dofCodim);

                    (*sN)[dofIdxGlobal]    = elemVolVars[scvIdx].saturation(nPhaseIdx);
                    (*sW)[dofIdxGlobal]    = elemVolVars[scvIdx].saturation(wPhaseIdx);
                    (*pn)[dofIdxGlobal]    = elemVolVars[scvIdx].pressure(nPhaseIdx);
                    (*pw)[dofIdxGlobal]    = elemVolVars[scvIdx].pressure(wPhaseIdx);
                    (*pc)[dofIdxGlobal]    = elemVolVars[scvIdx].capillaryPressure();
                    (*rhoW)[dofIdxGlobal]  = elemVolVars[scvIdx].density(wPhaseIdx);
                    (*rhoN)[dofIdxGlobal]  = elemVolVars[scvIdx].density(nPhaseIdx);
                    (*mobW)[dofIdxGlobal]  = elemVolVars[scvIdx].mobility(wPhaseIdx);
                    (*mobN)[dofIdxGlobal]  = elemVolVars[scvIdx].mobility(nPhaseIdx);
                    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                        {
                            (*massFrac[phaseIdx][compIdx])[dofIdxGlobal]
                                = elemVolVars[scvIdx].massFraction(phaseIdx, compIdx);

                            Valgrind::CheckDefined((*massFrac[phaseIdx][compIdx])[dofIdxGlobal][0]);
                        }
                    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                        {
                            (*moleFrac[phaseIdx][compIdx])[dofIdxGlobal]
                                = elemVolVars[scvIdx].moleFraction(phaseIdx, compIdx);

                            Valgrind::CheckDefined((*moleFrac[phaseIdx][compIdx])[dofIdxGlobal][0]);
                        }
                    (*poro)[dofIdxGlobal]  = elemVolVars[scvIdx].porosity();
                    (*temperature)[dofIdxGlobal] = elemVolVars[scvIdx].temperature();
                    (*phasePresence)[dofIdxGlobal]
                        = staticDat_[dofIdxGlobal].phasePresence;
                }

                // velocity output
                velocityOutput.calculateVelocity(*velocityW, elemVolVars, fvGeometry, element, wPhaseIdx);
                velocityOutput.calculateVelocity(*velocityN, elemVolVars, fvGeometry, element, nPhaseIdx);
            }

        } // loop over elements

        writer.attachDofData(*sN,     "Sn", isBox);
        writer.attachDofData(*sW,     "Sw", isBox);
        writer.attachDofData(*pn,     "pn", isBox);
        writer.attachDofData(*pw,     "pw", isBox);
        writer.attachDofData(*pc,     "pc", isBox);
        writer.attachDofData(*rhoW,   "rhoW", isBox);
        writer.attachDofData(*rhoN,   "rhoN", isBox);
        writer.attachDofData(*mobW,   "mobW", isBox);
        writer.attachDofData(*mobN,   "mobN", isBox);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                std::ostringstream oss;
                oss << "X_" << FluidSystem::phaseName(phaseIdx) << "^" << FluidSystem::componentName(compIdx);
                writer.attachDofData(*massFrac[phaseIdx][compIdx], oss.str(), isBox);
            }
        }
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                std::ostringstream oss;
                oss << "x_" << FluidSystem::phaseName(phaseIdx) << "^" << FluidSystem::componentName(compIdx);
                writer.attachDofData(*moleFrac[phaseIdx][compIdx], oss.str(), isBox);
            }
        }
        writer.attachDofData(*poro, "porosity", isBox);
        writer.attachDofData(*temperature,    "temperature", isBox);
        writer.attachDofData(*phasePresence,  "phase presence", isBox);

        if (velocityOutput.enableOutput()) // check if velocity output is demanded
        {
            writer.attachDofData(*velocityW,  "velocityW", isBox, dim);
            writer.attachDofData(*velocityN,  "velocityN", isBox, dim);
        }

        writer.attachCellData(*rank, "process rank");
    }

    /*!
     * \brief Write the current solution to a restart file.
     *
     * \param outStream The output stream of one vertex for the restart file
     * \param entity The entity, either a vertex or an element
     */
    template<class Entity>
    void serializeEntity(std::ostream &outStream, const Entity &entity)
    {
        // write primary variables
        ParentType::serializeEntity(outStream, entity);

        int dofIdxGlobal = this->dofMapper().index(entity);

        if (!outStream.good())
            DUNE_THROW(Dune::IOError, "Could not serialize entity " << dofIdxGlobal);

        outStream << staticDat_[dofIdxGlobal].phasePresence << " ";
    }

    /*!
     * \brief Reads the current solution from a restart file.
     *
     * \param inStream The input stream of one vertex from the restart file
     * \param entity The entity, either a vertex or an element
     */
    template<class Entity>
    void deserializeEntity(std::istream &inStream, const Entity &entity)
    {
        // read primary variables
        ParentType::deserializeEntity(inStream, entity);

        // read phase presence
        int dofIdxGlobal = this->dofMapper().index(entity);

        if (!inStream.good())
            DUNE_THROW(Dune::IOError,
                       "Could not deserialize entity " << dofIdxGlobal);

        inStream >> staticDat_[dofIdxGlobal].phasePresence;
        staticDat_[dofIdxGlobal].oldPhasePresence
            = staticDat_[dofIdxGlobal].phasePresence;

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
        int succeeded;
        try {
            for (unsigned i = 0; i < staticDat_.size(); ++i)
                staticDat_[i].visited = false;

            FVElementGeometry fvGeometry;
            static VolumeVariables volVars;
            for (const auto& element : Dune::elements(this->gridView_()))
            {
                fvGeometry.update(this->gridView_(), element);
                for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                {
                    int dofIdxGlobal = this->dofMapper().subIndex(element, scvIdx, dofCodim);

                    if (staticDat_[dofIdxGlobal].visited)
                        continue;

                    staticDat_[dofIdxGlobal].visited = true;
                    volVars.update(curGlobalSol[dofIdxGlobal],
                            this->problem_(),
                            element,
                            fvGeometry,
                            scvIdx,
                            false);
                    const GlobalPosition &globalPos = fvGeometry.subContVol[scvIdx].global;
                    if (primaryVarSwitch_(curGlobalSol,
                            volVars,
                            dofIdxGlobal,
                            globalPos))
                    {
                        this->jacobianAssembler().markDofRed(dofIdxGlobal);
                        wasSwitched = true;
                    }
                }
            }
            succeeded = 1;
        }
        catch (Dumux::NumericalProblem &e)
        {
            std::cout << "\n"
                      << "Rank " << this->problem_().gridView().comm().rank()
                      << " caught an exception while updating the static data." << e.what()
                      << "\n";
            succeeded = 0;
        }
        //make sure that all processes succeeded. If not throw a NumericalProblem to decrease the time step size.
        if (this->gridView_().comm().size() > 1)
            succeeded = this->gridView_().comm().min(succeeded);

        if (!succeeded) {
                DUNE_THROW(NumericalProblem,
                        "A process did not succeed in updating the static data.");
            return;
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
     * \brief Resets the current phase presence of all vertices to the old one.
     *
     * This is done after an update failed.
     */
    void resetPhasePresence_()
    {
        for (unsigned int idx = 0; idx < staticDat_.size(); ++idx)
        {
            staticDat_[idx].phasePresence
                = staticDat_[idx].oldPhasePresence;
            staticDat_[idx].wasSwitched = false;
        }
    }

    /*!
     * \brief Sets the phase presence of all vertices state to the current one.
     */
    void updateOldPhasePresence_()
    {
        for (unsigned int idx = 0; idx < staticDat_.size(); ++idx)
        {
            staticDat_[idx].oldPhasePresence
                = staticDat_[idx].phasePresence;
            staticDat_[idx].wasSwitched = false;
        }
    }

    /*!
     * \brief Sets whether there was a primary variable switch after
     *        the last timestep.
     */
    void setSwitched_(bool yesno)
    {
        switchFlag_ = yesno;
    }

    /*!
     * \brief Performs variable switch at a vertex, returns true if a
     *        variable switch was performed.
     */
    bool primaryVarSwitch_(SolutionVector &globalSol,
                           const VolumeVariables &volVars,
                           int dofIdxGlobal,
                           const GlobalPosition &globalPos)
    {
        // evaluate primary variable switch
        bool wouldSwitch = false;
        int phasePresence = staticDat_[dofIdxGlobal].phasePresence;
        int newPhasePresence = phasePresence;

        // check if a primary var switch is necessary
        if (phasePresence == nPhaseOnly)
        {
            // calculate mole fraction in the hypothetic wetting phase
            Scalar xww = volVars.moleFraction(wPhaseIdx, wCompIdx);
            Scalar xwn = volVars.moleFraction(wPhaseIdx, nCompIdx);

            Scalar xwMax = 1.0;
            if (xww + xwn > xwMax)
                wouldSwitch = true;
            if (staticDat_[dofIdxGlobal].wasSwitched)
                xwMax *= 1.02;

            // if the sum of the mole fractions is larger than
            // 100%, wetting phase appears
            if (xww + xwn > xwMax)
            {
                // wetting phase appears
                std::cout << "wetting phase appears at vertex " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", xww + xwn: "
                          << xww + xwn << std::endl;
                newPhasePresence = bothPhases;
                if (formulation == pnsw)
                    globalSol[dofIdxGlobal][switchIdx] = 0.0;
                else if (formulation == pwsn)
                    globalSol[dofIdxGlobal][switchIdx] = 1.0;
            }
        }
        else if (phasePresence == wPhaseOnly)
        {
            // calculate fractions of the partial pressures in the
            // hypothetic nonwetting phase
            Scalar xnw = volVars.moleFraction(nPhaseIdx, wCompIdx);
            Scalar xnn = volVars.moleFraction(nPhaseIdx, nCompIdx);

            Scalar xgMax = 1.0;
            if (xnw + xnn > xgMax)
                wouldSwitch = true;
            if (staticDat_[dofIdxGlobal].wasSwitched)
                xgMax *= 1.02;

            // if the sum of the mole fractions is larger than
            // 100%, nonwetting phase appears
            if (xnw + xnn > xgMax)
            {
                // nonwetting phase appears
                std::cout << "nonwetting phase appears at vertex " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", xnw + xnn: "
                          << xnw + xnn << std::endl;
                newPhasePresence = bothPhases;
                if (formulation == pnsw)
                    globalSol[dofIdxGlobal][switchIdx] = 0.999;
                else if (formulation == pwsn)
                    globalSol[dofIdxGlobal][switchIdx] = 0.001;
            }
        }
        else if (phasePresence == bothPhases)
        {
            Scalar Smin = 0.0;
            if (staticDat_[dofIdxGlobal].wasSwitched)
                Smin = -0.01;

            if (volVars.saturation(nPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // nonwetting phase disappears
                std::cout << "Nonwetting phase disappears at vertex " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", sn: "
                          << volVars.saturation(nPhaseIdx) << std::endl;
                newPhasePresence = wPhaseOnly;

                if(useMoles) // mole-fraction formulation
                {
                    globalSol[dofIdxGlobal][switchIdx]
                        = volVars.moleFraction(wPhaseIdx, nCompIdx);
                }
                else // mass-fraction formulation
                {
                    globalSol[dofIdxGlobal][switchIdx]
                        = volVars.massFraction(wPhaseIdx, nCompIdx);
                }
            }
            else if (volVars.saturation(wPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // wetting phase disappears
                std::cout << "Wetting phase disappears at vertex " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", sw: "
                          << volVars.saturation(wPhaseIdx) << std::endl;
                newPhasePresence = nPhaseOnly;

                if(useMoles) // mole-fraction formulation
                {
                    globalSol[dofIdxGlobal][switchIdx]
                        = volVars.moleFraction(nPhaseIdx, wCompIdx);
                }
                else // mass-fraction formulation
                {
                    globalSol[dofIdxGlobal][switchIdx]
                        = volVars.massFraction(nPhaseIdx, wCompIdx);
                }
            }
        }

        staticDat_[dofIdxGlobal].phasePresence = newPhasePresence;
        staticDat_[dofIdxGlobal].wasSwitched = wouldSwitch;
        return phasePresence != newPhasePresence;
    }

protected:
    // parameters given in constructor
    std::vector<StaticVars> staticDat_;
    bool switchFlag_;
};

}

#include "propertydefaults.hh"

#endif
