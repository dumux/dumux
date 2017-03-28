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
* \brief Adaption of the fully implicit model to the two-phase n-component flow model.
*/

#ifndef DUMUX_2PNC_MODEL_HH
#define DUMUX_2PNC_MODEL_HH

#include <dumux/porousmediumflow/implicit/velocityoutput.hh>

#include "properties.hh"
#include "indices.hh"
#include "localresidual.hh"

namespace Dumux
{
/*!
 * \ingroup TwoPNCModel
 * \brief Adaption of the fully implicit scheme to the
 *        two-phase n-component fully implicit model.
 *
 * This model implements two-phase n-component flow of two compressible and
 * partially miscible fluids \f$\alpha \in \{ w, n \}\f$ composed of the n components
 * \f$\kappa \in \{ w, n,\cdots \}\f$ in combination with mineral precipitation and dissolution.
 * The solid phases. The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum:
 * \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 \left(\text{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one gets one transport equation for each component
 * \f{eqnarray}
 && \frac{\partial (\sum_\alpha \varrho_\alpha X_\alpha^\kappa \phi S_\alpha )}
 {\partial t}
 - \sum_\alpha  \text{div} \left\{ \varrho_\alpha X_\alpha^\kappa
 \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 (\text{grad}\, p_\alpha - \varrho_{\alpha}  \mbox{\bf g}) \right\}
 \nonumber \\ \nonumber \\
    &-& \sum_\alpha \text{div} \left\{{\bf D_{\alpha, pm}^\kappa} \varrho_{\alpha} \text{grad}\, X^\kappa_{\alpha} \right\}
 - \sum_\alpha q_\alpha^\kappa = 0 \qquad \kappa \in \{w, a,\cdots \} \, ,
 \alpha \in \{w, g\}
 \f}
 *
 * The solid or mineral phases are assumed to consist of a single component.
 * Their mass balance consist only of a storage and a source term:
 *  \f$\frac{\partial \varrho_\lambda \phi_\lambda )} {\partial t}
 *  = q_\lambda\f$
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as
 * spatial and the implicit Euler method as time discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$ and \f$X^\kappa_w + X^\kappa_n = 1\f$, the number of
 * unknowns can be reduced to number of components.
 *
 * The used primary variables are, like in the two-phase model, either \f$p_w\f$ and \f$S_n\f$
 * or \f$p_n\f$ and \f$S_w\f$. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * TwoPTwoCIndices::pwsn or TwoPTwoCIndices::pnsw. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 *
 * Moreover, the second primary variable depends on the phase state, since a
 * primary variable switch is included. The phase state is stored for all nodes
 * of the system. The model is uses mole fractions.
 *Following cases can be distinguished:
 * <ul>
 *  <li> Both phases are present: The saturation is used (either \f$S_n\f$ or \f$S_w\f$, dependent on the chosen <tt>Formulation</tt>),
 *      as long as \f$ 0 < S_\alpha < 1\f$</li>.
 *  <li> Only wetting phase is present: The mole fraction of, e.g., air in the wetting phase \f$x^a_w\f$ is used,
 *      as long as the maximum mole fraction is not exceeded (\f$x^a_w<x^a_{w,max}\f$)</li>
 *  <li> Only non-wetting phase is present: The mole fraction of, e.g., water in the non-wetting phase, \f$x^w_n\f$, is used,
 *      as long as the maximum mole fraction is not exceeded (\f$x^w_n<x^w_{n,max}\f$)</li>
 * </ul>
 *
 * For the other components, the mole fraction \f$x^\kappa_w\f$ is the primary variable.
 */

template<class TypeTag>
class TwoPNCModel: public GET_PROP_TYPE(TypeTag, BaseModel)
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

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum {  numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum {  numMajorComponents = GET_PROP_VALUE(TypeTag, NumMajorComponents) };

    enum {
            pressureIdx = Indices::pressureIdx,
            switchIdx = Indices::switchIdx
    };
    enum {
            wPhaseIdx = Indices::wPhaseIdx,
            nPhaseIdx = Indices::nPhaseIdx
    };
    enum {
            wCompIdx = FluidSystem::wCompIdx,
            nCompIdx = FluidSystem::nCompIdx
    };
    enum {
            wPhaseOnly = Indices::wPhaseOnly,
            nPhaseOnly = Indices::nPhaseOnly,
            bothPhases = Indices::bothPhases
    };
    enum {
            pwsn = TwoPNCFormulation::pwsn,
            pnsw = TwoPNCFormulation::pnsw,
            formulation = GET_PROP_VALUE(TypeTag, Formulation)
    };

    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<0>::Entity Element;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldMatrix<CoordScalar, dimWorld, dimWorld> Tensor;

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

        auto numDofs = this->numDofs();

        staticDat_.resize(numDofs);

        setSwitched_(false);

        for (const auto& element : elements(this->gridView_()))
        {
            FVElementGeometry fvGeometry;
            fvGeometry.update(this->gridView_(), element);
            for (unsigned int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
            {
                auto dofIdxGlobal = this->dofMapper().subIndex(element, scvIdx, dofCodim);
                staticDat_[dofIdxGlobal].phasePresence = this->problem_().initialPhasePresence(element, fvGeometry, scvIdx);
                staticDat_[dofIdxGlobal].wasSwitched = false;
                staticDat_[dofIdxGlobal].oldPhasePresence = staticDat_[dofIdxGlobal].phasePresence;
            }
        }
    }

    /*!
     * \brief Compute the total storage of all conservation quantities in one phase
     *
     * \param storage Contains the storage of each component in one phase
     * \param phaseIdx The phase index
     */
    void globalPhaseStorage(PrimaryVariables &storage, int phaseIdx)
    {
        storage = 0;

        for (const auto& element : elements(this->gridView_(), Dune::Partitions::interior))
        {
            this->localResidual().evalPhaseStorage(element, phaseIdx);

            for (unsigned int i = 0; i < this->localResidual().storageTerm().size(); ++i)
                storage += this->localResidual().storageTerm()[i];
        }

        this->gridView_().comm().sum(storage);
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
     * \brief Return true if the primary variables were switched for
     *        at least one vertex after the last timestep.
     */
    bool switched() const
    {
        return switchFlag_;
    }

    /*!
     * \brief Returns the phase presence of the current or the old solution
     *
     * \param dofIdxGlobal The global DOF index
     * \param oldSol Evaluate function with solution of current or previous time step
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

        for (const auto& element : elements(this->gridView_(), Dune::Partitions::interior))
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
     * \param entity The Entity
     */
    template<class Entity>
    void serializeEntity(std::ostream &outStream, const Entity &entity)
    {
        // write primary variables
        ParentType::serializeEntity(outStream, entity);
        int dofIdxGlobal = this->dofMapper().index(entity);
        outStream << staticDat_[dofIdxGlobal].phasePresence << " ";
    }

    /*!
     * \brief Reads the current solution for a vertex from a restart
     *        file.
     *
     * \param inStream The input stream of one vertex from the restart file
     * \param entity The Entity
     */
    template<class Entity>
    void deserializeEntity(std::istream &inStream, const Entity &entity)
    {
        // read primary variables
        ParentType::deserializeEntity(inStream, entity);
        int dofIdxGlobal = this->dofMapper().index(entity);

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

        for (unsigned i = 0; i < staticDat_.size(); ++i)
            staticDat_[i].visited = false;

        for (const auto& element : elements(this->gridView_()))
        {
            FVElementGeometry fvGeometry;
            fvGeometry.update(this->gridView_(), element);
            for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
            {
                auto dofIdxGlobal = this->dofMapper().subIndex(element, scvIdx, dofCodim);

                if (staticDat_[dofIdxGlobal].visited)
                    continue;

                staticDat_[dofIdxGlobal].visited = true;
                VolumeVariables volVars;
                volVars.update(curGlobalSol[dofIdxGlobal],
                               this->problem_(),
                               element,
                               fvGeometry,
                               scvIdx,
                               false);
                const GlobalPosition &global = element.geometry().corner(scvIdx);
                if (primaryVarSwitch_(curGlobalSol, volVars, dofIdxGlobal, global))
                    wasSwitched = true;
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

    Tensor perm_(Scalar perm)
    {
        Tensor K(0.0);

        for(int i=0; i<dim; i++)
            K[i][i] = perm;

       return K;
    }

    Tensor perm_(Tensor perm)
    {
       return perm;
    }

    /*!
     * \brief Reset the current phase presence of all vertices to the old one.
     *
     * This is done after an update failed.
     */
    void resetPhasePresence_()
    {
        for (int i = 0; i < this->numDofs(); ++i)
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
        for (int i = 0; i < this->numDofs(); ++i)
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

    /*!
     * \brief  perform variable switch at a vertex; Returns true if a
     *         variable switch was performed.
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

            //check if a primary variable switch is necessary
            if (phasePresence == bothPhases)
            {
                Scalar Smin = 0; //saturation threshold
                if (staticDat_[dofIdxGlobal].wasSwitched)
                    Smin = -0.01;

                //if saturation of wetting phase is smaller 0 switch
                if (volVars.saturation(wPhaseIdx) <= Smin)
                {
                    wouldSwitch = true;
                    //wetting phase has to disappear
                    std::cout << "Wetting Phase disappears at vertex " << dofIdxGlobal
                                << ", coordinated: " << globalPos << ", Sw: "
                                << volVars.saturation(wPhaseIdx) << std::endl;
                    newPhasePresence = nPhaseOnly;

                    //switch saturation to xnH20 (not depending on formulation)
                    globalSol[dofIdxGlobal][switchIdx]
                            = volVars.moleFraction(nPhaseIdx, wCompIdx /*H2O*/);

                    //switch all secondary components to mole fraction in nonwetting phase
                    for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
                        globalSol[dofIdxGlobal][compIdx] = volVars.moleFraction(nPhaseIdx,compIdx);
                }
                //if saturation of nonwetting phase is smaller than 0 switch
                else if (volVars.saturation(nPhaseIdx) <= Smin)
                {
                    wouldSwitch = true;
                    //nonwetting phase has to disappear
                    std::cout << "Nonwetting Phase disappears at vertex " << dofIdxGlobal
                                << ", coordinated: " << globalPos << ", Sn: "
                                << volVars.saturation(nPhaseIdx) << std::endl;
                    newPhasePresence = wPhaseOnly;

                    //switch saturation to xwN2, not depending on formulation
                    globalSol[dofIdxGlobal][switchIdx]
                            = volVars.moleFraction(wPhaseIdx, nCompIdx /*N2*/);
                }
            }
            else if (phasePresence == nPhaseOnly)
            {
                Scalar xwmax = 1;
                Scalar sumxw = 0;
                //Calculate sum of mole fractions in the hypothetical wetting phase
                for (int compIdx = 0; compIdx < numComponents; compIdx++)
                {
                    sumxw += volVars.moleFraction(wPhaseIdx, compIdx);
                }
                if (sumxw > xwmax)
                    wouldSwitch = true;
                if (staticDat_[dofIdxGlobal].wasSwitched)
                    xwmax *=1.02;
                //wetting phase appears if sum is larger than one
                if (sumxw/*sum of mole fractions*/ > xwmax/*1*/)
                {
                    std::cout << "Wetting Phase appears at vertex " << dofIdxGlobal
                            << ", coordinated: " << globalPos << ", sumxw: "
                            << sumxw << std::endl;
                    newPhasePresence = bothPhases;

                    //saturation of the wetting phase set to 0.0001
                    if (formulation == pnsw)
                        globalSol[dofIdxGlobal][switchIdx] = 0.0001;
                    else if (formulation == pwsn)
                        globalSol[dofIdxGlobal][switchIdx] = 0.9999;

                    //switch all secondary components back to wetting mole fraction
                    for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
                        globalSol[dofIdxGlobal][compIdx] = volVars.moleFraction(wPhaseIdx,compIdx);
                }
            }
            else if (phasePresence == wPhaseOnly)
            {
                Scalar xnmax = 1;
                Scalar sumxn = 0;
                //Calculate sum of mole fractions in the hypothetical nonwetting phase
                for (int compIdx = 0; compIdx < numComponents; compIdx++)
                {
                    sumxn += volVars.moleFraction(nPhaseIdx, compIdx);
                }
                if (sumxn > xnmax)
                    wouldSwitch = true;
                if (staticDat_[dofIdxGlobal].wasSwitched)
                    xnmax *=1.02;
                //nonwetting phase appears if sum is larger than one
                if (sumxn > xnmax)
                {
                    std::cout << "Nonwetting Phase appears at vertex " << dofIdxGlobal
                            << ", coordinated: " << globalPos << ", sumxn: "
                            << sumxn << std::endl;
                    newPhasePresence = bothPhases;
                    //saturation of the wetting phase set to 0.9999
                    if (formulation == pnsw)
                        globalSol[dofIdxGlobal][switchIdx] = 0.9999;
                    else if (formulation == pwsn)
                        globalSol[dofIdxGlobal][switchIdx] = 0.0001;
                }
            }

            staticDat_[dofIdxGlobal].phasePresence = newPhasePresence;
            staticDat_[dofIdxGlobal].wasSwitched = wouldSwitch;
            return phasePresence != newPhasePresence;
        }

    // parameters given in constructor
    std::vector<StaticVars> staticDat_;
    bool switchFlag_;
};

}

#include "propertydefaults.hh"

#endif
