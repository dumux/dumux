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
 * \brief Adaption of the fully implicit scheme to the three-phase three-component
 *        flow model.
 *
 * The model is designed for simulating three fluid phases with water, gas, and
 * a liquid contaminant (NAPL - non-aqueous phase liquid)
 */
#ifndef DUMUX_3P3C_MODEL_HH
#define DUMUX_3P3C_MODEL_HH

#include <dumux/porousmediumflow/implicit/velocityoutput.hh>
#include "properties.hh"
#include "primaryvariableswitch.hh"

namespace Dumux
{
/*!
 * \ingroup ThreePThreeCModel
 * \brief Adaption of the fully implicit scheme to the three-phase three-component
 *        flow model.
 *
 * This model implements three-phase three-component flow of three fluid phases
 * \f$\alpha \in \{ water, gas, NAPL \}\f$ each composed of up to three components
 * \f$\kappa \in \{ water, air, contaminant \}\f$. The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum:
 * \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K}
 \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one transport equation for each component is obtained as
 * \f{eqnarray*}
 && \phi \frac{\partial (\sum_\alpha \varrho_{\alpha,mol} x_\alpha^\kappa
 S_\alpha )}{\partial t}
 - \sum\limits_\alpha \text{div} \left\{ \frac{k_{r\alpha}}{\mu_\alpha}
 \varrho_{\alpha,mol} x_\alpha^\kappa \mathbf{K}
 (\textbf{grad}\, p_\alpha - \varrho_{\alpha,mass} \mbox{\bf g}) \right\}
 \nonumber \\
 \nonumber \\
 && - \sum\limits_\alpha \text{div} \left\{ D_\text{pm}^\kappa \varrho_{\alpha,mol}
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
 * \f$x^w_\alpha + x^a_\alpha + x^c_\alpha = 1\f$ for the mole fractions.
 * Furthermore, the phase pressures are related to each other via
 * capillary pressures between the fluid phases, which are functions of
 * the saturation, e.g. according to the approach of Parker et al.
 *
 * The used primary variables are dependent on the locally present fluid phases.
 * An adaptive primary variable switch is included. The phase state is stored for all nodes
 * of the system. The following cases can be distinguished:
 * <ul>
 *  <li> All three phases are present: Primary variables are two saturations \f$(S_w\f$ and \f$S_n)\f$,
 *       and a pressure, in this case \f$p_g\f$. </li>
 *  <li> Only the water phase is present: Primary variables are now the mole fractions of air and
 *       contaminant in the water phase \f$(x_w^a\f$ and \f$x_w^c)\f$, as well as the gas pressure, which is,
 *       of course, in a case where only the water phase is present, just the same as the water pressure. </li>
 *  <li> Gas and NAPL phases are present: Primary variables \f$(S_n\f$, \f$x_g^w\f$, \f$p_g)\f$. </li>
 *  <li> Water and NAPL phases are present: Primary variables \f$(S_n\f$, \f$x_w^a\f$, \f$p_g)\f$. </li>
 *  <li> Only gas phase is present: Primary variables \f$(x_g^w\f$, \f$x_g^c\f$, \f$p_g)\f$. </li>
 *  <li> Water and gas phases are present: Primary variables \f$(S_w\f$, \f$x_w^g\f$, \f$p_g)\f$. </li>
 * </ul>
 */
template<class TypeTag>
class ThreePThreeCModel: public GET_PROP_TYPE(TypeTag, BaseModel)
{
    friend typename GET_PROP_TYPE(TypeTag, BaseModel);
    typedef typename GET_PROP_TYPE(TypeTag, BaseModel) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,

        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx,
        gCompIdx = Indices::gCompIdx,

        threePhases = Indices::threePhases,
        wPhaseOnly  = Indices::wPhaseOnly,
        gnPhaseOnly = Indices::gnPhaseOnly,
        wnPhaseOnly = Indices::wnPhaseOnly,
        gPhaseOnly  = Indices::gPhaseOnly,
        wgPhaseOnly = Indices::wgPhaseOnly

    };

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

public:
    /*!
     * \brief One Newton iteration was finished.
     * \param uCurrent The solution after the current Newton iteration
     */
    void newtonEndStep()
    {
        // \todo resize volvars vector if grid was adapted

        // update the variable switch
        switchFlag_ = priVarSwitch_().update(this->problem_(), this->curSol(),
                                             this->curVolVars_());

        // update the secondary variables
        // \note we only updated if phase presence changed as the volume variables
        //       are already updated once by the switch
        if (switchFlag_)
        {
            for (const auto& element : elements(this->problem_().gridView()))
            {
                // make sure FVElementGeometry & vol vars are bound to the element
                this->fvGeometries_().bindElement(element);
                this->curVolVars_().bindElement(element);

                const auto& fvGeometry = this->fvGeometries(element);
                for (const auto& scv : fvGeometry.scvs())
                {
                    auto dofIdxGlobal = scv.dofIndex();
                    if (priVarSwitch_().wasSwitched(dofIdxGlobal))
                    {
                        this->curVolVars_(dofIdxGlobal).update(this->curSol()[dofIdxGlobal],
                                                               this->problem_(),
                                                               element,
                                                               scv);
                    }
                }

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
     *        method was unsuccessful.
     */
    void updateFailed()
    {
        ParentType::updateFailed();

        switchFlag_ = false;
        priVarSwitch_().resetPhasePresence();
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
        priVarSwitch_().updateOldPhasePresence();
        switchFlag_ = false;
    }

    /*!
     * \brief Returns true if the primary variables were switched for
     *        at least one dof after the last timestep.
     */
    bool switched() const
    {
        return switchFlag_;
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
        // typedef Dune::BlockVector<Dune::FieldVector<double, dimWorld> > VectorField;

        // get the number of degrees of freedom
        unsigned numDofs = this->numDofs();

        // create the required scalar fields
        ScalarField *saturation[numPhases];
        ScalarField *pressure[numPhases];
        ScalarField *density[numPhases];

        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
        {
            saturation[phaseIdx] = writer.allocateManagedBuffer(numDofs);
            pressure[phaseIdx] = writer.allocateManagedBuffer(numDofs);
            density[phaseIdx] = writer.allocateManagedBuffer(numDofs);
        }

        ScalarField *phasePresence = writer.allocateManagedBuffer (numDofs);
        ScalarField *moleFraction[numPhases][numComponents];
        for (int i = 0; i < numPhases; ++i)
            for (int j = 0; j < numComponents; ++j)
                moleFraction[i][j] = writer.allocateManagedBuffer (numDofs);
        ScalarField *temperature = writer.allocateManagedBuffer (numDofs);
        ScalarField *poro = writer.allocateManagedBuffer(numDofs);
        // VectorField *velocityN = writer.template allocateManagedBuffer<double, dimWorld>(numDofs);
        // VectorField *velocityW = writer.template allocateManagedBuffer<double, dimWorld>(numDofs);
        // VectorField *velocityG = writer.template allocateManagedBuffer<double, dimWorld>(numDofs);
        // ImplicitVelocityOutput<TypeTag> velocityOutput(this->problem_());

        // if (velocityOutput.enableOutput()) // check if velocity output is demanded
        // {
        //     // initialize velocity fields
        //     for (unsigned int i = 0; i < numDofs; ++i)
        //     {
        //         (*velocityN)[i] = Scalar(0);
        //         (*velocityW)[i] = Scalar(0);
        //         (*velocityG)[i] = Scalar(0);
        //     }
        // }

        unsigned numElements = this->gridView_().size(0);
        ScalarField *rank = writer.allocateManagedBuffer (numElements);

        for (const auto& element : elements(this->gridView_(), Dune::Partitions::interior))
        {
            int eIdx = this->problem_().elementMapper().index(element);
            (*rank)[eIdx] = this->gridView_().comm().rank();

            // make sure FVElementGeometry & vol vars are bound to the element
            this->fvGeometries_().bindElement(element);
            this->curVolVars_().bindElement(element);

            const auto& fvGeometry = this->fvGeometries(element);
            for (const auto& scv : fvGeometry.scvs())
            {
                const auto& volVars = this->curVolVars(scv);
                int dofIdxGlobal = scv.dofIndex();

                for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
                {
                    (*saturation[phaseIdx])[dofIdxGlobal] = volVars.saturation(phaseIdx);
                    (*pressure[phaseIdx])[dofIdxGlobal] = volVars.pressure(phaseIdx);
                    (*density[phaseIdx])[dofIdxGlobal] = volVars.density(phaseIdx);
                }

                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                {
                    for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    {
                        (*moleFraction[phaseIdx][compIdx])[dofIdxGlobal] =
                            volVars.moleFraction(phaseIdx, compIdx);

                    }
                }

                (*poro)[dofIdxGlobal] = volVars.porosity();
                (*temperature)[dofIdxGlobal] = volVars.temperature();
                (*phasePresence)[dofIdxGlobal] = priVarSwitch().phasePresence(dofIdxGlobal);

            }

            // // velocity output
            // velocityOutput.calculateVelocity(*velocityW, elemVolVars, fvGeometry, element, wPhaseIdx);
            // velocityOutput.calculateVelocity(*velocityN, elemVolVars, fvGeometry, element, nPhaseIdx);
            // velocityOutput.calculateVelocity(*velocityN, elemVolVars, fvGeometry, element, gPhaseIdx);
        }

        writer.attachDofData(*saturation[wPhaseIdx], "Sw", isBox);
        writer.attachDofData(*saturation[nPhaseIdx], "Sn", isBox);
        writer.attachDofData(*saturation[gPhaseIdx], "Sg", isBox);
        writer.attachDofData(*pressure[wPhaseIdx], "pw", isBox);
        writer.attachDofData(*pressure[nPhaseIdx], "pn", isBox);
        writer.attachDofData(*pressure[gPhaseIdx], "pg", isBox);
        writer.attachDofData(*density[wPhaseIdx], "rhow", isBox);
        writer.attachDofData(*density[nPhaseIdx], "rhon", isBox);
        writer.attachDofData(*density[gPhaseIdx], "rhog", isBox);

        for (int i = 0; i < numPhases; ++i)
        {
            for (int j = 0; j < numComponents; ++j)
            {
                std::ostringstream oss;
                oss << "x^"
                    << FluidSystem::componentName(j)
                    << "_"
                    << FluidSystem::phaseName(i);
                writer.attachDofData(*moleFraction[i][j], oss.str().c_str(), isBox);
            }
        }
        writer.attachDofData(*poro, "porosity", isBox);
        writer.attachDofData(*temperature, "temperature", isBox);
        writer.attachDofData(*phasePresence, "phase presence", isBox);

        // if (velocityOutput.enableOutput()) // check if velocity output is demanded
        // {
        //     writer.attachDofData(*velocityW,  "velocityW", isBox, dim);
        //     writer.attachDofData(*velocityN,  "velocityN", isBox, dim);
        //     writer.attachDofData(*velocityG,  "velocityG", isBox, dim);
        // }

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

        int dofIdxGlobal = this->dofMapper().index(entity);

        if (!outStream.good())
            DUNE_THROW(Dune::IOError, "Could not serialize entity " << dofIdxGlobal);

        outStream << priVarSwitch().phasePresence(dofIdxGlobal) << " ";
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
        int dofIdxGlobal = this->dofMapper().index(entity);

        if (!inStream.good())
            DUNE_THROW(Dune::IOError, "Could not deserialize entity " << dofIdxGlobal);

        int phasePresence;
        inStream >> phasePresence;

        priVarSwitch_().setPhasePresence(dofIdxGlobal, phasePresence);
        priVarSwitch_().setOldPhasePresence(dofIdxGlobal, phasePresence);
    }

    const Dumux::ThreePThreeCPrimaryVariableSwitch<TypeTag>& priVarSwitch() const
    { return switch_; }

protected:

    Dumux::ThreePThreeCPrimaryVariableSwitch<TypeTag>& priVarSwitch_()
    { return switch_; }

    /*!
     * \brief Applies the initial solution for all vertices of the grid.
     *
     * \todo the initial condition needs to be unique for
     *       each vertex. we should think about the API...
     */
    void applyInitialSolution_()
    {
        ParentType::applyInitialSolution_();

        // initialize the primary variable switch
        priVarSwitch_().init(this->problem_());
    }

    //! the class handling the primary variable switch
    Dumux::ThreePThreeCPrimaryVariableSwitch<TypeTag> switch_;
    bool switchFlag_;
};

}

#include "propertydefaults.hh"

#endif
