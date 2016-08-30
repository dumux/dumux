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
 * \brief Adaption of the fully implicit scheme to the two-phase one-component flow model.
 *
 * Important note: The 2p1c model requires the use of the non-isothermal extension found in dumux/implicit/nonisothermal
 *
 * The model is designed for simulating two fluid phases with water as the only component.
 * This model is  particularly suitable for the simulation of steam injection in saturated conditions.
 */
#ifndef DUMUX_2P1C_MODEL_HH
#define DUMUX_2P1C_MODEL_HH

#include <dumux/porousmediumflow/implicit/velocityoutput.hh>
#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup TwoPOneCNIModel
 * \brief Adaption of the fully implicit scheme to the two-phase one-component flow model.
 *
 * This model implements two-phase one-component flow of three fluid phases
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
class TwoPOneCModel: public GET_PROP_TYPE(TypeTag, BaseModel)
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

        wPhaseIdx = Indices::wPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,

        twoPhases = Indices::twoPhases,
        wPhaseOnly  = Indices::wPhaseOnly,
        gPhaseOnly  = Indices::gPhaseOnly,
    };

    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;


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
            for(const auto& vertex: vertices(this->gridView_()))
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

        for (const auto& element : elements(this->gridView_())) {
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
        ScalarField *poro = writer.allocateManagedBuffer(numDofs);
        ScalarField *permXX = writer.allocateManagedBuffer(numDofs);
        ScalarField *permYY = writer.allocateManagedBuffer(numDofs);
        ScalarField *permZZ = writer.allocateManagedBuffer(numDofs);
        VectorField *velocityW = writer.template allocateManagedBuffer<double, dim>(numDofs);
        VectorField *velocityG = writer.template allocateManagedBuffer<double, dim>(numDofs);
        ImplicitVelocityOutput<TypeTag> velocityOutput(this->problem_());

        if (velocityOutput.enableOutput()) // check if velocity output is demanded
        {
            // initialize velocity fields
            for (unsigned int i = 0; i < numDofs; ++i)
            {
                (*velocityW)[i] = Scalar(0);
                (*velocityG)[i] = Scalar(0);
            }
        }

        unsigned numElements = this->gridView_().size(0);
        ScalarField *rank = writer.allocateManagedBuffer (numElements);

        for (const auto& element : elements(this->gridView_()))
        {
            int idx = this->dofMapper().index(element);
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

                (*poro)[globalIdx] = elemVolVars[scvIdx].porosity();
                (*phasePresence)[globalIdx] = staticDat_[globalIdx].phasePresence;

                FieldMatrix K = this->problem_().spatialParams().intrinsicPermeability(element,
                                                                                       fvGeometry,
                                                                                       scvIdx);
                (*permXX)[globalIdx] = K[0][0];
                (*permYY)[globalIdx] = K[1][1];
                (*permZZ)[globalIdx] = K[2][2];
            }

            // velocity output
            velocityOutput.calculateVelocity(*velocityW, elemVolVars, fvGeometry, element, wPhaseIdx);
            velocityOutput.calculateVelocity(*velocityG, elemVolVars, fvGeometry, element, gPhaseIdx);

        }

        writer.attachDofData(*saturation[wPhaseIdx], "sw", isBox);
        writer.attachDofData(*saturation[gPhaseIdx], "sg", isBox);
        writer.attachDofData(*pressure[wPhaseIdx], "pw", isBox);
        writer.attachDofData(*pressure[gPhaseIdx], "pg", isBox);
        writer.attachDofData(*density[wPhaseIdx], "rhow", isBox);
        writer.attachDofData(*density[gPhaseIdx], "rhog", isBox);

        writer.attachDofData(*mobility[wPhaseIdx], "MobW", isBox);
        writer.attachDofData(*mobility[gPhaseIdx], "MobG", isBox);

        writer.attachDofData(*viscosity[wPhaseIdx], "ViscosW", isBox);
        writer.attachDofData(*viscosity[gPhaseIdx], "ViscosG", isBox);

        writer.attachDofData(*poro, "porosity", isBox);
        writer.attachDofData(*permXX, "permeabilityXX", isBox);
        writer.attachDofData(*permYY, "permeabilityYY", isBox);
        writer.attachDofData(*permZZ, "permeabilityZZ", isBox);
        writer.attachDofData(*phasePresence, "phase presence", isBox);

        if (velocityOutput.enableOutput()) // check if velocity output is demanded
        {
            writer.attachDofData(*velocityW,  "velocityW", isBox, dim);
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
        int succeeded;
        try {
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
            updateFailed();

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

        globalSol[globalIdx][pressureIdx] = volVars.fluidState().pressure(gPhaseIdx);

        // check if a primary var switch is necessary
        if (phasePresence == twoPhases)
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
                newPhasePresence = wPhaseOnly;

                globalSol[globalIdx][switch1Idx] = volVars.fluidState().temperature();
            }
            else if (volVars.saturation(wPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // water phase disappears
                std::cout << "Wetting phase disappears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", sw: "
                          << volVars.saturation(wPhaseIdx) << std::endl;
                newPhasePresence = gPhaseOnly;

                globalSol[globalIdx][switch1Idx] = volVars.fluidState().temperature();
            }

        }
        else if (phasePresence == wPhaseOnly)
        {
            Scalar temp = volVars.fluidState().temperature();
            Scalar pg = volVars.fluidState().pressure(gPhaseIdx); //TODO: wPhaseIndex für Brooks Corey?
            Scalar tempVap = volVars.vaporTemperature();

            // if the the temperature would be larger than
            // the vapor temperature at the given pressure, gas phase appears
            if (temp >= tempVap)
            {
                wouldSwitch = true;
                // gas phase appears
                std::cout << "gas phase appears at vertex " << globalIdx
                          << ", coordinates: " << globalPos  << std::endl;

               newPhasePresence = twoPhases;
               globalSol[globalIdx][switch1Idx] = 0.9999; //wetting phase saturation
            }
        }


        else if (phasePresence == gPhaseOnly)
        {

            Scalar temp = volVars.fluidState().temperature();
            Scalar pg = volVars.fluidState().pressure(gPhaseIdx); //TODO: wPhaseIndex für Brooks Corey?
            Scalar tempVap = volVars.vaporTemperature();

            if (temp < tempVap)
            {
                wouldSwitch = true;
                // wetting phase appears
                std::cout << "wetting phase appears at vertex " << globalIdx
                          << ", coordinates: " << globalPos  << std::endl;


               newPhasePresence = twoPhases;
               globalSol[globalIdx][switch1Idx] = 0.0001; //arbitrary small value
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
