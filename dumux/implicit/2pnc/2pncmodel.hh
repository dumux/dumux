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
* \brief Adaption of the fully implicit box scheme to the two-phase n-component flow model.
*/

#ifndef DUMUX_2PNC_MODEL_HH
#define DUMUX_2PNC_MODEL_HH

#include <dune/common/version.hh>
// #include <dumux/material/constants.hh>

#include "2pncproperties.hh"
#include "2pncindices.hh"
#include "2pnclocalresidual.hh"
#include <dumux/implicit/common/implicitvelocityoutput.hh>

namespace Dumux
{
/*!
 * \ingroup TwoPNCModel
 * \brief Adaption of the fully implicit scheme to the
 *        two-phase n-component fully implicit model.
 *
 * This model implements two-phase n-component flow of two compressible and
 * partially miscible fluids \f$\alpha \in \{ w, n \}\f$ composed of the n components
 * \f$\kappa \in \{ w, a,\cdots \}\f$. The standard multiphase Darcy
 * approach is used as the equation for the conservation of momentum:
 * \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 \left(\text{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 * \f]
 *
 * By inserting this into the equations for the conservation of the
 * components, one gets one transport equation for each component
 * \f{eqnarray}
 && \phi \frac{\partial (\sum_\alpha \varrho_\alpha X_\alpha^\kappa S_\alpha )}
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
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme (this is not done for 2pnc approach yet, however possible) as  
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
 * TwoPTwoCIndices::pWsN or TwoPTwoCIndices::pNsW. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 *
 * Moreover, the second primary variable depends on the phase state, since a
 * primary variable switch is included. The phase state is stored for all nodes
 * of the system. The model is uses mole fractions.
 *Following cases can be distinguished:
 * <ul>
 *  <li> Both phases are present: The saturation is used (either \f$S_n\f$ or \f$S_w\f$, dependent on the chosen <tt>Formulation</tt>),
 *      as long as \f$ 0 < S_\alpha < 1\f$</li>.
 *  <li> Only wetting phase is present: The mass fraction of, e.g., air in the wetting phase \f$X^a_w\f$ is used,
 *      as long as the maximum mass fraction is not exceeded (\f$X^a_w<X^a_{w,max}\f$)</li>
 *  <li> Only non-wetting phase is present: The mass fraction of, e.g., water in the non-wetting phase, \f$X^w_n\f$, is used,
 *      as long as the maximum mass fraction is not exceeded (\f$X^w_n<X^w_{n,max}\f$)</li>
 * </ul>
 */

template<class TypeTag>
class TwoPNCModel: public GET_PROP_TYPE(TypeTag, BaseModel)
{
    typedef TwoPNCModel<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, BaseModel) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef Dumux::Constants<Scalar> Constant;

    enum {
            dim = GridView::dimension,
            dimWorld = GridView::dimensionworld,

            numEq = GET_PROP_VALUE(TypeTag, NumEq),
            numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
            numComponents = GET_PROP_VALUE(TypeTag, NumComponents),
            numMajorComponents = GET_PROP_VALUE(TypeTag, NumMajorComponents),

            pressureIdx = Indices::pressureIdx,
            switchIdx = Indices::switchIdx,

            wPhaseIdx = Indices::wPhaseIdx,
            nPhaseIdx = Indices::nPhaseIdx,

            wCompIdx = FluidSystem::wCompIdx,
            nCompIdx = FluidSystem::nCompIdx,

            wPhaseOnly = Indices::wPhaseOnly,
            nPhaseOnly = Indices::nPhaseOnly,
            bothPhases = Indices::bothPhases,

            plSg = TwoPNCFormulation::plSg,
            pgSl = TwoPNCFormulation::pgSl,
            formulation = GET_PROP_VALUE(TypeTag, Formulation),
            useElectrochem = GET_PROP_VALUE(TypeTag, useElectrochem)
          };

    typedef CompositionalFluidState<Scalar, FluidSystem> FluidState;

    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldMatrix<CoordScalar, dimWorld, dimWorld> Tensor;
    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;

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

        ElementIterator eIt = this->gridView_().template begin<0>();
        ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; eIt != elemEndIt; ++eIt)
        {
            if (!isBox) // i.e. cell-centered discretization
            {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
                int eIdxGlobal = this->dofMapper().index(*eIt);
#else
                int eIdxGlobal = this->dofMapper().map(*eIt);
#endif
                const GlobalPosition &globalPos = eIt->geometry().center();

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
            VertexIterator vIt = this->gridView_().template begin<dim> ();
            const VertexIterator &vEndIt = this->gridView_().template end<dim> ();
            for (; vIt != vEndIt; ++vIt)
            {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
                int vIdxGlobal = this->dofMapper().index(*vIt);
#else
                int vIdxGlobal = this->dofMapper().map(*vIt);
#endif
                const GlobalPosition &globalPos = vIt->geometry().corner(0);

                // initialize phase presence
                staticDat_[vIdxGlobal].phasePresence
                    = this->problem_().initialPhasePresence(*vIt, vIdxGlobal,
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
    void globalPhaseStorage(PrimaryVariables &storage, int phaseIdx)
    {
        storage = 0;

        ElementIterator eIt = this->gridView_().template begin<0>();
        const ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; eIt != elemEndIt; ++eIt)
        {
            if(eIt->partitionType() == Dune::InteriorEntity)
            {
              this->localResidual().evalPhaseStorage(*eIt, phaseIdx); 
              
              for (unsigned int i = 0; i < this->localResidual().storageTerm().size(); ++i)
                    storage += this->localResidual().storageTerm()[i];
            }
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
     * \brief Returns the phase presence of the current or the old solution of a vertex.
     *
     * \param globalVertexIdx The global vertex index
     * \param oldSol Evaluate function with solution of current or previous time step
     */
    int phasePresence(int globalVertexIdx, bool oldSol) const
    {
        return oldSol ? staticDat_[globalVertexIdx].oldPhasePresence
                : staticDat_[globalVertexIdx].phasePresence;
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

        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<double, dim> > VectorField;

        // get the number of degrees of freedom
        unsigned numDofs = this->numDofs();

        // create the required scalar fields

        ScalarField *currentDensity = writer.allocateManagedBuffer (numDofs);
        ScalarField *reactionSourceH2O = writer.allocateManagedBuffer (numDofs);
        ScalarField *reactionSourceO2 = writer.allocateManagedBuffer (numDofs);

        ScalarField *Sg            = writer.allocateManagedBuffer (numDofs);
        ScalarField *Sl            = writer.allocateManagedBuffer (numDofs);
        ScalarField *pg            = writer.allocateManagedBuffer (numDofs);
        ScalarField *pl            = writer.allocateManagedBuffer (numDofs);
        ScalarField *pc            = writer.allocateManagedBuffer (numDofs);
        ScalarField *rhoL          = writer.allocateManagedBuffer (numDofs);
        ScalarField *rhoG          = writer.allocateManagedBuffer (numDofs);
        ScalarField *mobL          = writer.allocateManagedBuffer (numDofs);
        ScalarField *mobG          = writer.allocateManagedBuffer (numDofs);
        ScalarField *temperature   = writer.allocateManagedBuffer (numDofs);
        ScalarField *poro          = writer.allocateManagedBuffer (numDofs);
        ScalarField *boxVolume     = writer.allocateManagedBuffer (numDofs);
        VectorField *velocityN = writer.template allocateManagedBuffer<double, dim>(numDofs);
        VectorField *velocityW = writer.template allocateManagedBuffer<double, dim>(numDofs);
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

        ScalarField *moleFraction[numPhases][numComponents];
        for (int i = 0; i < numPhases; ++i)
            for (int j = 0; j < numComponents; ++j)
                moleFraction[i][j] = writer.allocateManagedBuffer(numDofs);

        ScalarField *molarity[numComponents];
        for (int j = 0; j < numComponents ; ++j)
            molarity[j] = writer.allocateManagedBuffer(numDofs);

        ScalarField *Perm[dim];
        for (int j = 0; j < dim; ++j) //Permeability only in main directions xx and yy
            Perm[j] = writer.allocateManagedBuffer(numDofs);

        *boxVolume = 0;

        unsigned numElements = this->gridView_().size(0);
        ScalarField *rank =
                writer.allocateManagedBuffer (numElements);

        FVElementGeometry fvGeometry;
        VolumeVariables volVars;
        ElementVolumeVariables elemVolVars;

        ElementIterator eIt = this->gridView_().template begin<0>();
        ElementIterator eEndIt = this->gridView_().template end<0>();
        for (; eIt != eEndIt; ++eIt)
        {
            int idx = this->problem_().elementMapper().map(*eIt);
            (*rank)[idx] = this->gridView_().comm().rank();
            fvGeometry.update(this->gridView_(), *eIt);

            elemVolVars.update(this->problem_(),
                               *eIt,
                               fvGeometry,
                               false /* oldSol? */);

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
            int numVerts = eIt->subEntities(dim);
#else
            int numVerts = eIt->template count<dim> ();
#endif
            for (int i = 0; i < numVerts; ++i)
            {
                int globalIdx = this->vertexMapper().map(*eIt, i, dim);
                volVars.update(sol[globalIdx],
                               this->problem_(),
                               *eIt,
                               fvGeometry,
                               i,
                               false);

                GlobalPosition globalPos = fvGeometry.subContVol[i].global;
                //Scalar zmax = this->problem_().bBoxMax()[dim-1];
                (*Sg)[globalIdx]             = volVars.saturation(nPhaseIdx);
                (*Sl)[globalIdx]             = volVars.saturation(wPhaseIdx);
                (*pg)[globalIdx]             = volVars.pressure(nPhaseIdx);
                (*pl)[globalIdx]             = volVars.pressure(wPhaseIdx);
                (*pc)[globalIdx]             = volVars.capillaryPressure();
                (*rhoL)[globalIdx]           = volVars.fluidState().density(wPhaseIdx);
                (*rhoG)[globalIdx]           = volVars.fluidState().density(nPhaseIdx);
                (*mobL)[globalIdx]           = volVars.mobility(wPhaseIdx);
                (*mobG)[globalIdx]           = volVars.mobility(nPhaseIdx);
                (*boxVolume)[globalIdx]     += fvGeometry.subContVol[i].volume;
                (*poro)[globalIdx]           = volVars.porosity();

                if(useElectrochem){

                    //for electrochemistry output we need elemVolVars
                    elemVolVars.update(this->problem_(), *eIt, fvGeometry,/*oldSol=*/false);

                    //reactionSource Output
                    PrimaryVariables source;
                    this->problem_().solDependentSource(source, *eIt, fvGeometry, i, elemVolVars);

                    (*reactionSourceH2O)[globalIdx] = source[wPhaseIdx];
                    (*reactionSourceO2)[globalIdx] = source[numComponents-1];

                    //Current Output
                    (*currentDensity)[globalIdx] = -1.0*source[numComponents-1]*4*Constant::F;
                    //recorrection of the area for output
                    Scalar gridYMin = GET_RUNTIME_PARAM(TypeTag, Scalar, Grid.yMin);
                    Scalar gridYMax = GET_RUNTIME_PARAM(TypeTag, Scalar, Grid.yMax);
                    Scalar nCellsY  = GET_RUNTIME_PARAM(TypeTag, Scalar, Grid.CellsY);

                    Scalar lengthBox= (gridYMax - gridYMin)/nCellsY;

                    (*currentDensity)[globalIdx] = (*currentDensity)[globalIdx]/2*lengthBox/10000; //i in [A/cm^2], cf. paper Oliver
                }

                (*temperature)[globalIdx]    = volVars.temperature();
                
                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                    for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    {
                        (*moleFraction[phaseIdx][compIdx])[globalIdx]= volVars.fluidState().moleFraction(phaseIdx,compIdx);

                        Valgrind::CheckDefined((*moleFraction[phaseIdx][compIdx])[globalIdx]);

                    }
                for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    (*molarity[compIdx])[globalIdx] = (volVars.fluidState().molarity(wPhaseIdx, compIdx));

                Tensor K = perm_(this->problem_().spatialParams().intrinsicPermeability(*eIt, fvGeometry, i));

                for (int j = 0; j<dim; ++j)
                    (*Perm[j])[globalIdx] = K[j][j] /* volVars.permFactor()*/;
            };

            // velocity output
            if(velocityOutput.enableOutput()){
                velocityOutput.calculateVelocity(*velocityW, elemVolVars, fvGeometry, *eIt, wPhaseIdx);
                velocityOutput.calculateVelocity(*velocityN, elemVolVars, fvGeometry, *eIt, nPhaseIdx);
            }

        } // loop over element

        writer.attachVertexData(*Sg, "Sg");
        writer.attachVertexData(*Sl, "Sl");
        writer.attachVertexData(*pg, "pg");
        writer.attachVertexData(*pl, "pl");
        writer.attachVertexData(*pc, "pc");
        writer.attachVertexData(*rhoL, "rhoL");
        writer.attachVertexData(*rhoG, "rhoG");
        writer.attachVertexData(*mobL, "mobL");
        writer.attachVertexData(*mobG, "mobG");
        writer.attachVertexData(*poro, "porosity");
        writer.attachVertexData(*temperature, "temperature");
        writer.attachVertexData(*boxVolume, "boxVolume");

        if(useElectrochem){
            writer.attachVertexData(*reactionSourceH2O, "reactionSourceH2O [mol/(sm^2)]");
            writer.attachVertexData(*reactionSourceO2, "reactionSourceO2 [mol/(sm^2)]");
            writer.attachVertexData(*currentDensity, "currentDensity [A/cm^2]");
        }

        writer.attachVertexData(*Perm[0], "Kxx");
        if (dim >= 2)
            writer.attachVertexData(*Perm[1], "Kyy");
        if (dim == 3)
            writer.attachVertexData(*Perm[2], "Kzz");

        for (int i = 0; i < numPhases; ++i)
        {
            for (int j = 0; j < numComponents; ++j)
            {
                std::ostringstream oss;
                oss << "x"
                    << FluidSystem::componentName(j)
                    << FluidSystem::phaseName(i);
                writer.attachVertexData(*moleFraction[i][j], oss.str().c_str());
            }
        }

        for (int j = 0; j < numComponents; ++j)
        {
            std::ostringstream oss;
            oss << "m^w_"
                << FluidSystem::componentName(j);
            writer.attachVertexData(*molarity[j], oss.str().c_str());
        }

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

        int vertIdx = this->dofMapper().map(entity);
        if (!outStream.good())
            DUNE_THROW(Dune::IOError, "Could not serialize vertex " << vertIdx);

        outStream << staticDat_[vertIdx].phasePresence << " ";
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

        // read phase presence
        int vertIdx = this->dofMapper().map(entity);
        if (!inStream.good())
            DUNE_THROW(Dune::IOError,
                       "Could not deserialize vertex " << vertIdx);

        inStream >> staticDat_[vertIdx].phasePresence;
        staticDat_[vertIdx].oldPhasePresence
                = staticDat_[vertIdx].phasePresence;

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

        FVElementGeometry fvGeometry;
        static VolumeVariables volVars;
        ElementIterator it = this->gridView_().template begin<0> ();
        const ElementIterator &endit = this->gridView_().template end<0> ();
        for (; it != endit; ++it)
        {
            fvGeometry.update(this->gridView_(), *it);
            for (int i = 0; i < fvGeometry.numScv; ++i)
            {
                int globalIdx = this->vertexMapper().map(*it, i, dim);

                if (staticDat_[globalIdx].visited)
                    continue;

                staticDat_[globalIdx].visited = true;
                volVars.update(curGlobalSol[globalIdx],
                               this->problem_(),
                               *it,
                               fvGeometry,
                               i,
                               false);
                const GlobalPosition &global = it->geometry().corner(i);
                if (primaryVarSwitch_(curGlobalSol,
                                      volVars,
                                      globalIdx,
                                      global))
                { wasSwitched = true;
//                std::cout<<"Switch works :) "<<std::endl;
                }
            }
        }

        // make sure that if there was a variable switch in an
        // other partition we will also set the switch flag
        // for our partition.
//        if (this->gridView_().comm().size() > 1)
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
        int numDofs = this->gridView_().size(dim);
        for (int i = 0; i < numDofs; ++i)
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
        int numDofs = this->gridView_().size(dim);
        for (int i = 0; i < numDofs; ++i)
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
                           const VolumeVariables &volVars, int globalIdx,
                           const GlobalPosition &globalPos)
    {

            // evaluate primary variable switch
            bool wouldSwitch = false;
            int phasePresence = staticDat_[globalIdx].phasePresence;
            int newPhasePresence = phasePresence;

            //check if a primary variable switch is necessary
            if (phasePresence == bothPhases)
            {
                Scalar Smin = 0; //saturation threshold
                if (staticDat_[globalIdx].wasSwitched)
                    Smin = -0.01;

                //if saturation of liquid phase is smaller 0 switch
                if (volVars.saturation(wPhaseIdx) <= Smin)
                {
                    wouldSwitch = true;
                    //liquid phase has to disappear
                    std::cout << "Liquid Phase disappears at vertex " << globalIdx
                                << ", coordinated: " << globalPos << ", Sl: "
                                << volVars.saturation(wPhaseIdx) << std::endl;
                    newPhasePresence = nPhaseOnly;

                    //switch not depending on formulation
                    //switch "Sl" to "xgH20"
                    globalSol[globalIdx][switchIdx]
                            = volVars.fluidState().moleFraction(nPhaseIdx, wCompIdx /*H2O*/);

                    //switch all secondary components to mole fraction in gas phase
                    for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
                        globalSol[globalIdx][compIdx] = volVars.fluidState().moleFraction(nPhaseIdx,compIdx);
                }
                //if saturation of gas phase is smaller than 0 switch
                else if (volVars.saturation(nPhaseIdx) <= Smin)
                {
                    wouldSwitch = true;
                    //gas phase has to disappear
                    std::cout << "Gas Phase disappears at vertex " << globalIdx
                                << ", coordinated: " << globalPos << ", Sg: "
                                << volVars.saturation(nPhaseIdx) << std::endl;
                    newPhasePresence = wPhaseOnly;

                    //switch "Sl" to "xlN2"
                    globalSol[globalIdx][switchIdx]
                            = volVars.fluidState().moleFraction(wPhaseIdx, nCompIdx /*N2*/);
                }
            }
            else if (phasePresence == nPhaseOnly)
            {
                Scalar xlmax = 1;
                Scalar sumxl = 0;
                //Calculate sum of mole fractions in the hypothetical liquid phase
                for (int compIdx = 0; compIdx < numComponents; compIdx++)
                {
                    sumxl += volVars.fluidState().moleFraction(wPhaseIdx, compIdx);
                }
                if (sumxl > xlmax)
                    wouldSwitch = true;
                if (staticDat_[globalIdx].wasSwitched)
                    xlmax *=1.02;
                //liquid phase appears if sum is larger than one
                if (sumxl/*sum of mole fractions*/ > xlmax/*1*/)
                {
                    std::cout << "Liquid Phase appears at vertex " << globalIdx
                            << ", coordinated: " << globalPos << ", sumxl: "
                            << sumxl << std::endl;
                    newPhasePresence = bothPhases;

                    //saturation of the liquid phase set to 0.0001 (if formulation pgSl and vice versa)
                    if (formulation == pgSl)
                        globalSol[globalIdx][switchIdx] = 0.0001;
                    else if (formulation == plSg)
                        globalSol[globalIdx][switchIdx] = 0.9999;

                    //switch all secondary components back to liquid mole fraction
                    for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
                        globalSol[globalIdx][compIdx] = volVars.fluidState().moleFraction(wPhaseIdx,compIdx);
                }
            }
            else if (phasePresence == wPhaseOnly)
            {
                Scalar xgmax = 1;
                Scalar sumxg = 0;
                //Calculate sum of mole fractions in the hypothetical liquid phase
                for (int compIdx = 0; compIdx < numComponents; compIdx++)
                {
                    sumxg += volVars.fluidState().moleFraction(nPhaseIdx, compIdx);
                }
                if (sumxg > xgmax)
                    wouldSwitch = true;
                if (staticDat_[globalIdx].wasSwitched)
                    xgmax *=1.02;
                //liquid phase appears if sum is larger than one
                if (sumxg > xgmax)
                {
                    std::cout << "Gas Phase appears at vertex " << globalIdx
                            << ", coordinated: " << globalPos << ", sumxg: "
                            << sumxg << std::endl;
                    newPhasePresence = bothPhases;
                    //saturation of the liquid phase set to 0.9999 (if formulation pgSl and vice versa)
                    if (formulation == pgSl)
                        globalSol[globalIdx][switchIdx] = 0.9999;
                    else if (formulation == plSg)
                        globalSol[globalIdx][switchIdx] = 0.0001;

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

#include "2pncpropertydefaults.hh"

#endif
