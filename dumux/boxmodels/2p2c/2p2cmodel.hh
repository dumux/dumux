// $Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch     *
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

/*!
* \file
*
* \brief Adaption of the BOX scheme to the two-phase two-component flow model.
*/

#ifndef DUMUX_2P2C_MODEL_HH
#define DUMUX_2P2C_MODEL_HH

#include "2p2cproperties.hh"
#include "2p2clocalresidual.hh"
#include "2p2cproblem.hh"

namespace Dumux
{
/*!
 * \ingroup BoxModels
 * \defgroup TwoPTwoCModel Two-phase two-component box model
 */

/*!
 * \ingroup TwoPTwoCModel
 * \brief Adaption of the BOX scheme to the two-phase two-component flow model.
 *
 * This model implements two-phase two-component flow of two compressible and
 * partially miscible fluids \f$\alpha \in \{ w, n \}\f$ composed of the two components
 * \f$\kappa \in \{ w, a \}\f$. The standard multiphase Darcy
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
 - \sum_\alpha q_\alpha^\kappa = 0 \qquad \kappa \in \{w, a\} \, ,
 \alpha \in \{w, g\}
 \f}
 *
 * This is discretized using a fully-coupled vertex
 * centered finite volume (box) scheme as spatial and
 * the implicit Euler method as temporal discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$ and \f$X^\kappa_w + X^\kappa_n = 1\f$, the number of
 * unknowns can be reduced to two.
 * The used primary variables are, like in the two-phase model, either \f$p_w\f$ and \f$S_n\f$
 * or \f$p_n\f$ and \f$S_w\f$. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * TwoPTwoCIndices::pWsN or TwoPTwoCIndices::pNsW. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 * Moreover, the second primary variable depends on the phase state, since a
 * primary variable switch is included. The phase state is stored for all nodes
 * of the system. Following cases can be distinguished:
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
class TwoPTwoCModel: public BoxModel<TypeTag>
{
    typedef TwoPTwoCModel<TypeTag> ThisType;
    typedef BoxModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementBoundaryTypes)) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxVariables)) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexMapper)) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementMapper)) ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices)) Indices;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
        numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)),

        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        contiLEqIdx = Indices::contiLEqIdx,
        contiGEqIdx = Indices::contiGEqIdx,

        lPhaseIdx = Indices::lPhaseIdx,
        gPhaseIdx = Indices::gPhaseIdx,

        lCompIdx = Indices::lCompIdx,
        gCompIdx = Indices::gCompIdx,

        lPhaseOnly = Indices::lPhaseOnly,
        gPhaseOnly = Indices::gPhaseOnly,
        bothPhases = Indices::bothPhases,

        plSg = TwoPTwoCFormulation::plSg,
        pgSl = TwoPTwoCFormulation::pgSl,
        formulation = GET_PROP_VALUE(TypeTag, PTAG(Formulation))
    };

    typedef TwoPTwoCFluidState<TypeTag> FluidState;

    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    static const Scalar mobilityUpwindAlpha =
            GET_PROP_VALUE(TypeTag, PTAG(MobilityUpwindAlpha));

public:
    /*!
     * \brief Initialize the static data with the initial solution.
     *
     * \param problem The problem to be solved
     */
    void init(Problem &problem)
    {
        ParentType::init(problem);

        staticVertexDat_.resize(this->gridView_().size(dim));

        setSwitched_(false);

        VertexIterator it = this->gridView_().template begin<dim> ();
        const VertexIterator &endit = this->gridView_().template end<dim> ();
        for (; it != endit; ++it)
        {
            int globalIdx = this->dofMapper().map(*it);
            const GlobalPosition &globalPos = it->geometry().corner(0);

            // initialize phase presence
            staticVertexDat_[globalIdx].phasePresence
                = this->problem_().initialPhasePresence(*it, globalIdx,
                            globalPos);
            staticVertexDat_[globalIdx].wasSwitched = false;

            staticVertexDat_[globalIdx].oldPhasePresence
                    = staticVertexDat_[globalIdx].phasePresence;
        }
    }

    /*!
     * \brief Compute the total storage inside one phase of all
     *        conservation quantities.
     *
     * \param dest Contains the storage of each component for one phase
     * \param phaseIdx The phase index
     */
    void globalPhaseStorage(PrimaryVariables &dest, int phaseIdx)
    {
        dest = 0;

        ElementIterator elemIt = this->gridView_().template begin<0>();
        const ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            this->localResidual().evalPhaseStorage(*elemIt, phaseIdx);

            for (int i = 0; i < elemIt->template count<dim>(); ++i)
                dest += this->localResidual().residual(i);
        };

        if (this->gridView_().comm().size() > 1)
        	dest = this->gridView_().comm().sum(dest);
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
        /*this->localJacobian().updateStaticData(this->curSolFunction(),
         this->prevSolFunction());
         */
    };

    /*!
     * \brief Returns the relative weight of a primary variable for
     *        calculating relative errors.
     *
     * \param globalVertexIdx The global vertex index
     * \param pvIdx The primary variable index
     */
    Scalar primaryVarWeight(int globalVertexIdx, int pvIdx) const
    {
        if (Indices::pressureIdx == pvIdx)
            return std::min(1.0/this->prevSol()[globalVertexIdx][pvIdx], 1.0);
        return 1;
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
        return oldSol ? staticVertexDat_[globalVertexIdx].oldPhasePresence
                : staticVertexDat_[globalVertexIdx].phasePresence;
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

        // create the required scalar fields
        unsigned numVertices = this->problem_().gridView().size(dim);

        ScalarField *Sg = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *Sl = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *pg = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *pl = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *pc = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *rhoL =
                writer.template createField<Scalar, 1> (numVertices);
        ScalarField *rhoG =
                writer.template createField<Scalar, 1> (numVertices);
        ScalarField *mobL =
                writer.template createField<Scalar, 1> (numVertices);
        ScalarField *mobG =
                writer.template createField<Scalar, 1> (numVertices);
        ScalarField *phasePresence = writer.template createField<Scalar, 1> (
                numVertices);
        ScalarField *massFrac[numPhases][numComponents];
        for (int i = 0; i < numPhases; ++i)
            for (int j = 0; j < numComponents; ++j)
                massFrac[i][j] = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *temperature = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *poro = writer.template createField<Scalar, 1>(numVertices);

#ifdef VELOCITY_OUTPUT // check if velocity output is demanded
        ScalarField *velocityX = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *velocityY = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *velocityZ = writer.template createField<Scalar, 1>(numVertices);
//        Scalar maxV=0.; // variable to store the maximum face velocity

        // initialize velocity fields
        Scalar boxSurface[numVertices];
        for (int i = 0; i < numVertices; ++i)
        {
            (*velocityX)[i] = 0;
            if (dim > 1)
            (*velocityY)[i] = 0;
            if (dim > 2)
            (*velocityZ)[i] = 0;
            boxSurface[i] = 0.0; // initialize the boundary surface of the fv-boxes
        }
#endif

        unsigned numElements = this->gridView_().size(0);
        ScalarField *rank =
                writer.template createField<Scalar, 1> (numElements);

        FVElementGeometry fvElemGeom;
        VolumeVariables volVars;

        ElementIterator elemIt = this->gridView_().template begin<0>();
        ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            int idx = this->problem_().elementMapper().map(*elemIt);
            (*rank)[idx] = this->gridView_().comm().rank();
            fvElemGeom.update(this->gridView_(), *elemIt);

            int numVerts = elemIt->template count<dim> ();
            for (int i = 0; i < numVerts; ++i)
            {
                int globalIdx = this->vertexMapper().map(*elemIt, i, dim);
                volVars.update(sol[globalIdx],
                               this->problem_(),
                               *elemIt,
                               fvElemGeom,
                               i,
                               false);
                (*Sg)[globalIdx] = volVars.saturation(gPhaseIdx);
                (*Sl)[globalIdx] = volVars.saturation(lPhaseIdx);
                (*pg)[globalIdx] = volVars.pressure(gPhaseIdx);
                (*pl)[globalIdx] = volVars.pressure(lPhaseIdx);
                (*pc)[globalIdx] = volVars.capillaryPressure();
                (*rhoL)[globalIdx] = volVars.fluidState().density(lPhaseIdx);
                (*rhoG)[globalIdx] = volVars.fluidState().density(gPhaseIdx);
                (*mobL)[globalIdx] = volVars.mobility(lPhaseIdx);
                (*mobG)[globalIdx] = volVars.mobility(gPhaseIdx);
                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                    for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    {
                        (*massFrac[phaseIdx][compIdx])[globalIdx]
                                = volVars.fluidState().massFrac(phaseIdx,
                                        compIdx);

                        Valgrind::CheckDefined(
                            (*massFrac[phaseIdx][compIdx])[globalIdx][0]);
                    }
                (*poro)[globalIdx] = volVars.porosity();
                (*temperature)[globalIdx] = volVars.temperature();
                (*phasePresence)[globalIdx]
                        = staticVertexDat_[globalIdx].phasePresence;
            };

#ifdef VELOCITY_OUTPUT // check if velocity output is demanded
            // In the box method, the velocity is evaluated on the FE-Grid. However, to get an
            // average apparent velocity at the vertex, all contributing velocities have to be interpolated.
            GlobalPosition velocity(0.);
            ElementVolumeVariables elemVolVars;

            elemVolVars.update(this->problem_(),
                               *elemIt,
                               fvElemGeom,
                               false /* isOldSol? */);

            // loop over the phases
            for (int faceIdx = 0; faceIdx< fvElemGeom.numEdges; faceIdx++)
            {
                //prepare the flux calculations (set up and prepare geometry, FE gradients)
                FluxVariables fluxDat(this->problem_(),
                                 *elemIt,
                                 fvElemGeom,
                                 faceIdx,
                                 elemVolVars);

                // choose phase of interest. Alternatively, a loop over all phases would be possible.
                int phaseIdx = gPhaseIdx;

                // get darcy velocity
                velocity = fluxDat.Kmvp(phaseIdx); // mind the sign: vDarcy = kf grad p

                // up+downstream mobility
                const VolumeVariables &up = elemVolVars[fluxDat.upstreamIdx(phaseIdx)];
                const VolumeVariables &down = elemVolVars[fluxDat.downstreamIdx(phaseIdx)];
                Scalar scvfArea = fluxDat.face().normal.two_norm(); //get surface area to weight velocity at the IP with the surface area
                velocity *= (mobilityUpwindAlpha*up.mobility(phaseIdx) + (1-mobilityUpwindAlpha)*down.mobility(phaseIdx))* scvfArea;

                int vertIIdx = this->problem_().vertexMapper().map(*elemIt,
                        fluxDat.face().i,
                        dim);
                int vertJIdx = this->problem_().vertexMapper().map(*elemIt,
                        fluxDat.face().j,
                        dim);
                // add surface area for weighting purposes
                boxSurface[vertIIdx] += scvfArea;
                boxSurface[vertJIdx] += scvfArea;

                // Add velocity to upstream and downstream vertex.
                // Beware: velocity has to be substracted because of the (wrong) sign of vDarcy
                (*velocityX)[vertIIdx] -= velocity[0];
                (*velocityX)[vertJIdx] -= velocity[0];
                if (dim >= 2)
                {
                    (*velocityY)[vertIIdx] -= velocity[1];
                    (*velocityY)[vertJIdx] -= velocity[1];
                }
                if (dim == 3)
                {
                    (*velocityZ)[vertIIdx] -= velocity[2];
                    (*velocityZ)[vertJIdx] -= velocity[2];
                }
            }
#endif
        }

#ifdef VELOCITY_OUTPUT // check if velocity output is demanded
        // normalize the velocities at the vertices
        for (int i = 0; i < numVertices; ++i)
        {
            (*velocityX)[i] /= boxSurface[i];
            if (dim >= 2)
            (*velocityY)[i] /= boxSurface[i];
            if (dim == 3)
            (*velocityZ)[i] /= boxSurface[i];
        }
#endif

        writer.addVertexData(Sg, "Sg");
        writer.addVertexData(Sl, "Sl");
        writer.addVertexData(pg, "pg");
        writer.addVertexData(pl, "pl");
        writer.addVertexData(pc, "pc");
        writer.addVertexData(rhoL, "rhoL");
        writer.addVertexData(rhoG, "rhoG");
        writer.addVertexData(mobL, "mobL");
        writer.addVertexData(mobG, "mobG");
        for (int i = 0; i < numPhases; ++i)
        {
            for (int j = 0; j < numComponents; ++j)
            {
                std::string name = (boost::format("X_%s%s")
                        % ((i == lPhaseIdx) ? "l" : "g")
                        % FluidSystem::componentName(j)).str();
                writer.addVertexData(massFrac[i][j], name.c_str());
            }
        }
        writer.addVertexData(poro, "porosity");
        writer.addVertexData(temperature, "temperature");
        writer.addVertexData(phasePresence, "phase presence");

#ifdef VELOCITY_OUTPUT // check if velocity output is demanded
        writer.addVertexData(velocityX, "Vx");
        if (dim >= 2)
        writer.addVertexData(velocityY, "Vy");
        if (dim == 3)
        writer.addVertexData(velocityZ, "Vz");
#endif
        writer.addCellData(rank, "process rank");
    }

    /*!
     * \brief Write the current solution to a restart file.
     *
     * \param outStream The output stream of one vertex for the restart file
     * \param vert The vertex
     */
    void serializeEntity(std::ostream &outStream, const Vertex &vert)
    {
        // write primary variables
        ParentType::serializeEntity(outStream, vert);

        int vertIdx = this->dofMapper().map(vert);
        if (!outStream.good())
            DUNE_THROW(Dune::IOError, "Could not serialize vertex " << vertIdx);

        outStream << staticVertexDat_[vertIdx].phasePresence << " ";
    }

    /*!
     * \brief Reads the current solution for a vertex from a restart
     *        file.
     *
     * \param inStream The input stream of one vertex from the restart file
     * \param vert The vertex
     */
    void deserializeEntity(std::istream &inStream, const Vertex &vert)
    {
        // read primary variables
        ParentType::deserializeEntity(inStream, vert);

        // read phase presence
        int vertIdx = this->dofMapper().map(vert);
        if (!inStream.good())
            DUNE_THROW(Dune::IOError,
                       "Could not deserialize vertex " << vertIdx);

        inStream >> staticVertexDat_[vertIdx].phasePresence;
        staticVertexDat_[vertIdx].oldPhasePresence
                = staticVertexDat_[vertIdx].phasePresence;

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

        for (unsigned i = 0; i < staticVertexDat_.size(); ++i)
            staticVertexDat_[i].visited = false;

        FVElementGeometry fvElemGeom;
        static VolumeVariables volVars;
        ElementIterator it = this->gridView_().template begin<0> ();
        const ElementIterator &endit = this->gridView_().template end<0> ();
        for (; it != endit; ++it)
        {
            fvElemGeom.update(this->gridView_(), *it);
            for (int i = 0; i < fvElemGeom.numVertices; ++i)
            {
                int globalIdx = this->vertexMapper().map(*it, i, dim);

                if (staticVertexDat_[globalIdx].visited)
                    continue;

                staticVertexDat_[globalIdx].visited = true;
                volVars.update(curGlobalSol[globalIdx],
                               this->problem_(),
                               *it,
                               fvElemGeom,
                               i,
                               false);
                const GlobalPosition &global = it->geometry().corner(i);
                if (primaryVarSwitch_(curGlobalSol,
                                      volVars,
                                      globalIdx,
                                      global))
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

    /*!
     * \brief Reset the current phase presence of all vertices to the old one.
     *
     * This is done after an update failed.
     */
    void resetPhasePresence_()
    {
        int numVertices = this->gridView_().size(dim);
        for (int i = 0; i < numVertices; ++i)
        {
            staticVertexDat_[i].phasePresence
                    = staticVertexDat_[i].oldPhasePresence;
            staticVertexDat_[i].wasSwitched = false;
        }
    }

    /*!
     * \brief Set the old phase of all verts state to the current one.
     */
    void updateOldPhasePresence_()
    {
        int numVertices = this->gridView_().size(dim);
        for (int i = 0; i < numVertices; ++i)
        {
            staticVertexDat_[i].oldPhasePresence
                    = staticVertexDat_[i].phasePresence;
            staticVertexDat_[i].wasSwitched = false;
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
        int phasePresence = staticVertexDat_[globalIdx].phasePresence;
        int newPhasePresence = phasePresence;

        // check if a primary var switch is necessary
        if (phasePresence == gPhaseOnly)
        {
            // calculate mole fraction in the hypothetic liquid phase
            Scalar xll = volVars.fluidState().moleFrac(lPhaseIdx, lCompIdx);
            Scalar xlg = volVars.fluidState().moleFrac(lPhaseIdx, gCompIdx);

            Scalar xlMax = 1.0;
            if (xll + xlg > xlMax)
                wouldSwitch = true;
            if (staticVertexDat_[globalIdx].wasSwitched)
                xlMax *= 1.02;

            // if the sum of the mole fractions would be larger than
            // 100%, liquid phase appears
            if (xll + xlg > xlMax)
            {
                // liquid phase appears
                std::cout << "liquid phase appears at vertex " << globalIdx
                        << ", coordinates: " << globalPos << ", xll + xlg: "
                        << xll + xlg << std::endl;
                newPhasePresence = bothPhases;
                if (formulation == pgSl)
                    globalSol[globalIdx][switchIdx] = 0.0;
                else if (formulation == plSg)
                    globalSol[globalIdx][switchIdx] = 1.0;
            };
        }
        else if (phasePresence == lPhaseOnly)
        {
            // calculate fractions of the partial pressures in the
            // hypothetic gas phase
            Scalar xgl = volVars.fluidState().moleFrac(gPhaseIdx, lCompIdx);
            Scalar xgg = volVars.fluidState().moleFrac(gPhaseIdx, gCompIdx);

            Scalar xgMax = 1.0;
            if (xgl + xgg > xgMax)
                wouldSwitch = true;
            if (staticVertexDat_[globalIdx].wasSwitched)
                xgMax *= 1.02;

            // if the sum of the mole fractions would be larger than
            // 100%, gas phase appears
            if (xgl + xgg > xgMax)
            {
                // gas phase appears
                std::cout << "gas phase appears at vertex " << globalIdx
                        << ", coordinates: " << globalPos << ", x_gl + x_gg: "
                        << xgl + xgg << std::endl;
                newPhasePresence = bothPhases;
                if (formulation == pgSl)
                    globalSol[globalIdx][switchIdx] = 0.999;
                else if (formulation == plSg)
                    globalSol[globalIdx][switchIdx] = 0.001;
            }
        }
        else if (phasePresence == bothPhases)
        {
            Scalar Smin = 0.0;
            if (staticVertexDat_[globalIdx].wasSwitched)
                Smin = -0.01;

            if (volVars.saturation(gPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // gas phase disappears
                std::cout << "Gas phase disappears at vertex " << globalIdx
                        << ", coordinates: " << globalPos << ", Sg: "
                        << volVars.saturation(gPhaseIdx) << std::endl;
                newPhasePresence = lPhaseOnly;

                globalSol[globalIdx][switchIdx]
                        = volVars.fluidState().massFrac(lPhaseIdx, gCompIdx);
            }
            else if (volVars.saturation(lPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // liquid phase disappears
                std::cout << "Liquid phase disappears at vertex " << globalIdx
                        << ", coordinates: " << globalPos << ", Sl: "
                        << volVars.saturation(lPhaseIdx) << std::endl;
                newPhasePresence = gPhaseOnly;

                globalSol[globalIdx][switchIdx]
                        = volVars.fluidState().massFrac(gPhaseIdx, lCompIdx);
            }
        }

        staticVertexDat_[globalIdx].phasePresence = newPhasePresence;
        staticVertexDat_[globalIdx].wasSwitched = wouldSwitch;
        return phasePresence != newPhasePresence;
    }

    // parameters given in constructor
    std::vector<StaticVars> staticVertexDat_;
    bool switchFlag_;
};

}

#include "2p2cpropertydefaults.hh"

#endif
