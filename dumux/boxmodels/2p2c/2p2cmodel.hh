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

    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;
    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

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

        // check, if velocity output can be used (works only for cubes so far)
        velocityOutput_ = GET_PARAM(TypeTag, bool, EnableVelocityOutput);
        ElementIterator elemIt = this->gridView_().template begin<0>();
        ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            if(elemIt->geometry().type().isCube() == false){
                velocityOutput_ = false;
            };
        };
        if (velocityOutput_ != GET_PARAM(TypeTag, bool, EnableVelocityOutput))
            std::cout << "ATTENTION: Velocity output only works for cubes and is set to false for simplices\n";

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

        massUpwindWeight_ = GET_PARAM(TypeTag, Scalar, MassUpwindWeight);
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
                dest += this->localResidual().storageTerm()[i];
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
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
        typedef Dune::BlockVector<Dune::FieldVector<double, dim> > VectorField;

        // create the required scalar fields
        unsigned numVertices = this->problem_().gridView().size(dim);

        ScalarField *Sg    = writer.allocateManagedBuffer(numVertices);
        ScalarField *Sl    = writer.allocateManagedBuffer(numVertices);
        ScalarField *pg    = writer.allocateManagedBuffer(numVertices);
        ScalarField *pl    = writer.allocateManagedBuffer(numVertices);
        ScalarField *pc    = writer.allocateManagedBuffer(numVertices);
        ScalarField *rhoL  = writer.allocateManagedBuffer(numVertices);
        ScalarField *rhoG  = writer.allocateManagedBuffer(numVertices);
        ScalarField *mobL  = writer.allocateManagedBuffer(numVertices);
        ScalarField *mobG = writer.allocateManagedBuffer(numVertices);
        ScalarField *phasePresence = writer.allocateManagedBuffer(numVertices);
        ScalarField *massFrac[numPhases][numComponents];
        for (int i = 0; i < numPhases; ++i)
            for (int j = 0; j < numComponents; ++j)
                massFrac[i][j] = writer.allocateManagedBuffer(numVertices);
        ScalarField *temperature = writer.allocateManagedBuffer(numVertices);
        ScalarField *poro = writer.allocateManagedBuffer(numVertices);
        ScalarField *cellNum =writer.allocateManagedBuffer (numVertices);
        VectorField *velocityG = writer.template allocateManagedBuffer<double, dim>(numVertices);
        VectorField *velocityL = writer.template allocateManagedBuffer<double, dim>(numVertices);

        if(velocityOutput_) // check if velocity output is demanded
        {
            // initialize velocity fields
            for (int i = 0; i < numVertices; ++i)
            {
                (*velocityG)[i] = Scalar(0);
                (*velocityL)[i] = Scalar(0);
                (*cellNum)[i] = Scalar(0.0);
            }
        }

        unsigned numElements = this->gridView_().size(0);
        ScalarField *rank =
            writer.allocateManagedBuffer (numElements);

        FVElementGeometry fvElemGeom;
        VolumeVariables volVars;

        ElementIterator elemIt = this->gridView_().template begin<0>();
        ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            int idx = this->elementMapper().map(*elemIt);
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
                (*Sg)[globalIdx]    = volVars.saturation(gPhaseIdx);
                (*Sl)[globalIdx]    = volVars.saturation(lPhaseIdx);
                (*pg)[globalIdx]    = volVars.pressure(gPhaseIdx);
                (*pl)[globalIdx]    = volVars.pressure(lPhaseIdx);
                (*pc)[globalIdx]    = volVars.capillaryPressure();
                (*rhoL)[globalIdx]  = volVars.fluidState().density(lPhaseIdx);
                (*rhoG)[globalIdx]  = volVars.fluidState().density(gPhaseIdx);
                (*mobL)[globalIdx]  = volVars.mobility(lPhaseIdx);
                (*mobG)[globalIdx]  = volVars.mobility(gPhaseIdx);
                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                    for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    {
                        (*massFrac[phaseIdx][compIdx])[globalIdx]
                            = volVars.fluidState().massFraction(phaseIdx, compIdx);

                        Valgrind::CheckDefined(
                                               (*massFrac[phaseIdx][compIdx])[globalIdx][0]);
                    }
                (*poro)[globalIdx]  = volVars.porosity();
                (*temperature)[globalIdx] = volVars.temperature();
                (*phasePresence)[globalIdx]
                    = staticVertexDat_[globalIdx].phasePresence;
                if(velocityOutput_)
                {
                    (*cellNum)[globalIdx] += 1;
                }
            };

            if(velocityOutput_)
            {
                // calculate vertex velocities
                GlobalPosition tmpVelocity[numPhases];

                for(int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                {
                 tmpVelocity[phaseIdx]  = Scalar(0.0);
                }

                typedef Dune::BlockVector<Dune::FieldVector<Scalar, dim> > SCVVelocities;
                SCVVelocities scvVelocityL(8), scvVelocityG(8);

                scvVelocityL = 0;
                scvVelocityG = 0;

                ElementVolumeVariables elemVolVars;

                elemVolVars.update(this->problem_(),
                                  *elemIt,
                                  fvElemGeom,
                                  false /* oldSol? */);

                for (int faceIdx = 0; faceIdx< fvElemGeom.numEdges; faceIdx++)
                {

                    FluxVariables fluxDat(this->problem_(),
                                  *elemIt,
                                  fvElemGeom,
                                  faceIdx,
                                  elemVolVars);

                    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                    {

                         // data attached to upstream and the downstream vertices
                        // of the current phase
                        const VolumeVariables up =
                            elemVolVars[fluxDat.upstreamIdx(phaseIdx)];
                        const VolumeVariables dn =
                            elemVolVars[fluxDat.downstreamIdx(phaseIdx)];

                      // local position of integration point
                      const Dune::FieldVector<Scalar, dim>& localPosIP = fvElemGeom.subContVolFace[faceIdx].ipLocal;

                      // Transformation of the global normal vector to normal vector in the reference element
                      const Dune::FieldMatrix<CoordScalar, dim, dim> jacobianT1 = elemIt->geometry().jacobianTransposed(localPosIP);
                      const GlobalPosition globalNormal = fluxDat.face().normal;

                      GlobalPosition localNormal(0);
                      jacobianT1.mv(globalNormal, localNormal);
                        // note only works for cubes
                      const Scalar localArea = pow(2,-(dim-1));

                      localNormal /= localNormal.two_norm();

                      // Get the Darcy velocities. The Darcy velocities are divided by the area of the subcontrolvolume
                      // face in the reference element.
                      massUpwindWeight_ = GET_PARAM(TypeTag, Scalar, MassUpwindWeight);
                      PhasesVector q;
                      q[phaseIdx] = fluxDat.KmvpNormal(phaseIdx)
                                           * (massUpwindWeight_
                                           * up.mobility(phaseIdx)
                                           + (1- massUpwindWeight_)
                                           * dn.mobility(phaseIdx)) / localArea;

                      // transform the normal Darcy velocity into a vector
                      tmpVelocity[phaseIdx] = localNormal;
                      tmpVelocity[phaseIdx] *= q[phaseIdx];

                      if(phaseIdx == lPhaseIdx){
                      scvVelocityL[fluxDat.face().i] += tmpVelocity[phaseIdx];
                      scvVelocityL[fluxDat.face().j] += tmpVelocity[phaseIdx];
                      }
                      else if(phaseIdx == gPhaseIdx){
                      scvVelocityG[fluxDat.face().i] += tmpVelocity[phaseIdx];
                      scvVelocityG[fluxDat.face().j] += tmpVelocity[phaseIdx];
                      }
                   }
                }

                typedef Dune::GenericReferenceElements<Scalar, dim> ReferenceElements;
                const Dune::FieldVector<Scalar, dim>& localPos =
                    ReferenceElements::general(elemIt->geometry().type()).position(0, 0);

                 // get the transposed Jacobian of the element mapping
                const Dune::FieldMatrix<CoordScalar, dim, dim> jacobianT2 = elemIt->geometry().jacobianTransposed(localPos);

                // transform vertex velocities from local to global coordinates
                for (int i = 0; i < numVerts; ++i)
                {
                    int globalIdx = this->vertexMapper().map(*elemIt, i, dim);
                    // calculate the subcontrolvolume velocity by the Piola transformation
                    Dune::FieldVector<CoordScalar, dim> scvVelocity(0);

                    jacobianT2.mtv(scvVelocityL[i], scvVelocity);
                    scvVelocity /= elemIt->geometry().integrationElement(localPos);
                    // add up the wetting phase subcontrolvolume velocities for each vertex
                    (*velocityL)[globalIdx] += scvVelocity;

                    jacobianT2.mtv(scvVelocityG[i], scvVelocity);
                    scvVelocity /= elemIt->geometry().integrationElement(localPos);
                    // add up the nonwetting phase subcontrolvolume velocities for each vertex
                    (*velocityG)[globalIdx] += scvVelocity;
                }
            }
        }
            if(velocityOutput_)
            {
                // divide the vertex velocities by the number of adjacent scvs i.e. cells
                for(int globalIdx = 0; globalIdx<numVertices; ++globalIdx){
                (*velocityL)[globalIdx] /= (*cellNum)[globalIdx];
                (*velocityG)[globalIdx] /= (*cellNum)[globalIdx];
                }
            }


        writer.attachVertexData(*Sg,     "Sg");
        writer.attachVertexData(*Sl,     "Sl");
        writer.attachVertexData(*pg,     "pg");
        writer.attachVertexData(*pl,     "pl");
        writer.attachVertexData(*pc,     "pc");
        writer.attachVertexData(*rhoL,   "rhoL");
        writer.attachVertexData(*rhoG,   "rhoG");
        writer.attachVertexData(*mobL,   "mobL");
        writer.attachVertexData(*mobG,   "mobG");
        for (int i = 0; i < numPhases; ++i)
        {
            for (int j = 0; j < numComponents; ++j)
            {
                std::string name = (boost::format("X_%s%s")
                                    % ((i == lPhaseIdx) ? "l" : "g")
                                    % FluidSystem::componentName(j)).str();
                writer.attachVertexData(*massFrac[i][j], name.c_str());
            }
        }
        writer.attachVertexData(*poro, "porosity");
        writer.attachVertexData(*temperature,    "temperature");
        writer.attachVertexData(*phasePresence,  "phase presence");

        if(velocityOutput_) // check if velocity output is demanded
        {
            writer.attachVertexData(*velocityL,  "velocityL", dim);
            writer.attachVertexData(*velocityG,  "velocityG", dim);
        }
        writer.attachCellData(*rank, "process rank");
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
                {
                    this->jacobianAssembler().markVertexRed(globalIdx);
                    wasSwitched = true;
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
            Scalar xll = volVars.fluidState().moleFraction(lPhaseIdx, lCompIdx);
            Scalar xlg = volVars.fluidState().moleFraction(lPhaseIdx, gCompIdx);

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
            Scalar xgl = volVars.fluidState().moleFraction(gPhaseIdx, lCompIdx);
            Scalar xgg = volVars.fluidState().moleFraction(gPhaseIdx, gCompIdx);

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
                    = volVars.fluidState().massFraction(lPhaseIdx, gCompIdx);
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
                    = volVars.fluidState().massFraction(gPhaseIdx, lCompIdx);
            }
        }

        staticVertexDat_[globalIdx].phasePresence = newPhasePresence;
        staticVertexDat_[globalIdx].wasSwitched = wouldSwitch;
        return phasePresence != newPhasePresence;
    }

protected:
    // parameters given in constructor
    std::vector<StaticVars> staticVertexDat_;
    bool switchFlag_;
    bool velocityOutput_;
    Scalar massUpwindWeight_;
};

}

#include "2p2cpropertydefaults.hh"

#endif
