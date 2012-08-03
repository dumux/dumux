// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch     *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
 * \f{eqnarray*}
 && \phi \frac{\partial (\sum_\alpha \varrho_\alpha X_\alpha^\kappa S_\alpha )}
 {\partial t}
 - \sum_\alpha  \text{div} \left\{ \varrho_\alpha X_\alpha^\kappa
 \frac{k_{r\alpha}}{\mu_\alpha} \mbox{\bf K}
 (\text{grad}\, p_\alpha - \varrho_{\alpha}  \mbox{\bf g}) \right\}
 \nonumber \\ \nonumber \\
 &-& \sum_\alpha \text{div} \left\{{\bf D}_{\alpha, pm}^\kappa \varrho_{\alpha} \text{grad}\, X^\kappa_{\alpha} \right\}
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
 *      as long as the maximum mass fraction is not exceeded \f$(X^a_w<X^a_{w,max})\f$</li>
 *  <li> Only non-wetting phase is present: The mass fraction of, e.g., water in the non-wetting phase, \f$X^w_n\f$, is used,
 *      as long as the maximum mass fraction is not exceeded \f$(X^w_n<X^w_{n,max})\f$</li>
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
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    enum {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx,

        wPhaseOnly = Indices::wPhaseOnly,
        nPhaseOnly = Indices::nPhaseOnly,
        bothPhases = Indices::bothPhases,

        pwSn = TwoPTwoCFormulation::pwSn,
        pnSw = TwoPTwoCFormulation::pnSw,
        formulation = GET_PROP_VALUE(TypeTag, Formulation)
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;
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
        velocityOutput_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddVelocity);
        ElementIterator eIt = this->gridView_().template begin<0>();
        ElementIterator eEndIt = this->gridView_().template end<0>();
        for (; eIt != eEndIt; ++eIt)
        {
            if(eIt->geometry().type().isCube() == false){
                velocityOutput_ = false;
            }
        }
        if (velocityOutput_ != GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddVelocity))
            std::cout << "ATTENTION: Velocity output only works for cubes and is set to false for simplices\n";

        VertexIterator vIt = this->gridView_().template begin<dim> ();
        const VertexIterator &vEndIt = this->gridView_().template end<dim> ();
        for (; vIt != vEndIt; ++vIt)
        {
            int globalIdx = this->dofMapper().map(*vIt);
            const GlobalPosition &globalPos = vIt->geometry().corner(0);

            // initialize phase presence
            staticVertexDat_[globalIdx].phasePresence
                = this->problem_().initialPhasePresence(*vIt, globalIdx,
                                                        globalPos);
            staticVertexDat_[globalIdx].wasSwitched = false;

            staticVertexDat_[globalIdx].oldPhasePresence
                = staticVertexDat_[globalIdx].phasePresence;
        }

        massUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
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

        ElementIterator eIt = this->gridView_().template begin<0>();
        const ElementIterator eEndIt = this->gridView_().template end<0>();
        for (; eIt != eEndIt; ++eIt) {
            this->localResidual().evalPhaseStorage(*eIt, phaseIdx);

            for (int i = 0; i < eIt->template count<dim>(); ++i)
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
     * \brief Returns the relative weight of a primary variable for
     *        calculating relative errors.
     *
     * \param globalIdx The global vertex index
     * \param pvIdx The primary variable index
     */
    Scalar primaryVarWeight(const int globalIdx, const int pvIdx) const
    {
        if (Indices::pressureIdx == pvIdx)
            return std::min(1.0/this->prevSol()[globalIdx][pvIdx], 1.0);
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
     * \param globalIdx The global vertex index
     * \param oldSol Evaluate function with solution of current or previous time step
     */
    int phasePresence(int globalIdx, bool oldSol) const
    {
        return oldSol ? staticVertexDat_[globalIdx].oldPhasePresence
            : staticVertexDat_[globalIdx].phasePresence;
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

        ScalarField *sN    = writer.allocateManagedBuffer(numVertices);
        ScalarField *sW    = writer.allocateManagedBuffer(numVertices);
        ScalarField *pN    = writer.allocateManagedBuffer(numVertices);
        ScalarField *pW    = writer.allocateManagedBuffer(numVertices);
        ScalarField *pC    = writer.allocateManagedBuffer(numVertices);
        ScalarField *rhoW  = writer.allocateManagedBuffer(numVertices);
        ScalarField *rhoN  = writer.allocateManagedBuffer(numVertices);
        ScalarField *mobW  = writer.allocateManagedBuffer(numVertices);
        ScalarField *mobN = writer.allocateManagedBuffer(numVertices);
        ScalarField *phasePresence = writer.allocateManagedBuffer(numVertices);
        ScalarField *massFrac[numPhases][numComponents];
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                massFrac[phaseIdx][compIdx] = writer.allocateManagedBuffer(numVertices);
        ScalarField *temperature = writer.allocateManagedBuffer(numVertices);
        ScalarField *poro = writer.allocateManagedBuffer(numVertices);
        ScalarField *cellNum =writer.allocateManagedBuffer (numVertices);
        VectorField *velocityN = writer.template allocateManagedBuffer<double, dim>(numVertices);
        VectorField *velocityW = writer.template allocateManagedBuffer<double, dim>(numVertices);

        if(velocityOutput_) // check if velocity output is demanded
        {
            // initialize velocity fields
            for (unsigned int i = 0; i < numVertices; ++i)
            {
                (*velocityN)[i] = Scalar(0);
                (*velocityW)[i] = Scalar(0);
                (*cellNum)[i] = Scalar(0.0);
            }
        }

        unsigned numElements = this->gridView_().size(0);
        ScalarField *rank =
            writer.allocateManagedBuffer (numElements);

        FVElementGeometry fvGeometry;
        VolumeVariables volVars;

        ElementIterator eIt = this->gridView_().template begin<0>();
        ElementIterator eEndIt = this->gridView_().template end<0>();
        for (; eIt != eEndIt; ++eIt)
        {
            int idx = this->elementMapper().map(*eIt);
            (*rank)[idx] = this->gridView_().comm().rank();
            fvGeometry.update(this->gridView_(), *eIt);

            int numVerts = eIt->template count<dim> ();
            for (int scvIdx = 0; scvIdx < numVerts; ++scvIdx)
            {
                int globalIdx = this->vertexMapper().map(*eIt, scvIdx, dim);
                volVars.update(sol[globalIdx],
                               this->problem_(),
                               *eIt,
                               fvGeometry,
                               scvIdx,
                               false);
                (*sN)[globalIdx]    = volVars.saturation(nPhaseIdx);
                (*sW)[globalIdx]    = volVars.saturation(wPhaseIdx);
                (*pN)[globalIdx]    = volVars.pressure(nPhaseIdx);
                (*pW)[globalIdx]    = volVars.pressure(wPhaseIdx);
                (*pC)[globalIdx]    = volVars.capillaryPressure();
                (*rhoW)[globalIdx]  = volVars.fluidState().density(wPhaseIdx);
                (*rhoN)[globalIdx]  = volVars.fluidState().density(nPhaseIdx);
                (*mobW)[globalIdx]  = volVars.mobility(wPhaseIdx);
                (*mobN)[globalIdx]  = volVars.mobility(nPhaseIdx);
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
            }

            if(velocityOutput_)
            {
                // calculate vertex velocities
                GlobalPosition tmpVelocity[numPhases];

                for(int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                {
                    tmpVelocity[phaseIdx]  = Scalar(0.0);
                }

                typedef Dune::BlockVector<Dune::FieldVector<Scalar, dim> > SCVVelocities;
                SCVVelocities scvVelocityW(8), scvVelocityN(8);

                scvVelocityW = 0;
                scvVelocityN = 0;

                ElementVolumeVariables elemVolVars;

                elemVolVars.update(this->problem_(),
                                   *eIt,
                                   fvGeometry,
                                   false /* oldSol? */);

                for (int faceIdx = 0; faceIdx < fvGeometry.numEdges; faceIdx++)
                {

                    FluxVariables fluxVars(this->problem_(),
                                           *eIt,
                                           fvGeometry,
                                           faceIdx,
                                           elemVolVars);

                    for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                    {

                        // data attached to upstream and the downstream vertices
                        // of the current phase
                        const VolumeVariables up =
                            elemVolVars[fluxVars.upstreamIdx(phaseIdx)];
                        const VolumeVariables dn =
                            elemVolVars[fluxVars.downstreamIdx(phaseIdx)];

                        // local position of integration point
                        const Dune::FieldVector<Scalar, dim>& localPosIP = fvGeometry.subContVolFace[faceIdx].ipLocal;

                        // Transformation of the global normal vector to normal vector in the reference element
                        const Dune::FieldMatrix<CoordScalar, dim, dim> jacobianT1 = eIt->geometry().jacobianTransposed(localPosIP);
                        const GlobalPosition globalNormal = fluxVars.face().normal;

                        GlobalPosition localNormal(0);
                        jacobianT1.mv(globalNormal, localNormal);
                        // note only works for cubes
                        const Scalar localArea = pow(2,-(dim-1));

                        localNormal /= localNormal.two_norm();

                        // Get the Darcy velocities. The Darcy velocities are divided by the area of the subcontrolvolume
                        // face in the reference element.
                        massUpwindWeight_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Implicit, MassUpwindWeight);
                        PhasesVector flux;
                        flux[phaseIdx] = fluxVars.KmvpNormal(phaseIdx)
                            * (massUpwindWeight_
                               * up.mobility(phaseIdx)
                               + (1- massUpwindWeight_)
                               * dn.mobility(phaseIdx)) / localArea;

                        // transform the normal Darcy velocity into a vector
                        tmpVelocity[phaseIdx] = localNormal;
                        tmpVelocity[phaseIdx] *= flux[phaseIdx];

                        if(phaseIdx == wPhaseIdx){
                            scvVelocityW[fluxVars.face().i] += tmpVelocity[phaseIdx];
                            scvVelocityW[fluxVars.face().j] += tmpVelocity[phaseIdx];
                        }
                        else if(phaseIdx == nPhaseIdx){
                            scvVelocityN[fluxVars.face().i] += tmpVelocity[phaseIdx];
                            scvVelocityN[fluxVars.face().j] += tmpVelocity[phaseIdx];
                        }
                    }
                }

                typedef Dune::GenericReferenceElements<Scalar, dim> ReferenceElements;
                const Dune::FieldVector<Scalar, dim>& localPos =
                    ReferenceElements::general(eIt->geometry().type()).position(0, 0);

                // get the transposed Jacobian of the element mapping
                const Dune::FieldMatrix<CoordScalar, dim, dim> jacobianT2 = eIt->geometry().jacobianTransposed(localPos);

                // transform vertex velocities from local to global coordinates
                for (int scvIdx = 0; scvIdx < numVerts; ++scvIdx)
                {
                    int globalIdx = this->vertexMapper().map(*eIt, scvIdx, dim);
                    // calculate the subcontrolvolume velocity by the Piola transformation
                    Dune::FieldVector<CoordScalar, dim> scvVelocity(0);

                    jacobianT2.mtv(scvVelocityW[scvIdx], scvVelocity);
                    scvVelocity /= eIt->geometry().integrationElement(localPos);
                    // add up the wetting phase subcontrolvolume velocities for each vertex
                    (*velocityW)[globalIdx] += scvVelocity;

                    jacobianT2.mtv(scvVelocityN[scvIdx], scvVelocity);
                    scvVelocity /= eIt->geometry().integrationElement(localPos);
                    // add up the nonwetting phase subcontrolvolume velocities for each vertex
                    (*velocityN)[globalIdx] += scvVelocity;
                }
            }
        }
        if(velocityOutput_)
        {
            // divide the vertex velocities by the number of adjacent scvs i.e. cells
            for(unsigned int globalIdx = 0; globalIdx<numVertices; ++globalIdx){
                (*velocityW)[globalIdx] /= (*cellNum)[globalIdx];
                (*velocityN)[globalIdx] /= (*cellNum)[globalIdx];
            }
        }


        writer.attachVertexData(*sN,     "Sn");
        writer.attachVertexData(*sW,     "Sw");
        writer.attachVertexData(*pN,     "pN");
        writer.attachVertexData(*pW,     "pW");
        writer.attachVertexData(*pC,     "pC");
        writer.attachVertexData(*rhoW,   "rhoW");
        writer.attachVertexData(*rhoN,   "rhoN");
        writer.attachVertexData(*mobW,   "mobW");
        writer.attachVertexData(*mobN,   "mobN");
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                std::ostringstream oss;
                oss << "X_" << FluidSystem::phaseName(phaseIdx) << "^" << FluidSystem::componentName(compIdx);
                writer.attachVertexData(*massFrac[phaseIdx][compIdx], oss.str());
            }
        }
        writer.attachVertexData(*poro, "porosity");
        writer.attachVertexData(*temperature,    "temperature");
        writer.attachVertexData(*phasePresence,  "phase presence");

        if(velocityOutput_) // check if velocity output is demanded
        {
            writer.attachVertexData(*velocityW,  "velocityW", dim);
            writer.attachVertexData(*velocityN,  "velocityN", dim);
        }
        writer.attachCellData(*rank, "process rank");
    }

    /*!
     * \brief Write the current solution to a restart file.
     *
     * \param outStream The output stream of one vertex for the restart file
     * \param vertex The vertex
     */
    void serializeEntity(std::ostream &outStream, const Vertex &vertex)
    {
        // write primary variables
        ParentType::serializeEntity(outStream, vertex);

        int globalIdx = this->dofMapper().map(vertex);
        if (!outStream.good())
            DUNE_THROW(Dune::IOError, "Could not serialize vertex " << globalIdx);

        outStream << staticVertexDat_[globalIdx].phasePresence << " ";
    }

    /*!
     * \brief Reads the current solution for a vertex from a restart
     *        file.
     *
     * \param inStream The input stream of one vertex from the restart file
     * \param vertex The vertex
     */
    void deserializeEntity(std::istream &inStream, const Vertex &vertex)
    {
        // read primary variables
        ParentType::deserializeEntity(inStream, vertex);

        // read phase presence
        int globalIdx = this->dofMapper().map(vertex);
        if (!inStream.good())
            DUNE_THROW(Dune::IOError,
                       "Could not deserialize vertex " << globalIdx);

        inStream >> staticVertexDat_[globalIdx].phasePresence;
        staticVertexDat_[globalIdx].oldPhasePresence
            = staticVertexDat_[globalIdx].phasePresence;

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

        FVElementGeometry fvGeometry;
        static VolumeVariables volVars;
        ElementIterator eIt = this->gridView_().template begin<0> ();
        const ElementIterator &eEndIt = this->gridView_().template end<0> ();
        for (; eIt != eEndIt; ++eIt)
        {
            fvGeometry.update(this->gridView_(), *eIt);
            for (int scvIdx = 0; scvIdx < fvGeometry.numSCV; ++scvIdx)
            {
                int globalIdx = this->vertexMapper().map(*eIt, scvIdx, dim);

                if (staticVertexDat_[globalIdx].visited)
                    continue;

                staticVertexDat_[globalIdx].visited = true;
                volVars.update(curGlobalSol[globalIdx],
                               this->problem_(),
                               *eIt,
                               fvGeometry,
                               scvIdx,
                               false);
                const GlobalPosition &globalPos = eIt->geometry().corner(scvIdx);
                if (primaryVarSwitch_(curGlobalSol,
                                      volVars,
                                      globalIdx,
                                      globalPos))
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
        for (int idx = 0; idx < numVertices; ++idx)
        {
            staticVertexDat_[idx].phasePresence
                = staticVertexDat_[idx].oldPhasePresence;
            staticVertexDat_[idx].wasSwitched = false;
        }
    }

    /*!
     * \brief Set the old phase of all verts state to the current one.
     */
    void updateOldPhasePresence_()
    {
        int numVertices = this->gridView_().size(dim);
        for (int idx = 0; idx < numVertices; ++idx)
        {
            staticVertexDat_[idx].oldPhasePresence
                = staticVertexDat_[idx].phasePresence;
            staticVertexDat_[idx].wasSwitched = false;
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
        if (phasePresence == nPhaseOnly)
        {
            // calculate mole fraction in the hypothetic wetting phase
            Scalar xww = volVars.fluidState().moleFraction(wPhaseIdx, wCompIdx);
            Scalar xwn = volVars.fluidState().moleFraction(wPhaseIdx, nCompIdx);

            Scalar xwMax = 1.0;
            if (xww + xwn > xwMax)
                wouldSwitch = true;
            if (staticVertexDat_[globalIdx].wasSwitched)
                xwMax *= 1.02;

            // if the sum of the mole fractions would be larger than
            // 100%, wetting phase appears
            if (xww + xwn > xwMax)
            {
                // wetting phase appears
                std::cout << "wetting phase appears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", xww + xwn: "
                          << xww + xwn << std::endl;
                newPhasePresence = bothPhases;
                if (formulation == pnSw)
                    globalSol[globalIdx][switchIdx] = 0.0;
                else if (formulation == pwSn)
                    globalSol[globalIdx][switchIdx] = 1.0;
            }
        }
        else if (phasePresence == wPhaseOnly)
        {
            // calculate fractions of the partial pressures in the
            // hypothetic nonwetting phase
            Scalar xnw = volVars.fluidState().moleFraction(nPhaseIdx, wCompIdx);
            Scalar xnn = volVars.fluidState().moleFraction(nPhaseIdx, nCompIdx);

            Scalar xgMax = 1.0;
            if (xnw + xnn > xgMax)
                wouldSwitch = true;
            if (staticVertexDat_[globalIdx].wasSwitched)
                xgMax *= 1.02;

            // if the sum of the mole fractions would be larger than
            // 100%, nonwetting phase appears
            if (xnw + xnn > xgMax)
            {
                // nonwetting phase appears
                std::cout << "nonwetting phase appears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", xnw + xnn: "
                          << xnw + xnn << std::endl;
                newPhasePresence = bothPhases;
                if (formulation == pnSw)
                    globalSol[globalIdx][switchIdx] = 0.999;
                else if (formulation == pwSn)
                    globalSol[globalIdx][switchIdx] = 0.001;
            }
        }
        else if (phasePresence == bothPhases)
        {
            Scalar Smin = 0.0;
            if (staticVertexDat_[globalIdx].wasSwitched)
                Smin = -0.01;

            if (volVars.saturation(nPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // nonwetting phase disappears
                std::cout << "Nonwetting phase disappears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", Sn: "
                          << volVars.saturation(nPhaseIdx) << std::endl;
                newPhasePresence = wPhaseOnly;

                globalSol[globalIdx][switchIdx]
                    = volVars.fluidState().massFraction(wPhaseIdx, nCompIdx);
            }
            else if (volVars.saturation(wPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // wetting phase disappears
                std::cout << "Wetting phase disappears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << ", Sw: "
                          << volVars.saturation(wPhaseIdx) << std::endl;
                newPhasePresence = nPhaseOnly;

                globalSol[globalIdx][switchIdx]
                    = volVars.fluidState().massFraction(nPhaseIdx, wCompIdx);
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
