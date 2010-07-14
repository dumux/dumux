// $Id: 2p2cboxjacobian.hh 3795 2010-06-25 16:08:04Z melanie $
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
 *   Copyright (C) 2009-2010 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_NEW_2P2C_BOX_JACOBIAN_BASE_HH
#define DUMUX_NEW_2P2C_BOX_JACOBIAN_BASE_HH

#include <dumux/boxmodels/boxscheme/boxscheme.hh>
#include <dumux/common/math.hh>

#include "2p2cproperties.hh"

#include "2p2cvertexdata.hh"
#include "2p2celementdata.hh"
#include "2p2cfluxdata.hh"

#include "2p2cnewtoncontroller.hh"

#include <iostream>
#include <vector>

//#define VELOCITY_OUTPUT 1 // uncomment this line if an output of the velocity is needed

namespace Dumux
{
/*!
 * \ingroup TwoPTwoCBoxModel
 * \brief 2P-2C specific details needed to approximately calculate
 *        the local jacobian in the BOX scheme.
 *
 * This class is used to fill the gaps in BoxJacobian for the 2P-2C flow.
 */
template<class TypeTag>
class TwoPTwoCBoxJacobian: public BoxJacobian<TypeTag>
{
protected:
    typedef TwoPTwoCBoxJacobian<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalJacobian)) Implementation;
    typedef BoxJacobian<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GridView::Grid::ctype CoordScalar;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::SolutionVector SolutionVector;
    typedef typename SolutionTypes::SolutionOnElement SolutionOnElement;
    typedef typename SolutionTypes::PrimaryVarVector PrimaryVarVector;
    typedef typename SolutionTypes::BoundaryTypeVector BoundaryTypeVector;
    typedef typename SolutionTypes::JacobianAssembler JacobianAssembler;

    typedef TwoPTwoCFluidState<TypeTag> FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices)) Indices;

    enum
    {
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

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

    typedef typename GridView::CollectiveCommunication CollectiveCommunication;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexData)) VertexData;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementData)) ElementData;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxData)) FluxData;

    typedef std::vector<VertexData> VertexDataArray;

    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> Tensor;

    static const Scalar mobilityUpwindAlpha =
            GET_PROP_VALUE(TypeTag, PTAG(MobilityUpwindAlpha));

    /*!
     * \brief Data which is attached to each vertex and is not only
     *        stored locally.
     */
    struct StaticVertexData
    {
        int phasePresence;
        bool wasSwitched;

        int oldPhasePresence;
        bool visited;
    };

public:
    TwoPTwoCBoxJacobian(Problem &problem) :
        ParentType(problem), staticVertexDat_(this->gridView_.size(dim))
    {
        switchFlag_ = false;
    }
    ;

    /*!
     * \brief Evaluate the amount all conservation quantites
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     */
    void computeStorage(PrimaryVarVector &result, int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const VertexDataArray &elemDat = usePrevSol ? this->prevElemDat_
                : this->curElemDat_;
        const VertexData &vertDat = elemDat[scvIdx];

        // compute storage term of all components within all phases
        result = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                int eqIdx = (compIdx == lCompIdx) ? contiLEqIdx : contiGEqIdx;
                result[eqIdx] += vertDat.density(phaseIdx)
                        * vertDat.saturation(phaseIdx)
                        * vertDat.fluidState().massFrac(phaseIdx, compIdx);
            }
        }
        result *= vertDat.porosity();
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a subcontrol volume.
     */
    void computeFlux(PrimaryVarVector &flux, int faceIdx) const
    {
        FluxData vars(this->problem_, this->curElement_(),
                this->curElementGeom_, faceIdx, this->curElemDat_);

        flux = 0;
        asImp_()->computeAdvectiveFlux(flux, vars);
        asImp_()->computeDiffusiveFlux(flux, vars);

        // the face normal points into the outward direction, so we
        // have to multiply all fluxes with -1
        flux *= -1;
    }

    /*!
     * \brief Evaluates the advective mass flux of all components over
     *        a face of a subcontrol volume.
     */
    void computeAdvectiveFlux(PrimaryVarVector &flux, const FluxData &vars) const
    {
        ////////
        // advective fluxes of all components in all phases
        ////////
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // data attached to upstream and the downstream vertices
            // of the current phase
            const VertexData &up =
                    this->curElemDat_[vars.upstreamIdx(phaseIdx)];
            const VertexData &dn = this->curElemDat_[vars.downstreamIdx(
                    phaseIdx)];

            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            {
                int eqIdx = (compIdx == lCompIdx) ? contiLEqIdx : contiGEqIdx;
                // add advective flux of current component in current
                // phase
                if (mobilityUpwindAlpha > 0.0)
                    // upstream vertex
                    flux[eqIdx] += vars.KmvpNormal(phaseIdx)
                            * mobilityUpwindAlpha * (up.density(phaseIdx)
                            * up.mobility(phaseIdx) * up.fluidState().massFrac(
                            phaseIdx, compIdx));
                if (mobilityUpwindAlpha < 1.0)
                    // downstream vertex
                    flux[eqIdx] += vars.KmvpNormal(phaseIdx) * (1
                            - mobilityUpwindAlpha) * (dn.density(phaseIdx)
                            * dn.mobility(phaseIdx) * dn.fluidState().massFrac(
                            phaseIdx, compIdx));
            }
        }
    }

    /*!
     * \brief Adds the diffusive mass flux of all components over
     *        a face of a subcontrol volume.
     */
    void computeDiffusiveFlux(PrimaryVarVector &flux, const FluxData &vars) const
    {
        // add diffusive flux of gas component in liquid phase
        Scalar tmp = -vars.porousDiffCoeff(lPhaseIdx) * vars.molarDensityAtIP(
                lPhaseIdx) * (vars.molarConcGrad(lPhaseIdx)
                * vars.face().normal);
        flux[contiGEqIdx] += tmp * FluidSystem::molarMass(gCompIdx);
        flux[contiLEqIdx] -= tmp * FluidSystem::molarMass(lCompIdx);
        ;

        // add diffusive flux of liquid component in gas phase
        tmp = -vars.porousDiffCoeff(gPhaseIdx) * vars.molarDensityAtIP(
                gPhaseIdx) * (vars.molarConcGrad(gPhaseIdx)
                * vars.face().normal);
        flux[contiLEqIdx] += tmp * FluidSystem::molarMass(lCompIdx);
        flux[contiGEqIdx] -= tmp * FluidSystem::molarMass(gCompIdx);
    }

    /*!
     * \brief Calculate the source term of the equation
     */
    void computeSource(PrimaryVarVector &q, int localVertexIdx)
    {
        this->problem_.source(q, this->curElement_(), this->curElementGeom_,
                localVertexIdx);
    }

    /*!
     * \brief Initialize the static data with the initial solution.
     *
     * Called by TwoPTwoCBoxModel::initial()
     */
    void initStaticData()
    {
        setSwitched(false);

        VertexIterator it = this->gridView_.template begin<dim> ();
        const VertexIterator &endit = this->gridView_.template end<dim> ();
        for (; it != endit; ++it)
        {
            int globalIdx = this->problem_.model().dofEntityMapper().map(*it);
            const GlobalPosition &globalPos = it->geometry().corner(0);

            // initialize phase presence
            staticVertexDat_[globalIdx].phasePresence
                    = this->problem_.initialPhasePresence(*it, globalIdx,
                            globalPos);
            staticVertexDat_[globalIdx].wasSwitched = false;

            staticVertexDat_[globalIdx].oldPhasePresence
                    = staticVertexDat_[globalIdx].phasePresence;
        }
    }

    /*!
     * \brief Update the static data of all vertices in the grid.
     */
    void updateStaticData(SolutionVector &curGlobalSol,
            SolutionVector &oldGlobalSol)
    {

        bool wasSwitched = false;

        for (unsigned i = 0; i < staticVertexDat_.size(); ++i)
            staticVertexDat_[i].visited = false;

        static VertexData vertexData;
        ElementIterator it = this->gridView_.template begin<0> ();
        const ElementIterator &endit = this->gridView_.template end<0> ();
        for (; it != endit; ++it)
        {
            this->curElementGeom_.update(*it);
            for (int i = 0; i < this->curElementGeom_.numVertices; ++i)
            {
                int globalIdx = this->vertexMapper().map(*it, i, dim);

                if (staticVertexDat_[globalIdx].visited)
                    continue;

                staticVertexDat_[globalIdx].visited = true;
                vertexData.update(curGlobalSol[globalIdx], *it,
                        this->curElementGeom_, i, this->problem(), false);
                const GlobalPosition &global = it->geometry().corner(i);
                wasSwitched = primaryVarSwitch_(curGlobalSol, vertexData,
                        globalIdx, global) || wasSwitched;
            }
        }

        // make sure that if there was a variable switch in an
        // other partition we will also set the switch flag
        // for our partition.
        wasSwitched = this->gridView_.comm().max(wasSwitched);

        setSwitched(wasSwitched);
    }

    /*!
     * \brief Set the old phase of all verts state to the current one.
     */
    void updateOldPhasePresence()
    {
        int numVertices = this->gridView_.size(dim);
        for (int i = 0; i < numVertices; ++i)
        {
            staticVertexDat_[i].oldPhasePresence
                    = staticVertexDat_[i].phasePresence;
            staticVertexDat_[i].wasSwitched = false;
        }
    }

    /*!
     * \brief Returns the phase presence of the current or the old solution of a vertex.
     */
    int phasePresence(int globalVertexIdx, bool oldSol) const
    {
        return oldSol ? staticVertexDat_[globalVertexIdx].oldPhasePresence
                : staticVertexDat_[globalVertexIdx].phasePresence;
    }

    /*!
     * \brief Reset the current phase presence of all vertices to the old one.
     *
     * This is done after an update failed.
     */
    void resetPhasePresence()
    {
        int numVertices = this->gridView_.size(dim);
        for (int i = 0; i < numVertices; ++i)
        {
            staticVertexDat_[i].phasePresence
                    = staticVertexDat_[i].oldPhasePresence;
            staticVertexDat_[i].wasSwitched = false;
        }
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
     * \brief Set whether there was a primary variable switch after in the last
     *        timestep.
     */
    void setSwitched(bool yesno)
    {
        switchFlag_ = yesno;
    }

    /*!
     * \brief Calculate mass of both components in the whole model domain
     *         and get minimum and maximum values of primary variables
     *
     */
    void calculateMass(const SolutionVector &globalSol, Dune::FieldVector<
            Scalar, 2> &massGas, Dune::FieldVector<Scalar, 2> &massLiquid)
    {
        massGas = 0;
        massLiquid = 0;

        ElementIterator elementIt =
                this->model().gridView().template begin<0> ();
        ElementIterator endit = this->model().gridView().template end<0> ();

        SolutionOnElement curSol;
        VertexDataArray elemDat;

        Scalar minSat = 1e100;
        Scalar maxSat = -1e100;
        Scalar minP = 1e100;
        Scalar maxP = -1e100;
        Scalar minTe = 1e100;
        Scalar maxTe = -1e100;
        Scalar minX = 1e100;
        Scalar maxX = -1e100;

        // Loop over elements
        for (; elementIt != endit; ++elementIt)
        {
            if (elementIt->partitionType() != Dune::InteriorEntity)
                continue;

            setCurrentElement(*elementIt);

            int numLocalVerts = elementIt->template count<dim> ();
            curSol.resize(numLocalVerts);
            elemDat.resize(numLocalVerts);
            this->restrictToElement(curSol, globalSol);
            this->updateElementData_(elemDat, curSol, false);
            // get geometry type


            // Loop over element vertices
            for (int i = 0; i < numLocalVerts; ++i)
            {

                const VertexData &vdat = elemDat[i];
                const FluidState &fs = vdat.fluidState();
                Scalar vol = this->curElementGeom_.subContVol[i].volume;

                Scalar satN = fs.saturation(gPhaseIdx);
                Scalar xAW = fs.massFrac(lPhaseIdx, gCompIdx);
                Scalar pW = fs.phasePressure(lPhaseIdx);
                Scalar T = fs.temperature();

                // get minimum and maximum values of primary variables
                minSat = std::min(minSat, satN);
                maxSat = std::max(maxSat, satN);
                minP = std::min(minP, pW);
                maxP = std::max(maxP, pW);
                minX = std::min(minX, xAW);
                maxX = std::max(maxX, xAW);
                minTe = std::min(minTe, T);
                maxTe = std::max(maxTe, T);

                // calculate total mass
                Scalar mGas = vdat.porosity() * fs.saturation(gPhaseIdx)
                        * fs.density(gPhaseIdx) * vol;
                massGas[lCompIdx] += mGas * fs.massFrac(gPhaseIdx, lCompIdx);
                massGas[gCompIdx] += mGas * fs.massFrac(gPhaseIdx, gCompIdx);

                Scalar mLiquid = vdat.porosity() * fs.saturation(lPhaseIdx)
                        * fs.density(lPhaseIdx) * vol;
                massLiquid[lCompIdx] += mLiquid * fs.massFrac(lPhaseIdx,
                        lCompIdx);
                massLiquid[gCompIdx] += mLiquid * fs.massFrac(lPhaseIdx,
                        gCompIdx);
            }
        }

        massGas[lCompIdx] = this->gridView_.comm().sum(massGas[lCompIdx]);
        massGas[gCompIdx] = this->gridView_.comm().sum(massGas[gCompIdx]);
        massLiquid[lCompIdx] = this->gridView_.comm().sum(massLiquid[lCompIdx]);
        massLiquid[gCompIdx] = this->gridView_.comm().sum(massLiquid[gCompIdx]);

        Scalar minS = this->gridView_.comm().min(minSat);
        Scalar maxS = this->gridView_.comm().max(maxSat);
        Scalar minPr = this->gridView_.comm().min(minP);
        Scalar maxPr = this->gridView_.comm().max(maxP);
        Scalar minXw = this->gridView_.comm().min(minX);
        Scalar maxXw = this->gridView_.comm().max(maxX);
        Scalar minT = this->gridView_.comm().min(minTe);
        Scalar maxT = this->gridView_.comm().max(maxTe);

        // IF PARALLEL: mass calculation still needs to be adjusted
        //        massGas = this->gridView_.comm().sum(massGas);
        //        massLiquid = this->gridView_.comm().sum(massLiquid);

        if (this->gridView_.comm().rank() == 0) // IF PARALLEL: only print by processor with rank() == 0
        {
            // print minimum and maximum values
            std::cout << "gas phase saturation: min = " << minS << ", max = "
                    << maxS << std::endl;
            std::cout << "liquid phase pressure: min = " << minPr << ", max = "
                    << maxPr << std::endl;
            std::cout << "mass fraction gCompIdx: min = " << minXw
                    << ", max = " << maxXw << std::endl;
            std::cout << "temperature: min = " << minT << ", max = " << maxT
                    << std::endl;
        }
    }

    /*!
     * \brief Add the mass fraction of air in water to VTK output of
     *        the current timestep.
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer,
            const SolutionVector &globalSol)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->gridView_.size(dim);

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
        ScalarField *temperature = writer.template createField<Scalar, 1> (
                numVertices);
        ScalarField *phasePresence = writer.template createField<Scalar, 1> (
                numVertices);
        ScalarField *massFrac[numPhases][numComponents];
        for (int i = 0; i < numPhases; ++i)
            for (int j = 0; j < numComponents; ++j)
                massFrac[i][j] = writer.template createField<Scalar, 1> (
                        numVertices);

#ifdef VELOCITY_OUTPUT  // check if velocity output is demanded
        ScalarField *velocityX = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *velocityY = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *velocityZ = writer.template createField<Scalar, 1>(numVertices);
        Scalar maxV=0.; // variable to store the maximum face velocity

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

        unsigned numElements = this->gridView_.size(0);
        ScalarField *rank =
                writer.template createField<Scalar, 1> (numElements);

        SolutionOnElement tmpSol;
        VertexDataArray elemDat;

        ElementIterator elementIt = this->gridView_.template begin<0> ();
        ElementIterator endit = this->gridView_.template end<0> ();

        for (; elementIt != endit; ++elementIt)
        {
            int idx = this->problem_.model().elementMapper().map(*elementIt);
            (*rank)[idx] = this->gridView_.comm().rank();

            int numLocalVerts = elementIt->template count<dim> ();
            tmpSol.resize(numLocalVerts);

            setCurrentElement(*elementIt);
            this->restrictToElement(tmpSol, globalSol);
            updateElementData_(elemDat, tmpSol, false);

            for (int i = 0; i < numLocalVerts; ++i)
            {
                int globalIdx = this->problem_.model().vertexMapper().map(
                        *elementIt, i, dim);

                (*Sg)[globalIdx] = elemDat[i].saturation(gPhaseIdx);
                (*Sl)[globalIdx] = elemDat[i].saturation(lPhaseIdx);
                (*pg)[globalIdx] = elemDat[i].pressure(gPhaseIdx);
                (*pl)[globalIdx] = elemDat[i].pressure(lPhaseIdx);
                (*pc)[globalIdx] = elemDat[i].capillaryPressure();
                (*rhoL)[globalIdx] = elemDat[i].fluidState().density(lPhaseIdx);
                (*rhoG)[globalIdx] = elemDat[i].fluidState().density(gPhaseIdx);
                ;
                (*mobL)[globalIdx] = elemDat[i].mobility(lPhaseIdx);
                (*mobG)[globalIdx] = elemDat[i].mobility(gPhaseIdx);
                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                    for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                    {
                        (*massFrac[phaseIdx][compIdx])[globalIdx]
                                = elemDat[i].fluidState().massFrac(phaseIdx,
                                        compIdx);

                        Valgrind::CheckDefined(
                                (*massFrac[phaseIdx][compIdx])[globalIdx][0]);
                    }
                (*temperature)[globalIdx] = elemDat[i].temperature();
                (*phasePresence)[globalIdx]
                        = staticVertexDat_[globalIdx].phasePresence;
            };

#ifdef VELOCITY_OUTPUT        // check if velocity output is demanded
            // In the box method, the velocity is evaluated on the FE-Grid. However, to get an
            // average apparent velocity at the vertex, all contributing velocities have to be interpolated.
            GlobalPosition velocity(0.);
            // loop over the phases
            for (int faceIdx = 0; faceIdx< this->curElementGeom_.numEdges; faceIdx++)
            {
                //prepare the flux calculations (set up and prepare geometry, FE gradients)
                FluxData fluxDat(this->problem_,
                        this->curElement_(),
                        this->curElementGeom_,
                        faceIdx,
                        elemDat);

                // choose phase of interest. Alternatively, a loop over all phases would be possible.
                int phaseIdx = gPhaseIdx;

                // get darcy velocity
                velocity = fluxDat.KmvpNormal(phaseIdx); // mind the sign: vDarcy = kf grad p

                // up+downstream mobility
                const VertexData &up = this->curElemDat_[fluxDat.upstreamIdx(phaseIdx)];
                const VertexData &down = this->curElemDat_[fluxDat.downstreamIdx(phaseIdx)];
                Scalar scvfArea = fluxDat.face().normal.two_norm(); //get surface area to weight velocity at the IP with the surface area
                velocity *= (mobilityUpwindAlpha*up.mobility(phaseIdx) + (1-mobilityUpwindAlpha)*down.mobility(phaseIdx))* scvfArea;

                int vertIIdx = this->problem().model().vertexMapper().map(this->curElement_(),
                        fluxDat.face().i,
                        dim);
                int vertJIdx = this->problem().model().vertexMapper().map(this->curElement_(),
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

#ifdef VELOCITY_OUTPUT        // check if velocity output is demanded
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
        writer.addVertexData(temperature, "temperature");
        writer.addVertexData(phasePresence, "phase presence");

#ifdef VELOCITY_OUTPUT        // check if velocity output is demanded
        writer.addVertexData(velocityX, "Vx");
        if (dim >= 2)
        writer.addVertexData(velocityY, "Vy");
        if (dim == 3)
        writer.addVertexData(velocityZ, "Vz");
#endif
        writer.addCellData(rank, "process rank");
    }

    /*!
     * \brief Add the vector fields for analysing the convergence of
     *        the newton method to the a VTK multi writer.
     *
     * \param writer  The VTK multi writer where the fields should be added.
     * \param oldSol  The solution function before the Newton update
     * \param update  The delte of the solution function before and after the Newton update
     */
    template<class MultiWriter>
    void addConvergenceVtkFields(MultiWriter &writer,
            const SolutionVector &oldSol, const SolutionVector &update)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        SolutionVector globalDef(this->model(), 0.0);
        this->model().globalResidual(oldSol, globalDef);

        // create the required scalar fields
        unsigned numVertices = this->gridView_.size(dim);
        //unsigned numElements = this->gridView_.size(0);

        // global defect of the two auxiliary equations
        ScalarField* gd[numEq];
        ScalarField* delta[numEq];
        ScalarField* x[numEq];
        for (int i = 0; i < numEq; ++i)
        {
            x[i] = writer.template createField<Scalar, 1> (numVertices);
            delta[i] = writer.template createField<Scalar, 1> (numVertices);
            gd[i] = writer.template createField<Scalar, 1> (numVertices);
        }

        ElementIterator eIt = this->gridView_.template begin<0> ();
        ElementIterator eEndIt = this->gridView_.template end<0> ();

        for (; eIt != eEndIt; ++eIt)
        {
            this->curElementGeom_.update(*eIt);
            for (int scvIdx = 0; scvIdx < this->curElementGeom_.numVertices; ++scvIdx)
            {
                int globalIdx = this->problem().model().vertexMapper().map(
                        *eIt, scvIdx, dim);
                for (int i = 0; i < numEq; ++i)
                {
                    (*x[i])[globalIdx] = oldSol[globalIdx][i];
                    (*delta[i])[globalIdx] = -update[globalIdx][i];
                    (*gd[i])[globalIdx] = globalDef[globalIdx][i];
                }
            }
        }

        for (int i = 0; i < numEq; ++i)
        {
            writer.addVertexData(x[i],
                    (boost::format("x_%i") % i).str().c_str());
            writer.addVertexData(delta[i],
                    (boost::format("delta_%i") % i).str().c_str());
            writer.addVertexData(gd[i],
                    (boost::format("defect_%i") % i).str().c_str());
        }

        addOutputVtkFields(writer, oldSol);
    }

    /*!
     * \brief Reads the current solution for a vertex from a restart
     *        file.
     */
    void deserializeEntity(std::istream &inStream, const Vertex &vert)
    {
        int vertIdx = this->problem_.model().dofEntityMapper().map(vert);

        // read phase presence
        if (!inStream.good())
        {
            DUNE_THROW(Dune::IOError, "Could not deserialize vertex "
                    << vertIdx);
        }

        inStream >> staticVertexDat_[vertIdx].phasePresence;
        staticVertexDat_[vertIdx].oldPhasePresence
                = staticVertexDat_[vertIdx].phasePresence;
    }
    ;

    /*!
     * \brief Write the current phase presence of an vertex to a restart
     *        file.
     */
    void serializeEntity(std::ostream &outStream, const Vertex &vert)
    {
        int vertIdx = this->problem_.model().dofEntityMapper().map(vert);

        if (!outStream.good())
        {
            DUNE_THROW(Dune::IOError, "Could not serialize vertex " << vertIdx);
        }

        outStream << staticVertexDat_[vertIdx].phasePresence << " ";
    }
    ;

protected:
    Implementation *asImp_()
    {
        return static_cast<Implementation *> (this);
    }
    const Implementation *asImp_() const
    {
        return static_cast<const Implementation *> (this);
    }

    //  perform variable switch at a vertex; Returns true if a
    //  variable switch was performed.
    bool primaryVarSwitch_(SolutionVector &globalSol,
            const VertexData &vertexData, int globalIdx,
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
            Scalar xll = vertexData.fluidState().moleFrac(lPhaseIdx, lCompIdx);
            Scalar xlg = vertexData.fluidState().moleFrac(lPhaseIdx, gCompIdx);

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
            Scalar xgl = vertexData.fluidState().moleFrac(gPhaseIdx, lCompIdx);
            Scalar xgg = vertexData.fluidState().moleFrac(gPhaseIdx, gCompIdx);

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

            if (vertexData.saturation(gPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // gas phase disappears
                std::cout << "Gas phase disappears at vertex " << globalIdx
                        << ", coordinates: " << globalPos << ", Sg: "
                        << vertexData.saturation(gPhaseIdx) << std::endl;
                newPhasePresence = lPhaseOnly;

                globalSol[globalIdx][switchIdx]
                        = vertexData.fluidState().massFrac(lPhaseIdx, gCompIdx);
            }
            else if (vertexData.saturation(lPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // liquid phase disappears
                std::cout << "Liquid phase disappears at vertex " << globalIdx
                        << ", coordinates: " << globalPos << ", Sl: "
                        << vertexData.saturation(lPhaseIdx) << std::endl;
                newPhasePresence = gPhaseOnly;

                globalSol[globalIdx][switchIdx]
                        = vertexData.fluidState().massFrac(gPhaseIdx, lCompIdx);
            }
        }

        staticVertexDat_[globalIdx].phasePresence = newPhasePresence;
        staticVertexDat_[globalIdx].wasSwitched = wouldSwitch;
        return phasePresence != newPhasePresence;
    }

    // parameters given in constructor
    std::vector<StaticVertexData> staticVertexDat_;
    bool switchFlag_;
    int formulation_;
};

} // end namepace

#endif
