/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch     *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
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

#include <dumux/new_models/boxscheme/boxscheme.hh>
#include <dumux/new_models/boxscheme/p1boxtraits.hh>
#include <dumux/new_models/2p2c/2p2ctraits.hh>
#include <dumux/auxiliary/math.hh>

#include <dumux/new_models/2p2c/2p2celementdata.hh>
#include <dumux/new_models/2p2c/2p2cvertexdata.hh>
#include <dumux/new_models/2p2c/2p2cfluxdata.hh>

#include <dumux/auxiliary/apis.hh>
#include <dune/common/collectivecommunication.hh>
#include <vector>
#include <iostream>

namespace Dune
{
///////////////////////////////////////////////////////////////////////////
// TwoPTwoCBoxJacobian (evaluate the local jacobian for the newton method.)
///////////////////////////////////////////////////////////////////////////
/*!
 * \brief 2P-2C specific details needed to approximately calculate
 *        the local jacobian in the BOX scheme.
 *
 * This class is used to fill the gaps in BoxJacobian for the 2P-2C twophase flow.
 */
template<class ProblemT,
         class BoxTraitsT,
         class TwoPTwoCTraitsT,
         class ElementDataT,
         class VertexDataT,
         class FluxDataT,
         class Implementation>
class TwoPTwoCBoxJacobianBase : public BoxJacobian<ProblemT,
                                                   BoxTraitsT,
                                                   Implementation,
                                                   ElementDataT,
                                                   VertexDataT>

{
protected:
    typedef TwoPTwoCBoxJacobianBase<ProblemT,
                                    BoxTraitsT,
                                    TwoPTwoCTraitsT,
                                    ElementDataT,
                                    VertexDataT,
                                    FluxDataT,
                                    Implementation>   ThisType;
    typedef BoxJacobian<ProblemT, 
                        BoxTraitsT,
                        Implementation,
                        ElementDataT,
                        VertexDataT>    ParentType;
    
    typedef ProblemT                                Problem;
    typedef typename Problem::DomainTraits          DomTraits;
    typedef BoxTraitsT                              BoxTraits;
    typedef TwoPTwoCTraitsT                         TwoPTwoCTraits;

    enum {
        dim              = DomTraits::dim,
        dimWorld         = DomTraits::dimWorld,

        numEq            = BoxTraits::numEq,
        numPhases        = TwoPTwoCTraits::numPhases,
        numComponents    = TwoPTwoCTraits::numComponents,

        pressureIdx      = TwoPTwoCTraits::pressureIdx,
        switchIdx        = TwoPTwoCTraits::switchIdx,

        wPhase           = TwoPTwoCTraits::wPhase,
        nPhase           = TwoPTwoCTraits::nPhase,

        wComp            = TwoPTwoCTraits::wComp,
        nComp            = TwoPTwoCTraits::nComp,

        wPhaseOnly       = TwoPTwoCTraits::wPhaseOnly,
        nPhaseOnly       = TwoPTwoCTraits::nPhaseOnly,
        bothPhases       = TwoPTwoCTraits::bothPhases
    };
    static const int formulation  = TwoPTwoCTraits::formulation;
    enum {
        pWsN             = TwoPTwoCTraits::pWsN,
        pNsW             = TwoPTwoCTraits::pNsW,
    };


    typedef typename DomTraits::Scalar                Scalar;
    typedef typename DomTraits::CoordScalar           CoordScalar;
    typedef typename DomTraits::Grid                  Grid;
    typedef typename DomTraits::Vertex                Vertex;
    typedef typename DomTraits::Element               Element;
    typedef typename DomTraits::ElementIterator       ElementIterator;
    typedef typename Element::EntityPointer           ElementPointer;
    typedef typename DomTraits::LocalPosition         LocalPosition;
    typedef typename DomTraits::GlobalPosition        GlobalPosition;
    typedef typename DomTraits::VertexIterator        VertexIterator;

    typedef typename BoxTraits::SolutionVector      SolutionVector;
    typedef typename BoxTraits::FVElementGeometry   FVElementGeometry;
    typedef typename BoxTraits::SpatialFunction     SpatialFunction;
    typedef typename BoxTraits::LocalFunction       LocalFunction;
    typedef typename Grid::CollectiveCommunication  CollectiveCommunication;

    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;
    typedef ElementDataT         ElementData;
    typedef VertexDataT          VertexData;
    typedef FluxDataT            FluxData;

    typedef std::vector<VertexData> VertexDataArray;

    typedef FieldMatrix<Scalar, dim, dim>  Tensor;

    static const Scalar upwindAlpha = TwoPTwoCTraits::upwindAlpha;

    /*!
     * \brief Data which is attached to each vertex and is not only
     *        stored locally.
     */
    struct StaticVertexData {
        int phaseState;
        int oldPhaseState;
    };

public:
    TwoPTwoCBoxJacobianBase(ProblemT &problem)
        : ParentType(problem),
          staticVertexDat_(problem.numVertices())
    {
        switchFlag_ = false;
    };


    /*!
     * \brief Evaluate the amount all conservation quantites
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     */
    void computeStorage(SolutionVector &result, int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const VertexDataArray &elemDat = usePrevSol ? this->prevElemDat_  : this->curElemDat_;
        const VertexData  &vertDat = elemDat[scvIdx];
        
        // compute storage term of all components within all phases
        result = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
                result[compIdx] +=
                    vertDat.density[phaseIdx]*
                    vertDat.saturation[phaseIdx]*
                    vertDat.massfrac[compIdx][phaseIdx];
        result *= vertDat.porosity;
    }

    /*!
     * \brief Evaluates the total flux of all conservation quantities
     *        over a face of a subcontrol volume.
     */
    void computeFlux(SolutionVector &flux, int faceIdx) const
    {
        FluxData vars(this->problem_,
                      this->curElement_(),
                      this->curElementGeom_,
                      faceIdx,
                      this->curElemDat_);
        
        flux = 0;
        asImp_()->computeAdvectiveFlux(flux, vars);
        asImp_()->computeDiffusiveFlux(flux, vars);
    }
    
    /*!
     * \brief Evaluates the advective mass flux of all components over
     *        a face of a subcontrol volume.
     */
    void computeAdvectiveFlux(SolutionVector &flux, 
                              const FluxData &vars) const
    {
        ////////
        // advective fluxes of all components in all phases
        ////////
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // data attached to upstream and the downstream vertices
            // of the current phase
            const VertexData &up = this->curElemDat_[vars.upstreamIdx[phaseIdx]];
            const VertexData &dn = this->curElemDat_[vars.downstreamIdx[phaseIdx]];

            for (int  compIdx = 0; compIdx < numComponents; ++compIdx) {
                // add advective flux of current component in current
                // phase
                flux[compIdx] +=  
                    vars.vDarcyNormal[phaseIdx] * (
                        upwindAlpha* // upstream vertex
                        (  up.density[phaseIdx] *
                           up.mobility[phaseIdx] *
                           up.massfrac[compIdx][phaseIdx])
                        +
                        (1 - upwindAlpha)* // downstream vertex
                        (  dn.density[phaseIdx] *
                           dn.mobility[phaseIdx] *
                           dn.massfrac[compIdx][phaseIdx]));
            }
        }
    }

    /*!
     * \brief Adds the diffusive mass flux of all components over
     *        a face of a subcontrol volume.
     */
    void computeDiffusiveFlux(SolutionVector &flux, const FluxData &vars) const
    {        
        // add diffusive flux of non-wetting component in wetting phase
        Scalar tmp = 
            vars.diffCoeffPM[wPhase] * vars.densityAtIP[wPhase] * 
            (vars.concentrationGrad[wPhase]*vars.face->normal);
        flux[nComp] += tmp; 
        flux[wComp] -= tmp;
        
        // add diffusive flux of wetting component in non-wetting phase
        tmp = vars.diffCoeffPM[nPhase] * vars.densityAtIP[nPhase] * 
            (vars.concentrationGrad[nPhase]*vars.face->normal);;
        flux[wComp] += tmp;
        flux[nComp] -= tmp;

        // TODO: the diffusive flux of the wetting component in the
        // wetting phase does rarly exhibit the same mass as the flux
        // of the non-wetting component, which means that it is not
        // -tmp
    }

    /*!
     * \brief Calculate the source term of the equation
     */
    void computeSource(SolutionVector &q, int localVertexIdx)
    {
        this->problem_.source(q,
                              this->curElement_(),
                              this->curElementGeom_,
                              localVertexIdx);
    }

    /*!
     * \brief Return the temperature given the solution vector of a
     *        finite volume.
     */
    template <class SolutionVector>
    Scalar temperature(const SolutionVector &sol)
    { return this->problem_.temperature(); /* constant temperature */ }

    /*!
     * \brief Initialize the static data with the initial solution.
     *
     * Called by TwoPTwoCBoxModel::initial()
     */
    void initStaticData()
    {
        setSwitched(false);

        VertexIterator it = this->problem_.vertexBegin();
        VertexIterator endit = this->problem_.vertexEnd();
        for (; it != endit; ++it)
        {
            int globalIdx = this->problem_.vertexIdx(*it);
            const GlobalPosition &globalPos = it->geometry().corner(0);
            
            // initialize phase state
            staticVertexDat_[globalIdx].phaseState =
                this->problem_.initialPhaseState(*it, globalIdx, globalPos);
            staticVertexDat_[globalIdx].oldPhaseState =
                staticVertexDat_[globalIdx].phaseState;
        }
    }

    /*!
     * \brief Update the static data of all vertices in the grid.
     */
    void updateStaticData(SpatialFunction &curGlobalSol, SpatialFunction &oldGlobalSol)
    {
        bool wasSwitched = false;

        VertexIterator it = this->problem_.vertexBegin();
        for (; it != this->problem_.vertexEnd(); ++it)
        {
            int globalIdx = this->problem_.vertexIdx(*it);
            const GlobalPosition &global = it->geometry().corner(0);
            
            wasSwitched = primaryVarSwitch_(curGlobalSol,
                                            globalIdx,
                                            global)
                || wasSwitched;
        }

        // make sure that if there was a variable switch in an
        // other partition we will also set the switch flag
        // for our partition.
        wasSwitched = this->problem_.grid().comm().max(wasSwitched);

        setSwitched(wasSwitched);
    }

    /*!
     * \brief Set the old phase of all verts state to the current one.
     */
    void updateOldPhaseState()
    {
        int numVertices = this->problem_.numVertices();
        for (int i = 0; i < numVertices; ++i)
            staticVertexDat_[i].oldPhaseState = staticVertexDat_[i].phaseState;
    }

    /*!
     * \brief Returns the phase state of the current or the old solution of a vertex.
     */
    int phaseState(int globalVertexIdx, bool oldSol) const
    { 
        return
            oldSol?
            staticVertexDat_[globalVertexIdx].oldPhaseState :
            staticVertexDat_[globalVertexIdx].phaseState;
    }

    /*!
     * \brief Reset the current phase state of all vertices to the old one.
     *
     * This is done after an update failed.
     */
    void resetPhaseState()
    {
        int numVertices = this->problem_.numVertices();
        for (int i = 0; i < numVertices; ++i)
            staticVertexDat_[i].phaseState = staticVertexDat_[i].oldPhaseState;
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
    void calculateMass(const SpatialFunction &globalSol, Dune::FieldVector<Scalar, 4> &mass)
    {
        ElementIterator elementIt = this->problem_.elementBegin();
        ElementIterator endit = this->problem_.elementEnd();
        unsigned numVertices = this->problem_.numVertices();
        LocalFunction curSol(numVertices);
        ElementData   elemDat;
        VertexData tmp;
        int state;
        Scalar vol, poro, rhoN, rhoW, satN, satW, xAW, xWW, xWN, xAN, pW, Te;
        Scalar massNComp(0.), massNCompNPhase(0.), massWComp(0.), massWCompWPhase(0.);

        mass = 0;
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
            setCurrentElement(*elementIt);
            this->restrictToElement(curSol, globalSol);
            this->updateElementData_(elemDat, curSol, false);
            // get geometry type

            int numLocalVerts = elementIt->template count<dim>();

            // Loop over element vertices
            for (int i = 0; i < numLocalVerts; ++i)
            {
                int globalIdx = this->problem_.vertexIdx(*elementIt, i);

                vol = this->curElementGeom_.subContVol[i].volume;

                state =  staticVertexDat_[globalIdx].phaseState;
                poro = this->problem_.porosity(this->curElement_(), i);
                rhoN = elemDat[i].density[nPhase];
                rhoW = elemDat[i].density[wPhase];
                satN = elemDat[i].saturation[nPhase];
                satW = elemDat[i].saturation[wPhase];
                xAW = elemDat[i].massfrac[nComp][wPhase];
                xWW = elemDat[i].massfrac[wComp][wPhase];
                xWN = elemDat[i].massfrac[wComp][nPhase];
                xAN = elemDat[i].massfrac[nComp][nPhase];
                pW = elemDat[i].pressure[wPhase];
                Te = Implementation::temperature_((*globalSol)[globalIdx]);
                massNComp = vol * poro * (satN * rhoN * xAN + satW * rhoW * xAW);
                massNCompNPhase = vol * poro * satN * rhoN * xAN;
                massWComp = vol * poro * (satW * rhoW * xWW + satN * rhoN * xWN);
                massWCompWPhase = vol * poro * satW * rhoW * xWW;

                // get minimum and maximum values of primary variables
                minSat = std::min(minSat, satN);
                maxSat = std::max(maxSat, satN);
                minP = std::min(minP, pW);
                maxP = std::max(maxP, pW);
                minX = std::min(minX, xAW);
                maxX = std::max(maxX, xAW);
                minTe = std::min(minTe, Te);
                maxTe = std::max(maxTe, Te);

                // calculate total mass
                mass[0] += massNComp;       // total mass of nonwetting component
                mass[1] += massNCompNPhase; // mass of nonwetting component in nonwetting phase
                mass[2] += massWComp;       // total mass of wetting component
                mass[3] += massWCompWPhase; // mass of wetting component in wetting phase
            }
        }
        
        // IF PARALLEL: calculate total mass including all processors
        // also works for sequential calculation
        mass = this->problem_.grid().comm().sum(mass);
        
        if(this->problem_.grid().comm() == 0) // IF PARALLEL: only print by processor with rank() == 0
        {
            // print minimum and maximum values
            std::cout << "nonwetting phase saturation: min = "<< minSat
                      << ", max = "<< maxSat << std::endl;
            std::cout << "wetting phase pressure: min = "<< minP
                      << ", max = "<< maxP << std::endl;
            std::cout << "mass fraction nComp: min = "<< minX
                      << ", max = "<< maxX << std::endl;
            std::cout << "temperature: min = "<< minTe
                      << ", max = "<< maxTe << std::endl;
        }
    }
    
    /*!
     * \brief Add the mass fraction of air in water to VTK output of
     *        the current timestep.
     */
    template <class MultiWriter>
    void addVtkFields(MultiWriter &writer, const SpatialFunction &globalSol)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->problem_.numVertices();
        unsigned numElements = this->problem_.numElements();
        ScalarField *pW =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *pN =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *pC =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *Sw =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *Sn =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *rhoW =         writer.template createField<Scalar, 1>(numVertices);
        ScalarField *rhoN =         writer.template createField<Scalar, 1>(numVertices);
        ScalarField *mobW =         writer.template createField<Scalar, 1>(numVertices);
        ScalarField *mobN =         writer.template createField<Scalar, 1>(numVertices);
        ScalarField *massfracAinW = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *massfracAinN = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *massfracWinW = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *massfracWinN = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *temperature  = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *phaseState   = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *velocityX    = writer.template createField<Scalar, 1>(numElements);
        ScalarField *velocityY    = writer.template createField<Scalar, 1>(numElements);
        ScalarField *velocityZ    = writer.template createField<Scalar, 1>(numElements);

        LocalFunction   tmpSol;
        VertexDataArray elemDat(BoxTraits::ShapeFunctionSetContainer::maxsize);
        VertexData      tmp;

        ElementIterator elementIt = this->problem_.elementBegin();
        ElementIterator endit = this->problem_.elementEnd();

        for (; elementIt != endit; ++elementIt)
        {
            int numLocalVerts = elementIt->template count<dim>();
            tmpSol.resize(numLocalVerts);

            setCurrentElement(*elementIt);
            this->restrictToElement(tmpSol, globalSol);
            updateElementData_(elemDat, tmpSol, false);

            for (int i = 0; i < numLocalVerts; ++i)
            {
                int globalIdx = this->problem_.vertexIdx(*elementIt, i);

                (*pW)[globalIdx] = elemDat[i].pressure[wPhase];
                (*pN)[globalIdx] = elemDat[i].pressure[nPhase];
                (*pC)[globalIdx] = elemDat[i].pC;
                (*Sw)[globalIdx] = elemDat[i].saturation[wPhase];
                (*Sn)[globalIdx] = elemDat[i].saturation[nPhase];
                (*rhoW)[globalIdx] = elemDat[i].density[wPhase];
                (*rhoN)[globalIdx] = elemDat[i].density[nPhase];
                (*mobW)[globalIdx] = elemDat[i].mobility[wPhase];
                (*mobN)[globalIdx] = elemDat[i].mobility[nPhase];
                (*massfracAinW)[globalIdx] = elemDat[i].massfrac[nComp][wPhase];
                (*massfracAinN)[globalIdx] = elemDat[i].massfrac[nComp][nPhase];
                (*massfracWinW)[globalIdx] = elemDat[i].massfrac[wComp][wPhase];
                (*massfracWinN)[globalIdx] = elemDat[i].massfrac[wComp][nPhase];
                (*temperature)[globalIdx] = asImp_()->temperature((*globalSol)[globalIdx]);
                (*phaseState)[globalIdx] = staticVertexDat_[globalIdx].phaseState;
            };

            // Vector containing the velocity at the element
            GlobalPosition velocity[numPhases];
            GlobalPosition elementVelocity[numPhases];

            // loop over the phases
            for (int phase=0; phase < numPhases; phase++)
            {
                elementVelocity[phase] = 0;

                int elementIdx = this->problem_.elementIdx(*elementIt);
                for (int faceIdx = 0; faceIdx< this->curElementGeom_.numEdges; faceIdx++)
                {
                    velocity[phase] = 0;
                    /* asImp_().calculateDarcyVelocity(velocity[phase],
                                                     faceIdx);
                    */
                    elementVelocity[phase] += velocity[phase];
                }
                elementVelocity[phase] *= 1.0/this->curElementGeom_.numEdges;
                (*velocityX)[elementIdx] = elementVelocity[0][phase];
                if (dim >= 2)
                    (*velocityY)[elementIdx] = elementVelocity[1][phase];
                if (dim == 3)
                    (*velocityZ)[elementIdx] = elementVelocity[2][phase];
            }
        }

        writer.addVertexData(pW, "pW");
        writer.addVertexData(pN, "pN");
        writer.addVertexData(pC, "pC");
        writer.addVertexData(Sw, "SW");
        writer.addVertexData(Sn, "SN");
        writer.addVertexData(rhoW, "rhoW");
        writer.addVertexData(rhoN, "rhoN");
        writer.addVertexData(mobW, "mobW");
        writer.addVertexData(mobN, "mobN");
        writer.addVertexData(massfracAinW, "XaW");
        writer.addVertexData(massfracAinN, "XaN");
        writer.addVertexData(massfracWinW, "XwW");
        writer.addVertexData(massfracWinN, "XwN");
        writer.addVertexData(temperature, "T");
        writer.addVertexData(phaseState, "phase state");
        writer.addCellData(velocityX, "Vx");
        if (dim >= 2)
            writer.addCellData(velocityY, "Vy");
        if (dim == 3)
            writer.addCellData(velocityZ, "Vz");
    }

    /*!
     * \brief Reads the current solution for a vertex from a restart
     *        file.
     */
    void deserializeEntity(std::istream &inStream,
                           const Vertex &vert)
    {
        int vertIdx = this->problem_.vertexIdx(vert);

        // read phase state
        if (!inStream.good()) {
            DUNE_THROW(IOError,
                       "Could not deserialize vertex "
                       << vertIdx);
        }

        inStream >> staticVertexDat_[vertIdx].phaseState;
        staticVertexDat_[vertIdx].oldPhaseState
            = staticVertexDat_[vertIdx].phaseState;
    };

    /*!
     * \brief Write the current phase state of an vertex to a restart
     *        file.
     */
    void serializeEntity(std::ostream &outStream,
                         const Vertex &vert)
    {
        int vertIdx = this->problem_.vertexIdx(vert);

        if (!outStream.good()) {
            DUNE_THROW(IOError,
                       "Could not serialize vertex "
                       << vertIdx);
        }
        
        outStream << staticVertexDat_[vertIdx].phaseState
                  << " ";
    };


protected:
    Implementation *asImp_()
    { return static_cast<Implementation *>(this); }
    const Implementation *asImp_() const
    { return static_cast<const Implementation *>(this); }


    //  perform variable switch at a vertex; Returns true if a
    //  variable switch was performed.
    bool primaryVarSwitch_(SpatialFunction &globalSol,
                           int globalIdx,
                           const GlobalPosition &globalPos)
    {
        // evaluate primary variable switch
        int phaseState    = staticVertexDat_[globalIdx].phaseState;
        int newPhaseState = phaseState;
        Scalar temperature = asImp_()->temperature((*globalSol)[globalIdx]);
        // calculate saturations, phase pressures and mass fractions
        static VertexData vertexData;
        static LocalPosition localPos(0.0);
        vertexData.updateSaturations((*globalSol)[globalIdx], phaseState);
        Scalar pC = this->problem_.materialLaw().pC(vertexData.saturation[wPhase],
                                                    globalPos,
                                                    * (this->problem_.elementBegin()), // HACK
                                                    localPos, 
                                                    temperature);
        vertexData.updatePressures((*globalSol)[globalIdx], pC);
        vertexData.updateMassFracs((*globalSol)[globalIdx],
                                   this->problem_.multicomp(),
                                   phaseState,
                                   temperature);

        // check if a primary var switch is necessary
        if (phaseState == nPhaseOnly)
        {
            Scalar xWNmax = this->problem_.multicomp().xWN(vertexData.pressure[nPhase], temperature);
            if (vertexData.massfrac[wComp][nPhase] > xWNmax)
            {
                // wetting phase appears
                std::cout << "wetting phase appears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << std::endl;
                newPhaseState = bothPhases;
                if (formulation == pNsW) 
                    (*globalSol)[globalIdx][switchIdx] = 0.0;
                else if (formulation == pWsN)
                    (*globalSol)[globalIdx][switchIdx] = 1.0;
            };
        }
        else if (phaseState == wPhaseOnly)
        {
            Scalar xAWmax = this->problem_.multicomp().xAW(vertexData.pressure[wPhase], temperature);
            if (vertexData.massfrac[nComp][wPhase] > xAWmax)
            {
                // non-wetting phase appears
                std::cout << "Non-wetting phase appears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << std::endl;
                newPhaseState = bothPhases;
                if (formulation == pNsW)
                    (*globalSol)[globalIdx][switchIdx] = 1.0;
                else if (formulation == pWsN)
                    (*globalSol)[globalIdx][switchIdx] = 0.0;
            }
        }
        else if (phaseState == bothPhases) {
            if (vertexData.saturation[nPhase] <= 0) {
                // non-wetting phase disappears
                std::cout << "Non-wetting phase disappears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << std::endl;
                newPhaseState = wPhaseOnly;
                (*globalSol)[globalIdx][switchIdx]
                    = this->problem_.multicomp().xAW(vertexData.pressure[nPhase], temperature);
            }
            else if (vertexData.saturation[wPhase] <= 0) {
                // wetting phase disappears
                std::cout << "Wetting phase disappears at vertex " << globalIdx
                          << ", coordinates: " << globalPos << std::endl;
                newPhaseState = nPhaseOnly;
                (*globalSol)[globalIdx][switchIdx]
                    = this->problem_.multicomp().xWN(vertexData.pressure[nPhase], temperature);
            }
        }

        staticVertexDat_[globalIdx].phaseState = newPhaseState;

        return phaseState != newPhaseState;
    }

    // parameters given in constructor
    std::vector<StaticVertexData> staticVertexDat_;
    bool                          switchFlag_;
    int                           formulation_;
};


} // end namepace

#endif
