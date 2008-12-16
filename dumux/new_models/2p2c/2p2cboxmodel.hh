/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
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
#ifndef DUMUX_NEW_2P2C_BOX_MODEL_HH
#define DUMUX_NEW_2P2C_BOX_MODEL_HH

#include <dumux/new_models/boxscheme/boxscheme.hh>
#include <dumux/new_models/boxscheme/p1boxtraits.hh>

#include <dumux/auxiliary/apis.hh>

#include <vector>

namespace Dune
{
    ///////////////////////////////////////////////////////////////////////////
    // two-phase two-component traits (central place for names and
    // indices required by the TwoPTwoCBoxJacobian and TwoPTwoCBoxModel)
    ///////////////////////////////////////////////////////////////////////////
    /*!
     * \brief The 2P-2C specific traits.
     */
    template <class Scalar>
    class TwoPTwoCTraits
    {
    public:
        enum {
            numEq		  = 2,  //!< Number of primary variables / equations
            numPhases     = 2,  //!< Number of fluid phases
            numComponents = 2   //!< Number of fluid components within a phase
        };
        enum { // Primary variable indices
            pWIdx     = 0,       //!< Idx for the wetting phase pressure in a field vector
            switchIdx = 1,       //!< Idx for the non-wetting phase quantity
        };
        enum { // Phase Indices
            wPhase = 0,  //!< Idx of the wetting phase
            nPhase = 1   //!< Idx of the non-wetting phase
        };
        enum { // Component indices
            wComp = 0,  //!< Idx of the major component of the wetting phase
            nComp = 1   //!< Idx of the major component of the non-wetting phase
        };
        enum { // present phases
            nPhaseOnly = 0, //!< Only the non-wetting phase is present
            wPhaseOnly = 1, //!< Only the wetting phase is present
            bothPhases = 2  //!< Both phases are present
        };

        typedef FieldVector<Scalar, numPhases>         PhasesVector;

        /*!
         * \brief Data which is attached to each vert of the and can
         *        be shared between multiple calculations and should
         *        thus be cached in order to increase efficency.
         */
        struct VariableVertexData
        {
            Scalar satN;
            Scalar satW;
            Scalar pW;
            Scalar pC;
            Scalar pN;
            PhasesVector mobility;  //FieldVector with the number of phases
            PhasesVector density;
            Dune::FieldMatrix<Scalar, numComponents, numPhases> massfrac;
            int phaseState;
        };

    };


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
             class Implementation>
    class TwoPTwoCBoxJacobianBase : public BoxJacobian<ProblemT,
                                                       BoxTraitsT,
                                                       Implementation >
    {
    protected:
        typedef TwoPTwoCBoxJacobianBase<ProblemT,
                                        BoxTraitsT,
                                        TwoPTwoCTraitsT,
                                        Implementation>              ThisType;
        typedef BoxJacobian<ProblemT, BoxTraitsT, Implementation>    ParentType;

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

            pWIdx            = TwoPTwoCTraits::pWIdx,
            switchIdx        = TwoPTwoCTraits::switchIdx,

            wPhase           = TwoPTwoCTraits::wPhase,
            nPhase           = TwoPTwoCTraits::nPhase,

            wComp            = TwoPTwoCTraits::wComp,
            nComp            = TwoPTwoCTraits::nComp,

            wPhaseOnly       = TwoPTwoCTraits::wPhaseOnly,
            nPhaseOnly       = TwoPTwoCTraits::nPhaseOnly,
            bothPhases       = TwoPTwoCTraits::bothPhases
        };


        typedef typename DomTraits::Scalar                Scalar;
        typedef typename DomTraits::CoordScalar           CoordScalar;
        typedef typename DomTraits::Grid                  Grid;
        typedef typename DomTraits::Element               Element;
        typedef typename DomTraits::ElementIterator       ElementIterator;
        typedef typename Element::EntityPointer           ElementPointer;
        typedef typename DomTraits::LocalPosition         LocalPosition;
        typedef typename DomTraits::GlobalPosition            GlobalPosition;
        typedef typename DomTraits::VertexIterator        VertexIterator;

        typedef typename BoxTraits::SolutionVector      SolutionVector;
        typedef typename BoxTraits::FVElementGeometry   FVElementGeometry;
        typedef typename BoxTraits::SpatialFunction     SpatialFunction;
        typedef typename BoxTraits::LocalFunction       LocalFunction;

        typedef typename TwoPTwoCTraits::PhasesVector        PhasesVector;
        typedef typename TwoPTwoCTraits::VariableVertexData  VariableVertexData;
        typedef FieldMatrix<Scalar, dim, dim>  Tensor;

        /*!
         * \brief Cached data for the each vert of the element.
         */
        struct ElementData
        {
            VariableVertexData vertex[BoxTraits::ShapeFunctionSetContainer::maxsize];
        };

        /*!
         * \brief Data which is attached to each vert and is not only
         *        locally.
         */
        struct StaticVertexData {
            int phaseState;
            int oldPhaseState;
        };

        void updateVarVertexData_(VariableVertexData &d,
                                const SolutionVector &vertSol,
                                int phaseState,
                                const Element &element,
                                int localIdx,
                                Problem &problem,
                                Scalar temperature) const
                {
                    const GlobalPosition &global = element.geometry()[localIdx];
                    const LocalPosition &local =
                        DomTraits::referenceElement(element.type()).position(localIdx,
                                                                          dim);
                    d.pW = vertSol[pWIdx];
                    if (phaseState == bothPhases) d.satN = vertSol[switchIdx];
                    else if (phaseState == wPhaseOnly) d.satN = 0.0;
                    else if (phaseState == nPhaseOnly) d.satN = 1.0;
                    else DUNE_THROW(Dune::InvalidStateException, "Phase state " << phaseState << " is invalid.");

                    d.satW = 1.0 - d.satN;
                    d.pC = problem.materialLaw().pC(d.satW,
                                                    global,
                                                    element,
                                                    local);
                    d.pN = d.pW + d.pC;

                    // Solubilities of components in phases
                    if (phaseState == bothPhases) {
                        d.massfrac[nComp][wPhase] = problem.multicomp().xAW(d.pN, temperature);
                        d.massfrac[wComp][nPhase] = problem.multicomp().xWN(d.pN, temperature);
                    }
                    else if (phaseState == wPhaseOnly) {
                        d.massfrac[wComp][nPhase] = 0.0;
                        d.massfrac[nComp][wPhase] =  vertSol[switchIdx];
                    }
                    else if (phaseState == nPhaseOnly){
                        d.massfrac[wComp][nPhase] = vertSol[switchIdx];
                        d.massfrac[nComp][wPhase] = 0.0;
                    }
                    else DUNE_THROW(Dune::InvalidStateException, "Phase state " << phaseState << " is invalid.");

                    d.massfrac[wComp][wPhase] = 1.0 - d.massfrac[nComp][wPhase];
                    d.massfrac[nComp][nPhase] = 1.0 - d.massfrac[wComp][nPhase];
                    d.phaseState = phaseState;

                    // Density of Water is set constant here!
                    d.density[wPhase] = problem.wettingPhase().density(temperature,
                                                                       d.pW,
                                                                       d.massfrac[nComp][wPhase]);
                    d.density[nPhase] = problem.nonwettingPhase().density(temperature,
                                                                          d.pN,
                                                                          d.massfrac[wComp][nPhase]);

                    // Mobilities
                    d.mobility[wPhase] = problem.materialLaw().mobW(d.satW, 
                                                                    global,
                                                                    element, 
                                                                    local, 
                                                                    temperature, 
                                                                    d.pW);
                    d.mobility[nPhase] = problem.materialLaw().mobN(d.satN,
                                                                    global, 
                                                                    element,
                                                                    local,
                                                                    temperature,
                                                                    d.pN);
                }

    public:
        TwoPTwoCBoxJacobianBase(ProblemT &problem)
            : ParentType(problem),
              staticVertexDat_(problem.numVertices())
            {
                switchFlag_ = false;
            };

        /*!
         * \brief Set the current grid element.
         */
        void setCurrentElement(const Element &element)
            {
                ParentType::setCurrentElement_(element);
            };

        /*!
         * \brief Set the parameters for the calls to the remaining
         *        members.
         */
        void setParams(const Element &element, LocalFunction &curSol, LocalFunction &prevSol)
            {
                setCurrentElement(element);

                // TODO: scheme which allows not to copy curSol and
                // prevSol all the time
                curSol_ = curSol;
                updateElementData_(curElemDat_, curSol_, false);
                curSolDeflected_ = false;

                prevSol_ = prevSol;
                updateElementData_(prevElemDat_, prevSol_, true);
            };

        /*!
         * \brief Vary a single component of a single vert of the
         *        local solution for the current element.
         *
         * This method is a optimization, since if varying a single
         * component at a degree of freedom not the whole element cache
         * needs to be recalculated. (Updating the element cache is very
         * expensive since material laws need to be evaluated.)
         */
        void deflectCurSolution(int vert, int component, Scalar value)
            {
                // make sure that the orignal state can be restored
                if (!curSolDeflected_) {
                    curSolDeflected_ = true;

                    curSolOrigValue_ = curSol_[vert][component];
                    curSolOrigVarData_ = curElemDat_.vertex[vert];
                }

                int globalIdx = ParentType::problem_.vertIdx(ParentType::curElement_(),
                                                               vert);

                curSol_[vert][component] = value;
                asImp_()->updateVarVertexData_(curElemDat_.vertex[vert],
                                             curSol_[vert],
                                             staticVertexDat_[globalIdx].phaseState,
                                             this->curElement_(),
                                             vert,
                                             this->problem_,
                                             Implementation::temperature_(curSol_[vert]));
            }

        /*!
         * \brief Restore the local jacobian to the state before
         *        deflectCurSolution() was called.
         *
         * This only works if deflectSolution was only called with
         * (vert, component) as arguments.
         */
        void restoreCurSolution(int vert, int component)
            {
                curSolDeflected_ = false;
                curSol_[vert][component] = curSolOrigValue_;
                curElemDat_.vertex[vert] = curSolOrigVarData_;
            };

        /*!
         * \brief Evaluate the rate of change of all conservation
         *        quantites (e.g. phase mass) within a sub control
         *        volume of a finite volume element in the 2P-2C
         *        model.
         *
         * This function should not include the source and sink terms.
         */
        void computeStorage(SolutionVector &result, int scvId, bool usePrevSol) const
            {
                result = Scalar(0);

                const LocalFunction &sol   = usePrevSol ? this->prevSol_ : this->curSol_;
                const ElementData &elementCache = usePrevSol ? prevElemDat_  : curElemDat_;

                Scalar satN = elementCache.vertex[scvId].satN;
                Scalar satW = elementCache.vertex[scvId].satW;

                // assume porosity defined at verts
                Scalar porosity =
                    this->problem_.porosity(this->curElement_(), scvId);

                // storage of component water
                result[pWIdx] =
                    porosity*(elementCache.vertex[scvId].density[wPhase]*
                              satW*
                              elementCache.vertex[scvId].massfrac[wComp][wPhase]
                              + elementCache.vertex[scvId].density[nPhase]*
                              satN*
                              elementCache.vertex[scvId].massfrac[wComp][nPhase]);

                // storage of component air
                result[switchIdx] =
                    porosity*(elementCache.vertex[scvId].density[nPhase]*
                              satN*
                              elementCache.vertex[scvId].massfrac[nComp][nPhase]
                              + elementCache.vertex[scvId].density[wPhase]*
                              satW*
                              elementCache.vertex[scvId].massfrac[nComp][wPhase]);

                // storage of energy
                asImp_()->heatStorage(result, scvId, sol, elementCache);
            }

        /*!
         * \brief Evaluates the mass flux over a face of a subcontrol
         *        volume.
         */
        void computeFlux(SolutionVector &flux, int faceId) const
            {
                // set flux vector to zero
                int i = this->curElementGeom_.subContVolFace[faceId].i;
                int j = this->curElementGeom_.subContVolFace[faceId].j;

                // normal vector, value of the area of the scvf
                const GlobalPosition &normal(this->curElementGeom_.subContVolFace[faceId].normal);

                // get global coordinates of verts i,j
                const GlobalPosition &global_i = this->curElementGeom_.subContVol[i].global;
                const GlobalPosition &global_j = this->curElementGeom_.subContVol[j].global;

                // get local coordinates of verts i,j
                const LocalPosition &local_i = this->curElementGeom_.subContVol[i].local;
                const LocalPosition &local_j = this->curElementGeom_.subContVol[j].local;

                GlobalPosition pGrad[numPhases];
                GlobalPosition xGrad[numPhases];
                GlobalPosition tempGrad(0.0);
                for (int k = 0; k < numPhases; ++k) {
                    pGrad[k] = Scalar(0);
                    xGrad[k] = Scalar(0);
                }

                GlobalPosition tmp(0.0);
                PhasesVector pressure(0.0), massfrac(0.0);

                // calculate FE gradient (grad p for each phase)
                for (int k = 0; k < this->curElementGeom_.numVertices; k++) // loop over adjacent verts
                {
                    // FEGradient at vert k
                    const LocalPosition &feGrad = this->curElementGeom_.subContVolFace[faceId].grad[k];

                    pressure[wPhase] = this->curElemDat_.vertex[k].pW;
                    pressure[nPhase] = this->curElemDat_.vertex[k].pN;

                    // compute sum of pressure gradients for each phase
                    for (int phase = 0; phase < numPhases; phase++)
                    {
                        tmp = feGrad;
                        tmp *= pressure[phase];
                        pGrad[phase] += tmp;
                    }

                    // for diffusion of air in wetting phase
                    tmp = feGrad;
                    tmp *= this->curElemDat_.vertex[k].massfrac[nComp][wPhase];
                    xGrad[wPhase] += tmp;

                    // for diffusion of water in nonwetting phase
                    tmp = feGrad;
                    tmp *= this->curElemDat_.vertex[k].massfrac[wComp][nPhase];
                    xGrad[nPhase] += tmp;

                    // temperature gradient
                    asImp_()->updateTempGrad(tempGrad, feGrad, this->curSol_, k);
                }

                // correct the pressure gradients by the hydrostatic
                // pressure due to gravity
                for (int phase=0; phase < numPhases; phase++)
                {
                    tmp = this->problem_.gravity();
                    tmp *= this->curElemDat_.vertex[i].density[phase];
                    pGrad[phase] -= tmp;
                }


                // calculate the permeability tensor
                Tensor K         = this->problem_.soil().K(global_i, ParentType::curElement_(), local_i);
                const Tensor &Kj = this->problem_.soil().K(global_j, ParentType::curElement_(), local_j);
                harmonicMeanK_(K, Kj);

                // magnitute of darcy velocity of each phase projected
                // on the normal of the sub-control volume's face
                PhasesVector vDarcyOut;
                // temporary vector for the Darcy velocity
                GlobalPosition vDarcy;
                for (int phase=0; phase < numPhases; phase++)
                {
                    K.mv(pGrad[phase], vDarcy);  // vDarcy = K * grad p
                    vDarcyOut[phase] = vDarcy*normal;
                }

                // find upsteam and downstream verts
                const VariableVertexData *upW = &(this->curElemDat_.vertex[i]);
                const VariableVertexData *dnW = &(this->curElemDat_.vertex[j]);
                const VariableVertexData *upN = &(this->curElemDat_.vertex[i]);
                const VariableVertexData *dnN = &(this->curElemDat_.vertex[j]);
                if (vDarcyOut[wPhase] > 0) {
                    std::swap(upW, dnW);
                };
                if (vDarcyOut[nPhase] > 0)  {
                    std::swap(upN, dnN);
                };

                // Upwind parameter
                Scalar alpha = 1.0; // -> use only the upstream vert

                ////////
                // advective flux of the wetting component
                ////////

                // flux in the wetting phase
                flux[pWIdx] =  vDarcyOut[wPhase] * (
                    alpha* // upstream verts
                    (  upW->density[wPhase] *
                       upW->mobility[wPhase] *
                       upW->massfrac[wComp][wPhase])
                    +
                    (1-alpha)* // downstream verts
                    (  dnW->density[wPhase] *
                       dnW->mobility[wPhase] *
                       dnW->massfrac[wComp][wPhase]));

                // flux in the non-wetting phase
                flux[pWIdx] += vDarcyOut[nPhase] * (
                    alpha* // upstream vert
                    (  upN->density[nPhase] *
                       upN->mobility[nPhase] *
                       upN->massfrac[wComp][nPhase])
                    +
                    (1-alpha)* // downstream vert
                    (  dnN->density[nPhase] *
                       dnN->mobility[nPhase] *
                       dnN->massfrac[wComp][nPhase]) );

                ////////
                // advective flux of the non-wetting component
                ////////

                // flux in the wetting phase
                flux[switchIdx]   = vDarcyOut[nPhase] * (
                    alpha * // upstream verts
                    (  upN->density[nPhase] *
                       upN->mobility[nPhase] *
                       upN->massfrac[nComp][nPhase])
                    +
                    (1-alpha) * // downstream vert
                    (  dnN->density[nPhase] *
                       dnN->mobility[nPhase] *
                       dnN->massfrac[nComp][nPhase]) );

                // flux in the non-wetting phase
                flux[switchIdx]  += vDarcyOut[wPhase] * (
                    alpha * // upstream vert
                    (  upW->density[wPhase] *
                       upW->mobility[wPhase] *
                       upW->massfrac[nComp][wPhase])
                    +
                    (1-alpha) * // downstream vert
                    (  dnW->density[wPhase] *
                       dnW->mobility[wPhase] *
                       dnW->massfrac[nComp][wPhase]) );

                ////////
                // advective flux of energy
                ////////
                asImp_()->advectiveHeatFlux(flux, vDarcyOut, alpha, upW, dnW, upN, dnN);

                asImp_()->diffusiveHeatFlux(flux, faceId, tempGrad);

                return;

                // DIFFUSION
                SolutionVector normDiffGrad;

                // get local to global id map
                int state_i = this->curElemDat_.vertex[i].phaseState;
                int state_j = this->curElemDat_.vertex[j].phaseState;

                Scalar diffusionWW(0.0), diffusionWN(0.0); // diffusion of water
                Scalar diffusionAW(0.0), diffusionAN(0.0); // diffusion of air
                SolutionVector avgDensity, avgDpm;

                // Diffusion coefficent
                // TODO: needs to be continuously dependend on the phase saturations
                avgDpm[wPhase]=2e-9;
                avgDpm[nPhase]=2.25e-5;
                if (state_i == nPhaseOnly || state_j == nPhaseOnly)
                {
                    // only the nonwetting phase is present in at
                    // least one element -> no diffusion within the
                    // wetting phase
                    avgDpm[wPhase] = 0;
                }
                if (state_i == wPhaseOnly || state_j == wPhaseOnly)
                {
                    // only the wetting phase is present in at least
                    // one element -> no diffusion within the non wetting
                    // phase
                    avgDpm[nPhase] = 0;
                }

                // length of the diffusion gradient
                normDiffGrad[wPhase] = xGrad[wPhase]*normal;
                normDiffGrad[nPhase] = xGrad[nPhase]*normal;

                // calculate the arithmetic mean of densities
                avgDensity[wPhase] = 0.5*(this->curElemDat_.vertex[i].density[wPhase] + this->curElemDat_.vertex[j].density[wPhase]);
                avgDensity[nPhase] = 0.5*(this->curElemDat_.vertex[i].density[nPhase] + this->curElemDat_.vertex[j].density[nPhase]);

                diffusionAW = avgDpm[wPhase] * avgDensity[wPhase] * normDiffGrad[wPhase];
                diffusionWW = - diffusionAW;
                diffusionWN = avgDpm[nPhase] * avgDensity[nPhase] * normDiffGrad[nPhase];
                diffusionAN = - diffusionWN;

                // add diffusion of water to water flux
                flux[wComp] += (diffusionWW + diffusionWN);

                // add diffusion of air to air flux
                flux[nComp] += (diffusionAN + diffusionAW);
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
                    int globalIdx = this->problem_.vertIdx(*it);
                    const GlobalPosition &globalPos = it->geometry()[0];

                    // initialize phase state
                    staticVertexDat_[globalIdx].phaseState =
                        this->problem_.initialPhaseState(*it, globalIdx, globalPos);
                    staticVertexDat_[globalIdx].oldPhaseState =
                        staticVertexDat_[globalIdx].phaseState;
                }
            }

        /*!
         * \brief Update the static data of a single vert and do a
         *        variable switch if necessary.
         */
        void updateStaticData(SpatialFunction &curGlobalSol, SpatialFunction &oldGlobalSol)
            {
                bool wasSwitched = false;

                VertexIterator it = this->problem_.vertexBegin();
                for (; it != this->problem_.vertexEnd(); ++it)
                {
                    int globalIdx = this->problem_.vertIdx(*it);
                    const GlobalPosition &global = it->geometry()[0];

                    wasSwitched = primaryVarSwitch_(curGlobalSol,
                                                    globalIdx,
                                                    global)
                        || wasSwitched;
                }

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
         * \brief Reset the current phase state of all verts to the old one after an update failed
         */
        void resetPhaseState()
            {
                int numVertices = this->problem_.numVertices();
                for (int i = 0; i < numVertices; ++i)
                    staticVertexDat_[i].phaseState = staticVertexDat_[i].oldPhaseState;
            }

        /*!
         * \brief Return true if the primary variables were switched
         *        after the last timestep.
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
         * \brief Add the mass fraction of air in water to VTK output of
         *        the current timestep.
         */
        template <class MultiWriter>
        void addVtkFields(MultiWriter &writer, const SpatialFunction &globalSol) const
            {
                typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

                // create the required scalar fields
                unsigned numVertices = this->problem_.numVertices();
                ScalarField *pW =           writer.template createField<Scalar, 1>(numVertices);
                ScalarField *pN =           writer.template createField<Scalar, 1>(numVertices);
                ScalarField *pC =           writer.template createField<Scalar, 1>(numVertices);
                ScalarField *Sw =           writer.template createField<Scalar, 1>(numVertices);
                ScalarField *Sn =           writer.template createField<Scalar, 1>(numVertices);
                ScalarField *mobW =         writer.template createField<Scalar, 1>(numVertices);
                ScalarField *mobN =         writer.template createField<Scalar, 1>(numVertices);
                ScalarField *massfracAinW = writer.template createField<Scalar, 1>(numVertices);
                ScalarField *massfracAinA = writer.template createField<Scalar, 1>(numVertices);
                ScalarField *massfracWinW = writer.template createField<Scalar, 1>(numVertices);
                ScalarField *massfracWinA = writer.template createField<Scalar, 1>(numVertices);
                ScalarField *temperature  = writer.template createField<Scalar, 1>(numVertices);
                ScalarField *phaseState   = writer.template createField<Scalar, 1>(numVertices);

                VariableVertexData tmp;
                ElementIterator it = this->problem_.elementBegin();
                ElementIterator endit = this->problem_.elementEnd();
                for (; it != endit; ++it) {
                    for (int i = 0; i < it->template count<dim>(); ++i) {
                        int globalI = this->problem_.vertIdx(*it, i);
                        asImp_()->updateVarVertexData_(tmp,
                                                     (*globalSol)[globalI],
                                                     staticVertexDat_[globalI].phaseState,
                                                     *it,
                                                     i,
                                                     this->problem_,
                                                     Implementation::temperature_((*globalSol)[globalI]));

                        (*pW)[globalI] = tmp.pW;
                        (*pN)[globalI] = tmp.pN;
                        (*pC)[globalI] = tmp.pC;
                        (*Sw)[globalI] = tmp.satW;
                        (*Sn)[globalI] = tmp.satN;
                        (*mobW)[globalI] = tmp.mobility[wPhase];
                        (*mobN)[globalI] = tmp.mobility[nPhase];
                        (*massfracAinW)[globalI] = tmp.massfrac[nComp][wPhase];
                        (*massfracAinA)[globalI] = tmp.massfrac[nComp][nPhase];
                        (*massfracWinW)[globalI] = tmp.massfrac[wComp][wPhase];
                        (*massfracWinA)[globalI] = tmp.massfrac[wComp][nPhase];
                        (*temperature)[globalI] = Implementation::temperature_((*globalSol)[globalI]);
                        (*phaseState)[globalI] = staticVertexDat_[globalI].phaseState;
                    };
                }


                writer.addVertexData(pW, "pW");
                writer.addVertexData(pN, "pN");
                writer.addVertexData(pC, "pC");
                writer.addVertexData(Sw, "Sw");
                writer.addVertexData(Sn, "Sn");
                writer.addVertexData(mobW, "mobW");
                writer.addVertexData(mobN, "mobN");
                writer.addVertexData(massfracAinW, "Xaw");
                writer.addVertexData(massfracAinA, "Xaa");
                writer.addVertexData(massfracWinW, "Xww");
                writer.addVertexData(massfracWinA, "Xwa");
                writer.addVertexData(temperature, "T");
                writer.addVertexData(phaseState, "phase state");
            }


    protected:
        Implementation *asImp_()
        { return static_cast<Implementation *>(this); }
        const Implementation *asImp_() const
        { return static_cast<const Implementation *>(this); }


        void updateElementData_(ElementData &dest, const LocalFunction &sol, bool isOldSol)
            {
                int phaseState;
                int numVertices = this->curElement_().template count<dim>();
                for (int i = 0; i < numVertices; i++) {
                    int iGlobal = ParentType::problem_.vertIdx(ParentType::curElement_(), i);
                    phaseState = isOldSol?staticVertexDat_[iGlobal].oldPhaseState:staticVertexDat_[iGlobal].phaseState;
                    asImp_()->updateVarVertexData_(dest.vertex[i],
                                                 sol[i],
                                                 phaseState,
                                                 this->curElement_(),
                                                 i,
                                                 this->problem_,
                                                 Implementation::temperature_(sol[i]));
                }
            }


        //  perform variable switch at a vert. Retrurns true iff a
        //  variable switch was performed.
        bool primaryVarSwitch_(SpatialFunction &globalSol,
                               int globalIdx,
                               const GlobalPosition &globalPos)
            {
                // evaluate primary variable switch
                int phaseState    = staticVertexDat_[globalIdx].phaseState;
                int newPhaseState = phaseState;

                // Evaluate saturation and pressures
                Scalar pW = (*globalSol)[globalIdx][pWIdx];
                Scalar satW = 0.0;
                Scalar temperature = Implementation::temperature_((*globalSol)[globalIdx]);
                if      (phaseState == bothPhases) satW = 1.0 - (*globalSol)[globalIdx][switchIdx];
                else if (phaseState == wPhaseOnly) satW = 1.0;
                else if (phaseState == nPhaseOnly) satW = 0.0;

                Scalar pC = this->problem_.pC(satW, globalIdx, globalPos);
                Scalar pN = pW + pC;

                if (phaseState == nPhaseOnly)
                {
                    Scalar xWN = (*globalSol)[globalIdx][switchIdx];
                    Scalar xWNmax = this->problem_.multicomp().xWN(pN, temperature);

                    if (xWN > xWNmax*(1 + 2e-5))
                    {
                        // wetting phase appears
                        std::cout << "wetting phase appears at vert " << globalIdx
                                  << ", coordinates: " << globalPos << std::endl;
                        newPhaseState = bothPhases;
                        (*globalSol)[globalIdx][switchIdx] = 1.0 - 2e-5;
                    };
                }
                else if (phaseState == wPhaseOnly)
                {
                    Scalar xAW = (*globalSol)[globalIdx][switchIdx];
                    Scalar xAWmax = this->problem_.multicomp().xAW(pN, temperature);

                    if (xAW > xAWmax*(1 + 2e-5))
                    {
                        // non-wetting phase appears
                        std::cout << "Non-wetting phase appears at vert " << globalIdx
                                  << ", coordinates: " << globalPos << std::endl;
                        (*globalSol)[globalIdx][switchIdx] = 2e-5;
                        newPhaseState = bothPhases;
                    }
                }
                else if (phaseState == bothPhases) {
                    Scalar satN = 1 - satW;

                    if (satN < -1e-5) {
                        // non-wetting phase disappears
                        std::cout << "Non-wetting phase disappears at vert " << globalIdx
                                  << ", coordinates: " << globalPos << std::endl;
                        (*globalSol)[globalIdx][switchIdx]
                            = this->problem_.multicomp().xAW(pN, temperature);
                        newPhaseState = wPhaseOnly;
                    }
                    else if (satW < -1e-5) {
                        // wetting phase disappears
                        std::cout << "Wetting phase disappears at vert " << globalIdx
                                  << ", coordinates: " << globalPos << std::endl;
                        (*globalSol)[globalIdx][switchIdx]
                            = this->problem_.multicomp().xWN(pN, temperature);
                        newPhaseState = nPhaseOnly;
                    }
                }

                staticVertexDat_[globalIdx].phaseState = newPhaseState;

                return phaseState != newPhaseState;
            }

        // harmonic mean of the permeability computed directly.  the
        // first parameter is used to store the result.
        static void harmonicMeanK_(Tensor &Ki, const Tensor &Kj)
            {
                double eps = 1e-20;

                for (int kx=0; kx < Tensor::rows; kx++){
                    for (int ky=0; ky< Tensor::cols; ky++){
                        if (Ki[kx][ky] != Kj[kx][ky]) {
                            Ki[kx][ky] = 2 / (1/(Ki[kx][ky]+eps) + (1/(Kj[kx][ky]+eps)));
                        }
                    }
                }
            }

        // parameters given in constructor
        std::vector<StaticVertexData> staticVertexDat_;
        bool                        switchFlag_;

        // current solution
        LocalFunction       curSol_;
        ElementData        curElemDat_;

        // needed for restoreCurSolution()
        bool               curSolDeflected_;
        Scalar                 curSolOrigValue_;
        VariableVertexData curSolOrigVarData_;

        // previous solution
        LocalFunction      prevSol_;
        ElementData        prevElemDat_;
    };


    template<class ProblemT,
             class BoxTraitsT,
             class TwoPTwoCTraitsT>
    class TwoPTwoCBoxJacobian : public TwoPTwoCBoxJacobianBase<ProblemT,
                                                               BoxTraitsT,
                                                               TwoPTwoCTraitsT,
                                                               TwoPTwoCBoxJacobian<ProblemT,
                                                                                   BoxTraitsT,
                                                                                   TwoPTwoCTraitsT>
                                                               >
    {
        typedef TwoPTwoCBoxJacobian<ProblemT,
                                      BoxTraitsT,
                                      TwoPTwoCTraitsT> ThisType;
        typedef TwoPTwoCBoxJacobianBase<ProblemT,
                                        BoxTraitsT,
                                        TwoPTwoCTraitsT,
                                        ThisType>      ParentType;

        typedef ProblemT                               Problem;

        typedef typename Problem::DomainTraits         DomTraits;
        typedef TwoPTwoCTraitsT                        TwoPTwoCTraits;
        typedef BoxTraitsT                             BoxTraits;

        typedef typename DomTraits::Scalar             Scalar;

        typedef typename DomTraits::GlobalPosition         GlobalPosition;
        typedef typename DomTraits::LocalPosition         LocalPosition;

        typedef typename BoxTraits::SolutionVector     SolutionVector;
        typedef typename BoxTraits::FVElementGeometry  FVElementGeometry;

        typedef typename ParentType::VariableVertexData  VariableVertexData;
        typedef typename ParentType::LocalFunction     LocalFunction;
        typedef typename ParentType::ElementData         ElementData;

        typedef typename TwoPTwoCTraits::PhasesVector PhasesVector;


    public:
        TwoPTwoCBoxJacobian(ProblemT &problem)
            : ParentType(problem)
            {
            };


        /*!
         * \brief The storage term of heat
         */
        void heatStorage(SolutionVector &result,
                         int scvId,
                         const LocalFunction &sol,
                         const ElementData &elementCache) const
            {
                // only relevant for the non-isothermal model!
            }

        /*!
         * \brief Update the temperature gradient at a face of a FV
         *        element.
         */
        void updateTempGrad(GlobalPosition &tempGrad,
                            const GlobalPosition &feGrad,
                            const LocalFunction &sol,
                            int vertIdx) const
            {
                // only relevant for the non-isothermal model!
            }

        /*!
         * \brief Sets the temperature term of the flux vector to the
         *        heat flux due to advection of the fluids.
         */
        void advectiveHeatFlux(SolutionVector &flux,
                               const PhasesVector &darcyOut,
                               Scalar alpha, // upwind parameter
                               const VariableVertexData *upW, // up/downstream verts
                               const VariableVertexData *dnW,
                               const VariableVertexData *upN,
                               const VariableVertexData *dnN) const
            {
                // only relevant for the non-isothermal model!
            }

        /*!
         * \brief Adds the diffusive heat flux to the flux vector over
         *        the face of a sub-control volume.
         */
        void diffusiveHeatFlux(SolutionVector &flux,
                               int faceIdx,
                               const GlobalPosition &tempGrad) const
            {
                // only relevant for the non-isothermal model!
            }

        // internal method!
        static Scalar temperature_(const SolutionVector &sol)
            { return 283.15; }
    };

    ///////////////////////////////////////////////////////////////////////////
    // TwoPTwoCBoxModel (The actual numerical model.)
    ///////////////////////////////////////////////////////////////////////////
    /**
     * \brief Isothermal two phase two component model with Pw and
     *        Sn/X as primary unknowns.
     *
     * This implements an isothermal two phase two component model
     * with Pw and Sn/X as primary unknowns.
     */
    template<class ProblemT>
    class TwoPTwoCBoxModel
        : public BoxScheme<TwoPTwoCBoxModel<ProblemT>, // Implementation of the box scheme

                           // The Traits for the BOX method
                           P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                       typename ProblemT::DomainTraits::Grid,
                                       TwoPTwoCTraits<typename ProblemT::DomainTraits::Scalar>::numEq>,

                           // The actual problem we would like to solve
                           ProblemT,

                           // The local jacobian operator
                           TwoPTwoCBoxJacobian<ProblemT,
                                               P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                                           typename ProblemT::DomainTraits::Grid,
                                                           TwoPTwoCTraits<typename ProblemT::DomainTraits::Scalar>::numEq>,
                                               TwoPTwoCTraits<typename ProblemT::DomainTraits::Scalar> > >
    {
        typedef typename ProblemT::DomainTraits::Grid   Grid;
        typedef typename ProblemT::DomainTraits::Scalar Scalar;
        typedef TwoPTwoCBoxModel<ProblemT>              ThisType;

    public:
        typedef Dune::TwoPTwoCTraits<Scalar>                           TwoPTwoCTraits;
        typedef P1BoxTraits<Scalar, Grid, TwoPTwoCTraits::numEq> BoxTraits;

    private:
        typedef TwoPTwoCBoxJacobian<ProblemT, BoxTraits, TwoPTwoCTraits>  TwoPTwoCLocalJacobian;
        typedef BoxScheme<ThisType,
                          BoxTraits,
                          ProblemT,
                          TwoPTwoCLocalJacobian>        ParentType;

        typedef typename ProblemT::DomainTraits           DomTraits;
        typedef typename DomTraits::Element                  Element;
        typedef typename DomTraits::ElementIterator          ElementIterator;
        typedef typename DomTraits::LocalPosition            LocalPosition;
        typedef typename DomTraits::GlobalPosition            GlobalPosition;

        enum {
            dim              = DomTraits::dim,
            dimWorld         = DomTraits::dimWorld
        };

    public:
        typedef NewNewtonMethod<ThisType> NewtonMethod;

        TwoPTwoCBoxModel(ProblemT &prob)
            : ParentType(prob, twoPTwoCLocalJacobian_),
              twoPTwoCLocalJacobian_(prob)
            {
                Api::require<Api::BasicDomainTraits, typename ProblemT::DomainTraits>();
            }


        /*!
         * \brief Called by the update() method if applying the newton
         *         method was unsuccessful.
         */
        void updateFailedTry()
            {
                ParentType::updateFailedTry();

                twoPTwoCLocalJacobian_.setSwitched(false);
                twoPTwoCLocalJacobian_.resetPhaseState();
                twoPTwoCLocalJacobian_.updateStaticData(this->currentSolution(),
                                                        this->previousSolution());
            };

        /*!
         * \brief Called by the BoxScheme's update method.
         */
        void updateSuccessful()
            {
                ParentType::updateSuccessful();

                twoPTwoCLocalJacobian_.updateOldPhaseState();
                twoPTwoCLocalJacobian_.setSwitched(false);
            }


        /*!
         * \brief Add the mass fraction of air in water to VTK output of
         *        the current timestep.
         */
        template <class MultiWriter>
        void addVtkFields(MultiWriter &writer) const
            {
                twoPTwoCLocalJacobian_.addVtkFields(writer, this->currentSolution());
            }

        /*!
         * \brief Returns true if there was a primary variable switch
         *        after the last time step.
         */
        bool switched() const
            { return twoPTwoCLocalJacobian_.switched(); }


    private:
        // calculates the jacobian matrix at a given position
        TwoPTwoCLocalJacobian  twoPTwoCLocalJacobian_;
    };
}

#endif
