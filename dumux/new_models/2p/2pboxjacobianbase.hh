//$Id:$
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
#ifndef DUMUX_TWOP_BOX_JACOBIAN_BASE_HH
#define DUMUX_TWOP_BOX_JACOBIAN_BASE_HH

#include <dumux/new_models/boxscheme/boxscheme.hh>
#include <dumux/new_models/boxscheme/p1boxtraits.hh>
#include <dumux/new_models/2p/2ptraits.hh>
#include <dumux/auxiliary/math.hh>

#include <dumux/new_models/2p/2pvertexdata.hh>
#include <dumux/new_models/2p/2pfluxdata.hh>

#include <dumux/auxiliary/apis.hh>
#include <dune/common/collectivecommunication.hh>
#include <vector>
#include <iostream>

namespace Dune
{
///////////////////////////////////////////////////////////////////////////
// TwoPBoxJacobian (evaluate the local jacobian for the newton method.)
///////////////////////////////////////////////////////////////////////////
/*!
 * \brief Calculate the local Jacobian for the two phase model in the
 *        BOX scheme.
 */
template<class ProblemT,
		 class BoxTraitsT,
		 class TwoPTraitsT,
		 class VertexDataT,
		 class FluxDataT,
		 class Implementation>
class TwoPBoxJacobianBase : public BoxJacobian<ProblemT,
                                           BoxTraitsT,
                                           Implementation,
                                           VertexDataT>
{
    protected:
    typedef TwoPBoxJacobianBase<ProblemT,
                            BoxTraitsT,
                            TwoPTraitsT,
                            VertexDataT,
                            FluxDataT,
                            Implementation>            ThisType;

    typedef BoxJacobian<ProblemT,
                        BoxTraitsT,
                        Implementation,
                        VertexDataT>   ParentType;

    typedef ProblemT                       Problem;
    typedef typename Problem::DomainTraits DomTraits;
    typedef BoxTraitsT                     BoxTraits;
    typedef TwoPTraitsT                    TwoPTraits;

    enum {
        dim            = DomTraits::dim,
        dimWorld       = DomTraits::dimWorld,

        numEq          = BoxTraits::numEq,
        numPhases      = TwoPTraits::numPhases,

        pressureIdx    = TwoPTraits::pressureIdx,
        saturationIdx  = TwoPTraits::saturationIdx,

        wMassIdx       = TwoPTraits::wMassIdx,
        nMassIdx       = TwoPTraits::nMassIdx,

        formulation    = TwoPTraits::formulation,

        wPhase         = TwoPTraits::wPhase,
        nPhase         = TwoPTraits::nPhase,
    };

    typedef typename DomTraits::Scalar              Scalar;
    typedef typename DomTraits::CoordScalar         CoordScalar;
    typedef typename DomTraits::Grid                Grid;
    typedef typename DomTraits::Element             Element;
    typedef typename DomTraits::ElementIterator     ElementIterator;
    typedef typename Element::EntityPointer         ElementPointer;

    typedef typename DomTraits::LocalPosition       LocalPosition;
    typedef typename DomTraits::GlobalPosition      GlobalPosition;

    typedef typename BoxTraits::SolutionVector      SolutionVector;
    typedef typename BoxTraits::FVElementGeometry   FVElementGeometry;
    typedef typename BoxTraits::SpatialFunction     SpatialFunction;
    typedef typename BoxTraits::LocalFunction       LocalFunction;
    typedef typename Grid::CollectiveCommunication  CollectiveCommunication;

    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;
    typedef VertexDataT          VertexData;
    typedef FluxDataT            FluxData;

    typedef std::vector<VertexData>        VertexDataArray;
    typedef FieldMatrix<Scalar, dim, dim>  Tensor;
    static const Scalar upwindAlpha = TwoPTraits::upwindAlpha;

public:
    TwoPBoxJacobianBase(ProblemT &problem)
        : ParentType(problem)
    {};

    /*!
     * \brief Evaluate the amount all conservation quantites
     *        (e.g. phase mass) within a finite volume.
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

        // wetting phase mass
        result[wMassIdx] =
            vertDat.density[wPhase]
            * vertDat.porosity
            * vertDat.satW;

        // non-wetting phase mass
        result[nMassIdx] =
            vertDat.density[nPhase]
            * vertDat.porosity
            * vertDat.satN;
    }

    /*!
     * \brief Evaluates the mass flux over a face of a subcontrol
     *        volume.
     */
    void computeFlux(SolutionVector &flux, int faceId) const
    {
		FluxData vars(this->problem_,
                      this->curElement_(),
                      this->curElementGeom_,
                      faceId,
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

                // add advective flux of current component in current
                // phase
                flux[phaseIdx] +=
                    vars.vDarcyNormal[phaseIdx] * (
                        upwindAlpha* // upstream vertex
                        (  up.density[phaseIdx] *
                           up.mobility[phaseIdx])
                        +
                        (1 - upwindAlpha)* // downstream vertex
                        (  dn.density[phaseIdx] *
                           dn.mobility[phaseIdx]));
        }
    }

    /*!
     * \brief Adds the diffusive flux to the flux vector over
     *        the face of a sub-control volume.
     *        Is not needed in two phase model but needs to
     *        be used in the non-isothermal two-phase model
     *        to calculate diffusive heat fluxes
     */
    void computeDiffusiveFlux(SolutionVector &flux,
                              const FluxData &fluxData) const
    {
        // diffusive fluxes
        flux += 0.0;
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
    {
        return this->problem_.temperature(); /* constant temperature */
    }

    /*!
     * \brief Calculate mass of both components in the whole model domain
     *         and get minimum and maximum values of primary variables
     *
     */
    void calculateMass(const SpatialFunction &globalSol, Dune::FieldVector<Scalar, 2> &mass)
    {
        ElementIterator elementIt = this->problem_.elementBegin();
        ElementIterator endit = this->problem_.elementEnd();
        unsigned numVertices = this->problem_.numVertices();
        LocalFunction tmpSol;
        VertexDataArray elemDat(BoxTraits::ShapeFunctionSetContainer::maxsize);
        VertexData tmp;
        Scalar vol, poro, rhoN, rhoW, satN, satW, pW, Te;
        Scalar massNPhase(0.), massWPhase(0.);

        mass = 0;
        Scalar minSat = 1e100;
        Scalar maxSat = -1e100;
        Scalar minP = 1e100;
        Scalar maxP = -1e100;
        Scalar minTe = 1e100;
        Scalar maxTe = -1e100;

        // Loop over elements
        for (; elementIt != endit; ++elementIt)
        {
            setCurrentElement(*elementIt);
            this->restrictToElement(tmpSol, globalSol);
            this->updateElementData_(elemDat, tmpSol, false);
            // get geometry type

            int numLocalVerts = elementIt->template count<dim>();
            tmpSol.resize(numLocalVerts);
            // Loop over element vertices
            for (int i = 0; i < numLocalVerts; ++i)
            {
                int globalIdx = this->problem_.vertexIdx(*elementIt, i);

                vol = this->curElementGeom_.subContVol[i].volume;
                poro = this->problem_.porosity(this->curElement_(), i);
                rhoN = elemDat[i].density[nPhase];
                rhoW = elemDat[i].density[wPhase];
                satN = elemDat[i].saturation[nPhase];
                satW = elemDat[i].saturation[wPhase];
                pW = elemDat[i].pressure[wPhase];
                Te = Implementation::temperature_((*globalSol)[globalIdx]);

                massNPhase = vol * poro * satN * rhoN;
                massWPhase = vol * poro * satW * rhoW;

                // get minimum and maximum values of primary variables
                minSat = std::min(minSat, satN);
                maxSat = std::max(maxSat, satN);
                minP = std::min(minP, pW);
                maxP = std::max(maxP, pW);
                minTe = std::min(minTe, Te);
                maxTe = std::max(maxTe, Te);

                // calculate total mass
                mass[0] += massNPhase; // mass nonwetting phase
                mass[1] += massWPhase; // mass in wetting phase
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
        ScalarField *pW =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *pN =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *pC =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *Sw =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *Sn =           writer.template createField<Scalar, 1>(numVertices);
        ScalarField *Te =           writer.template createField<Scalar, 1>(numVertices);


        LocalFunction tmpSol;
        ElementIterator elementIt = this->problem_.elementBegin();
        ElementIterator endit = this->problem_.elementEnd();
        for (; elementIt != endit; ++elementIt)
        {
            int numLocalVerts = elementIt->template count<dim>();
            tmpSol.resize(numLocalVerts);

            setCurrentElement(*elementIt);
            this->restrictToElement(tmpSol, globalSol);
            this->setCurrentSolution(tmpSol);

            for (int i = 0; i < numLocalVerts; ++i)
            {
                int globalIdx = this->problem_.vertexIdx(*elementIt, i);

                (*pW)[globalIdx] = this->curElemDat_[i].pressure[wPhase];
                (*pN)[globalIdx] = this->curElemDat_[i].pressure[nPhase];
                (*pC)[globalIdx] = this->curElemDat_[i].pC;
                (*Sw)[globalIdx] = this->curElemDat_[i].satW;
                (*Sn)[globalIdx] = this->curElemDat_[i].satN;
                (*Te)[globalIdx] = asImp_()->temperature((*globalSol)[globalIdx]);
            };
        }

        writer.addVertexData(pW, "pW");
        writer.addVertexData(pN, "pN");
        writer.addVertexData(pC, "pC");
        writer.addVertexData(Sw, "SW");
        writer.addVertexData(Sn, "SN");
        writer.addVertexData(Te, "Te");
    }
    protected:
    Implementation *asImp_()
    { return static_cast<Implementation *>(this); }
    const Implementation *asImp_() const
    { return static_cast<const Implementation *>(this); }
};
}

#endif
