// $Id:$
/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
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
#ifndef DUMUX_TWOP_BOX_JACOBIAN_BASE_HH
#define DUMUX_TWOP_BOX_JACOBIAN_BASE_HH

#include <dumux/boxmodels/boxscheme/boxscheme.hh>

#include "2pproperties.hh"

#include "2pvertexdata.hh"
#include "2pelementdata.hh"
#include "2pfluxdata.hh"

#include <dune/common/collectivecommunication.hh>
#include <vector>
#include <iostream>

namespace Dune
{
/*!
 * \ingroup TwoPBoxModel
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase box model.
 *
 * This class is also used for the non-isothermal model, which means
 * that it uses static polymorphism.
 */
template<class TypeTag, class Implementation>
class TwoPBoxJacobianBase : public BoxJacobian<TypeTag, Implementation>
{
protected:
    typedef TwoPBoxJacobianBase<TypeTag, Implementation> ThisType;
    typedef BoxJacobian<TypeTag, Implementation>         ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))      Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices))  Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))       Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))     GridView;

    enum {
        dim            = GridView::dimension,
        dimWorld       = GridView::dimensionworld,

        numEq          = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        numPhases      = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),

        pressureIdx    = Indices::pressureIdx,
        saturationIdx  = Indices::saturationIdx,

        wPhase         = Indices::wPhase,
        nPhase         = Indices::nPhase,
    };

    typedef typename GridView::template Codim<0>::Entity   Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef FieldVector<Scalar, dim>       LocalPosition;
    typedef FieldVector<Scalar, dimWorld>  GlobalPosition;

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVarVector        PrimaryVarVector;
    typedef typename SolutionTypes::SolutionFunction        SolutionFunction;
    typedef typename SolutionTypes::DofEntityMapper         DofEntityMapper;
    typedef typename SolutionTypes::SolutionOnElement       SolutionOnElement;

    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexData))   VertexData;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementData))  ElementData;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxData))     FluxData;

    typedef std::vector<VertexData>        VertexDataArray;
    typedef FieldMatrix<Scalar, dim, dim>  Tensor;

    static const Scalar mobilityUpwindAlpha = GET_PROP_VALUE(TypeTag, PTAG(MobilityUpwindAlpha));

public:
    TwoPBoxJacobianBase(Problem &problem)
        : ParentType(problem)
    {
    };

    /*!
     * \brief Evaluate the amount all conservation quantites
     *        (e.g. phase mass) within a finite sub-control volume.
     */
    void computeStorage(PrimaryVarVector &result, int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const VertexDataArray &elemDat = usePrevSol ? this->prevElemDat_  : this->curElemDat_;
        const VertexData  &vertDat = elemDat[scvIdx];

        // wetting phase mass
        result[Indices::phase2Mass(wPhase)] =
            vertDat.density[wPhase]
            * vertDat.porosity
            * vertDat.satW;

        // non-wetting phase mass
        result[Indices::phase2Mass(nPhase)] =
            vertDat.density[nPhase]
            * vertDat.porosity
            * vertDat.satN;
    }

    /*!
     * \brief Evaluates the mass flux over a face of a sub-control
     *        volume.
     */
    void computeFlux(PrimaryVarVector &flux, int faceId) const
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
     *
     * This method is called by compute flux and is mainly there for
     * derived models to ease adding equations selectively.
     */
    void computeAdvectiveFlux(PrimaryVarVector &flux,
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
            flux[Indices::phase2Mass(phaseIdx)] +=
                vars.vDarcyNormal[phaseIdx] * (
                    mobilityUpwindAlpha* // upstream vertex
                    (  up.density[phaseIdx] *
                       up.mobility[phaseIdx])
                    +
                    (1 - mobilityUpwindAlpha)* // downstream vertex
                    (  dn.density[phaseIdx] *
                       dn.mobility[phaseIdx]));
        }
    }

    /*!
     * \brief Adds the diffusive flux to the flux vector over
     *        the face of a sub-control volume.
     
     * It doesn't do anything in two-phase model but is used by the
     * non-isothermal two-phase models to calculate diffusive heat
     * fluxes
     */
    void computeDiffusiveFlux(PrimaryVarVector &flux,
                              const FluxData &fluxData) const
    {
        // diffusive fluxes
        //flux += 0.0;
    }

    /*!
     * \brief Calculate the source term of the equation
     */
    void computeSource(PrimaryVarVector &q, int localVertexIdx)
    {
        // retrieve the source term intrinsic to the problem
        this->problem_.source(q,
                              this->curElement_(),
                              this->curElementGeom_,
                              localVertexIdx);
    }

    /*!
     * \brief Return the temperature given the solution vector of a
     *        finite volume.
     */
    template <class PrimaryVarVector>
    Scalar temperature(const PrimaryVarVector &sol)
    {
        return this->problem_.temperature(); /* constant temperature */
    }

    /*!
     * \brief Calculate mass of both components in the whole model domain
     *         and get minimum and maximum values of primary variables
     *
     */
    void calculateMass(const SolutionFunction &globalSol, Dune::FieldVector<Scalar, 2> &mass)
    {
        ElementIterator elementIt = this->problem_.elementBegin();
        ElementIterator endit = this->problem_.elementEnd();
        unsigned numVertices = this->problem_.numVertices();
        SolutionOnElement tmpSol;
        VertexDataArray elemDat;
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

            // Loop over element vertices
            int numLocalVerts = elementIt->template count<dim>();
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
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     */
    template <class MultiWriter>
    void addVtkFields(MultiWriter &writer, const SolutionFunction &globalSol)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->problem_.gridView().size(dim);
        ScalarField *pW = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *pN = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *pC = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *Sw = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *Sn = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *Te = writer.template createField<Scalar, 1>(numVertices);

        const DofEntityMapper &dofMapper = this->problem_.model().dofEntityMapper();
        SolutionOnElement tmpSol;
        ElementIterator elementIt = this->problem_.gridView().template begin<0>();
        ElementIterator endit = this->problem_.gridView().template end<0>();
        for (; elementIt != endit; ++elementIt)
        {
            setCurrentElement(*elementIt);
            this->restrictToElement(tmpSol, globalSol);
            this->setCurrentSolution(tmpSol);

            int numVerts = elementIt->template count<dim>();
            for (int i = 0; i < numVerts; ++i)
            {
                int globalIdx = dofMapper.map(*elementIt, 
                                              i,
                                              dim);
                
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

/*!
 * \brief Calculate the local Jacobian for the two phase model in the
 *        BOX scheme.
 *
 * This is just a wrapper around TwoPBoxJacobianBase.
 */
template<class TypeTag>
class TwoPBoxJacobian : public TwoPBoxJacobianBase<TypeTag,
                                                   TwoPBoxJacobian<TypeTag> >
{
    typedef TwoPBoxJacobian<TypeTag>                       ThisType;
    typedef TwoPBoxJacobianBase<TypeTag, ThisType>         ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

public:
    TwoPBoxJacobian(Problem &problem)
        : ParentType(problem)
    {
    };
};

}

#endif
