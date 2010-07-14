// $Id: 2pboxjacobian.hh 3794 2010-06-25 16:04:52Z melanie $
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
#include "2pfluidstate.hh"

#include <dune/common/collectivecommunication.hh>
#include <vector>
#include <iostream>

namespace Dumux
{
/*!
 * \ingroup TwoPBoxModel
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the two-phase box model.
 *
 * This class is also used for the non-isothermal model, which means
 * that it uses static polymorphism.
 */
template<class TypeTag>
class TwoPBoxJacobian: public BoxJacobian<TypeTag>
{
protected:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(LocalJacobian)) Implementation;
    typedef TwoPBoxJacobian<TypeTag> ThisType;
    typedef BoxJacobian<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    enum
    {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),

        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,

        contiWEqIdx = Indices::contiWEqIdx,
        contiNEqIdx = Indices::contiNEqIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVarVector PrimaryVarVector;
    typedef typename SolutionTypes::SolutionVector SolutionVector;
    typedef typename SolutionTypes::DofEntityMapper DofEntityMapper;
    typedef typename SolutionTypes::SolutionOnElement SolutionOnElement;

    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexData)) VertexData;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementData)) ElementData;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxData)) FluxData;

    typedef std::vector<VertexData> VertexDataArray;
    typedef Dune::FieldMatrix<Scalar, dim, dim> Tensor;

    static const Scalar mobilityUpwindAlpha =
            GET_PROP_VALUE(TypeTag, PTAG(MobilityUpwindAlpha));

public:
    TwoPBoxJacobian(Problem &problem) :
        ParentType(problem)
    {
    }
    ;

    /*!
     * \internal
     * \brief Returns the epsilon value which is added and removed
     *        from the current solution.
     *
     * \param elemSol    The current solution on the element
     * \param vertexIdx  The local index of the element's vertex for
     *                   which the local derivative ought to be calculated.
     * \param pvIdx      The index of the primary variable which gets varied
     */
    Scalar numericEpsilon_(const SolutionOnElement &elemSol, int vertIdx,
            int pvIdx) const
    {
        if (pvIdx == 0)
            return 1e-2; // pressure
        else if (pvIdx == 1)
            return 1e-7; // saturation

        return ParentType::numericEpsilon_(elemSol, vertIdx, pvIdx);
    }

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
        const VertexDataArray &elemDat = usePrevSol ? this->prevElemDat_
                : this->curElemDat_;
        const VertexData &vertDat = elemDat[scvIdx];

        // wetting phase mass
        result[contiWEqIdx] = vertDat.density(wPhaseIdx) * vertDat.porosity()
                * vertDat.saturation(wPhaseIdx);

        // non-wetting phase mass
        result[contiNEqIdx] = vertDat.density(nPhaseIdx) * vertDat.porosity()
                * vertDat.saturation(nPhaseIdx);
        ;
    }

    /*!
     * \brief Evaluates the mass flux over a face of a sub-control
     *        volume.
     */
    void computeFlux(PrimaryVarVector &flux, int faceId) const
    {
        FluxData vars(this->problem_, this->curElement_(),
                this->curElementGeom_, faceId, this->curElemDat_);
        flux = 0;
        asImp_()->computeAdvectiveFlux(flux, vars);
        asImp_()->computeDiffusiveFlux(flux, vars);
        flux *= -1;
    }

    /*!
     * \brief Evaluates the advective mass flux of all components over
     *        a face of a subcontrol volume.
     *
     * This method is called by compute flux and is mainly there for
     * derived models to ease adding equations selectively.
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

            // add advective flux of current component in current
            // phase
            int eqIdx = (phaseIdx == wPhaseIdx) ? contiWEqIdx : contiNEqIdx;
            flux[eqIdx] += vars.KmvpNormal(phaseIdx) * (mobilityUpwindAlpha * // upstream vertex
                    (up.density(phaseIdx) * up.mobility(phaseIdx)) + (1
                    - mobilityUpwindAlpha) * // downstream vertex
                    (dn.density(phaseIdx) * dn.mobility(phaseIdx)));
        }
    }

    /*!
     * \brief Adds the diffusive flux to the flux vector over
     *        the face of a sub-control volume.

     * It doesn't do anything in two-phase model but is used by the
     * non-isothermal two-phase models to calculate diffusive heat
     * fluxes
     */
    void computeDiffusiveFlux(PrimaryVarVector &flux, const FluxData &fluxData) const
    {
        // diffusive fluxes
        flux += 0.0;
    }

    /*!
     * \brief Calculate the source term of the equation
     */
    void computeSource(PrimaryVarVector &q, int localVertexIdx)
    {
        // retrieve the source term intrinsic to the problem
        this->problem_.source(q, this->curElement_(), this->curElementGeom_,
                localVertexIdx);
    }

    /*!
     * \brief Calculate the fluid phases flux across a certain layer in the domain.
     * The layer is situated perpendicular to the coordinate axis "coord" and cuts
     * the axis at the value "coordValue"
     *
     */
    void calculateFluxAcrossLayer(const SolutionVector &globalSol,
            Dune::FieldVector<Scalar, 2> &flux, int coord, Scalar coordVal)
    {
        SolutionOnElement tmpSol;
        ElementIterator elementIt =
                this->problem_.gridView().template begin<0> ();
        ElementIterator endit = this->problem_.gridView().template end<0> ();
        GlobalPosition globalI, globalJ;
        PrimaryVarVector tmpFlux(0.0);
        int sign;

        // Loop over elements
        for (; elementIt != endit; ++elementIt)
        {
            if (elementIt->partitionType() != Dune::InteriorEntity)
                continue;

            setCurrentElement(*elementIt);
            this->restrictToElement(tmpSol, globalSol);
            this->setCurrentSolution(tmpSol);

            for (int faceId = 0; faceId < this->curElementGeom_.numEdges; faceId++)
            {
                int idxI = this->curElementGeom_.subContVolFace[faceId].i;

                int idxJ = this->curElementGeom_.subContVolFace[faceId].j;

                int flagI, flagJ;

                globalI = this->curElementGeom_.subContVol[idxI].global;
                globalJ = this->curElementGeom_.subContVol[idxJ].global;
                // 2D case: give y or x value of the line over which flux is to be
                //            calculated.
                // up to now only flux calculation to lines or planes (3D) parallel to
                // x, y and z axis possible

                // Flux across plane with z = 80 numEq
                if (globalI[coord] < coordVal)
                    flagI = 1;
                else
                    flagI = -1;

                if (globalJ[coord] < coordVal)
                    flagJ = 1;
                else
                    flagJ = -1;

                if (flagI == flagJ)
                {
                    sign = 0;
                }
                else
                {
                    if (flagI > 0)
                        sign = -1;
                    else
                        sign = 1;
                }

                // get variables

                if (flagI != flagJ)
                {
                    computeFlux(tmpFlux, faceId);
                    tmpFlux *= sign;
                    flux += tmpFlux;
                }
            }
        }

        flux = this->problem_.gridView().comm().sum(flux);
    }

    /*!
     * \brief Calculate mass of both components in the whole model domain
     *         and get minimum and maximum values of primary variables
     *
     */
    void calculateMass(const SolutionVector &globalSol, Dune::FieldVector<
            Scalar, 2> &mass)
    {
        mass = 0;

        ElementIterator elementIt =
                this->model().gridView().template begin<0> ();
        ElementIterator endit = this->model().gridView().template end<0> ();

        SolutionOnElement curSol;
        VertexDataArray elemDat;

        Scalar minSat = 1e100;
        Scalar maxSat = -1e100;
        Scalar minP = 1e100;
        Scalar maxP = -1e100;
#if !ISOTHERMAL
        Scalar minTe = 1e100;
        Scalar maxTe = -1e100;
#endif

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

            for (int i = 0; i < numLocalVerts; ++i)
            {
                const VertexData &vdat = elemDat[i];
                Scalar vol = this->curElementGeom_.subContVol[i].volume;

                Scalar satN = vdat.saturation(nPhaseIdx);
                Scalar pW = vdat.pressure(wPhaseIdx);
                Scalar T = vdat.temperature();

                // get minimum and maximum values of primary variables
                minSat = std::min(minSat, satN);
                maxSat = std::max(maxSat, satN);
                minP = std::min(minP, pW);
                maxP = std::max(maxP, pW);
#if !ISOTHERMAL
                minTe = std::min(minTe, T);
                maxTe = std::max(maxTe, T);
#endif

                mass[nPhaseIdx] += vdat.porosity() * vdat.saturation(nPhaseIdx)
                        * vdat.density(nPhaseIdx) * vol;

                mass[wPhaseIdx] += vdat.porosity() * vdat.saturation(wPhaseIdx)
                        * vdat.density(wPhaseIdx) * vol;
            }
        }

        mass = this->gridView_.comm().sum(mass);

        Scalar minS = this->gridView_.comm().min(minSat);
        Scalar maxS = this->gridView_.comm().max(maxSat);
        Scalar minPr = this->gridView_.comm().min(minP);
        Scalar maxPr = this->gridView_.comm().max(maxP);
        Scalar minT = this->gridView_.comm().min(minTe);
        Scalar maxT = this->gridView_.comm().max(maxTe);

        if (this->problem_.gridView().comm().rank() == 0) // IF PARALLEL: only print by processor with rank() == 0
        {
            // print minimum and maximum values
            std::cout << "nonwetting phase saturation: min = " << minS
                    << ", max = " << maxS << std::endl;
            std::cout << "wetting phase pressure: min = " << minPr
                    << ", max = " << maxPr << std::endl;
#if !ISOTHERMAL
            std::cout << "temperature: min = " << minT << ", max = " << maxT
                    << std::endl;
#endif
        }
    }

    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer,
            const SolutionVector &globalSol)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->problem_.gridView().size(dim);
        ScalarField *pW = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *pN = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *pC = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *Sw = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *Sn = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *Te = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *rhoW =
                writer.template createField<Scalar, 1> (numVertices);
        ScalarField *rhoN =
                writer.template createField<Scalar, 1> (numVertices);
        ScalarField *mobW =
                writer.template createField<Scalar, 1> (numVertices);
        ScalarField *mobN =
                writer.template createField<Scalar, 1> (numVertices);

        unsigned numElements = this->gridView_.size(0);
        ScalarField *rank =
                writer.template createField<Scalar, 1> (numElements);

        const DofEntityMapper &dofMapper =
                this->problem_.model().dofEntityMapper();
        SolutionOnElement tmpSol;
        ElementIterator elementIt =
                this->problem_.gridView().template begin<0> ();
        ElementIterator endit = this->problem_.gridView().template end<0> ();
        for (; elementIt != endit; ++elementIt)
        {
            int idx = this->problem_.model().elementMapper().map(*elementIt);
            (*rank)[idx] = this->gridView_.comm().rank();

            setCurrentElement(*elementIt);
            this->restrictToElement(tmpSol, globalSol);
            this->setCurrentSolution(tmpSol);

            int numVerts = elementIt->template count<dim> ();
            for (int i = 0; i < numVerts; ++i)
            {
                int globalIdx = dofMapper.map(*elementIt, i, dim);

                (*pW)[globalIdx] = this->curElemDat_[i].pressure(wPhaseIdx);
                (*pN)[globalIdx] = this->curElemDat_[i].pressure(nPhaseIdx);
                (*pC)[globalIdx] = this->curElemDat_[i].capillaryPressure();
                (*Sw)[globalIdx] = this->curElemDat_[i].saturation(wPhaseIdx);
                (*Sn)[globalIdx] = this->curElemDat_[i].saturation(nPhaseIdx);
                (*rhoW)[globalIdx] = this->curElemDat_[i].density(wPhaseIdx);
                (*rhoN)[globalIdx] = this->curElemDat_[i].density(nPhaseIdx);
                (*mobW)[globalIdx] = this->curElemDat_[i].mobility(wPhaseIdx);
                (*mobN)[globalIdx] = this->curElemDat_[i].mobility(nPhaseIdx);
                (*Te)[globalIdx] = this->curElemDat_[i].temperature();
            };
        }

        writer.addVertexData(Sn, "Sn");
        writer.addVertexData(Sw, "Sw");
        writer.addVertexData(pW, "pg");
        writer.addVertexData(pN, "pn");
        writer.addVertexData(pC, "pc");
        writer.addVertexData(rhoW, "rhoW");
        writer.addVertexData(rhoN, "rhoN");
        writer.addVertexData(mobW, "mobW");
        writer.addVertexData(mobN, "mobN");
        writer.addVertexData(Te, "temperature");
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

        asImp_()->addOutputVtkFields(writer, oldSol);
    }

    /*!
     * \brief Returns the epsilon value which is added and removed
     *        from the current solution.
     *
     * \param elemSol    The current solution on the element
     * \param vertexIdx  The local index of the element's vertex for
     *                   which the local derivative ought to be calculated.
     * \param pvIdx      The index of the primary variable which gets varied
     */
    //    Scalar numericEpsilon_(const SolutionOnElement &elemSol,
    //                           int vertIdx,
    //                           int pvIdx) const
    //    {
    //        if (pvIdx == pressureIdx)
    //            return 1e1;
    //        else
    //            return 1e-6;
    //    }

protected:
    Implementation *asImp_()
    {
        return static_cast<Implementation *> (this);
    }
    const Implementation *asImp_() const
    {
        return static_cast<const Implementation *> (this);
    }
};

}

#endif
