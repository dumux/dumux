// $Id: 2pboxmodel.hh 3738 2010-06-15 14:01:09Z lauser $
/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2007-2007 by Bernd Flemisch                               *
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
#ifndef DUMUX_TWOP_MODEL_HH
#define DUMUX_TWOP_MODEL_HH

#include "2plocalresidual.hh"
#include "2pnewtoncontroller.hh"
#include "2pproblem.hh"

namespace Dumux
{

/*!
 * \ingroup BoxProblems
 * \defgroup TwoPBoxProblems Two-phase box problems
 */

/*!
 * \ingroup BoxModels
 * \defgroup TwoPBoxModel Two-phase box model
 */

/*!
 * \ingroup TwoPBoxModel
 * \brief Adaption of the BOX scheme to the twophase flow model.
 *
 * This model implements two-phase flow of two completely immiscible fluids
 * \f$\alpha \in \{ w, n \}\f$ using a standard multiphase Darcy
 * approach as the equation for the conservation of momentum:
 \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} K
 \left(\text{grad} p_\alpha - \varrho_{\alpha} \boldsymbol{g} \right)
 \f]
 *
 * By inserting this into the equation for the conservation of the
 * phase mass, one gets
 \f[
 \phi \frac{\partial \varrho_\alpha S_\alpha}{\partial t}
 -
 \text{div} \left\{
 \varrho_\alpha \frac{k_{r\alpha}}{\mu_\alpha} K \left(\text{grad} p_\alpha - \varrho_{\alpha} \boldsymbol{g} \right)
 \right\} = q_\alpha \;,
 \f]
 * discretized by a fully-coupled vertex
 * centered finite volume (box) scheme as spatial and
 * the implicit Euler method as time discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$, the number of
 * unknowns can be reduced to two. Currently the model supports
 * choosing either \f$p_w\f$ and \f$S_n\f$ or \f$p_n\f$ and \f$S_w\f$
 * as primary variables. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * <tt>TwoPCommonIndices::pWsN</tt> or <tt>TwoPCommonIndices::pNsW</tt>. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 */
template<class TypeTag >
class TwoPModel : public BoxModel<TypeTag>
{
    typedef TwoPModel<TypeTag> ThisType;
    typedef BoxModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SecondaryVars)) SecondaryVars;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementSecondaryVars)) ElementSecondaryVars;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementBoundaryTypes)) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexMapper)) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementMapper)) ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    enum {
        dim = GridView::dimension,
        nPhaseIdx = Indices::nPhaseIdx,
        wPhaseIdx = Indices::wPhaseIdx,
    };

public:   
    /*!
     * \brief Calculate mass of both components in the whole model domain
     *         and get minimum and maximum values of primary variables
     *
     */
    void calculateMass(const SolutionVector &curSol, Dune::FieldVector<Scalar, 2> &mass)
    {
        mass = 0;

        ElementIterator elementIt =
                this->model().gridView().template begin<0> ();
        ElementIterator endit = this->model().gridView().template end<0> ();

        ElementSecondaryVars elemDat;

        Scalar minSat = 1e100;
        Scalar maxSat = -1e100;
        Scalar minP = 1e100;
        Scalar maxP = -1e100;
#if !ISOTHERMAL
        Scalar minTe = 1e100;
        Scalar maxTe = -1e100;
#endif

        FVElementGeometry fvElemGeom;
        SecondaryVars secVars;
        ElementBoundaryTypes elemBcTypes;

        ElementIterator elemIt = this->gridView_().template begin<0>();
        ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            int idx = this->problem_.model().elementMapper().map(*elemIt);
            fvElemGeom.update(this->gridView_(), *elemIt);
            elemBcTypes.update(this->problem_(), *elemIt, fvElemGeom);

            int numLocalVerts = elementIt->template count<dim> ();
            for (int i = 0; i < numLocalVerts; ++i)
            {
                int globalIdx = this->vertexMapper().map(*elemIt, i, dim);
                secVars.update(sol[globalIdx], 
                               this->problem_(),
                               *elemIt,
                               fvElemGeom, 
                               i,
                               false);

                Scalar vol = fvElemGeom.subContVol[i].volume;

                Scalar satN = secVars.saturation(nPhaseIdx);
                Scalar pW = secVars.pressure(wPhaseIdx);
                Scalar T = secVars.temperature();

                // get minimum and maximum values of primary variables
                minSat = std::min(minSat, satN);
                maxSat = std::max(maxSat, satN);
                minP = std::min(minP, pW);
                maxP = std::max(maxP, pW);
#if !ISOTHERMAL
                minTe = std::min(minTe, T);
                maxTe = std::max(maxTe, T);
#endif

                mass[nPhaseIdx] += secVars.porosity() * secVars.saturation(nPhaseIdx)
                        * secVars.density(nPhaseIdx) * vol;

                mass[wPhaseIdx] += secVars.porosity() * secVars.saturation(wPhaseIdx)
                        * secVars.density(wPhaseIdx) * vol;
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
    void addOutputVtkFields(const SolutionVector &sol, 
                            MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->problem_().gridView().size(dim);
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

        unsigned numElements = this->gridView_().size(0);
        ScalarField *rank =
                writer.template createField<Scalar, 1> (numElements);

        FVElementGeometry fvElemGeom;
        SecondaryVars secVars;
        ElementBoundaryTypes elemBcTypes;

        ElementIterator elemIt = this->gridView_().template begin<0>();
        ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            int idx = this->problem_().model().elementMapper().map(*elemIt);
            (*rank)[idx] = this->gridView_().comm().rank();

            fvElemGeom.update(this->gridView_(), *elemIt);
            elemBcTypes.update(this->problem_(), *elemIt, fvElemGeom);

            int numVerts = elemIt->template count<dim> ();
            for (int i = 0; i < numVerts; ++i)
            {
                int globalIdx = this->vertexMapper().map(*elemIt, i, dim);
                secVars.update(sol[globalIdx], 
                               this->problem_(),
                               *elemIt,
                               fvElemGeom, 
                               i,
                               false);
                
                (*pW)[globalIdx] = secVars.pressure(wPhaseIdx);
                (*pN)[globalIdx] = secVars.pressure(nPhaseIdx);
                (*pC)[globalIdx] = secVars.capillaryPressure();
                (*Sw)[globalIdx] = secVars.saturation(wPhaseIdx);
                (*Sn)[globalIdx] = secVars.saturation(nPhaseIdx);
                (*rhoW)[globalIdx] = secVars.density(wPhaseIdx);
                (*rhoN)[globalIdx] = secVars.density(nPhaseIdx);
                (*mobW)[globalIdx] = secVars.mobility(wPhaseIdx);
                (*mobN)[globalIdx] = secVars.mobility(nPhaseIdx);
                (*Te)[globalIdx] = secVars.temperature();
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
     * \brief Calculate the flux of the nonwetting phase across a given
     * layer for the current timestep
     */
    void calculateFluxAcrossLayer(Dune::FieldVector<Scalar, 2> &flux, int coord, Scalar coordVal)
    {
        this->localJacobian().calculateFluxAcrossLayer(this->curSolFunction(), flux, coord, coordVal);
    }
    
    /*!
     * \brief Calculate the phase masses in the system for the current
     *        timestep.
     */
    void calculateMass(Dune::FieldVector<Scalar, 2> &mass)
    {
        this->localJacobian().calculateMass(this->curSolFunction(), mass);
    }
};
}

#endif
