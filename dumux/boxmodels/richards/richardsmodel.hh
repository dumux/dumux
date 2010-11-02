// $Id$
/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2009 by Onur Dogan                                        *
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

/*!
* \file
*
* \brief Adaption of the box scheme to the Richards model.
*/

#ifndef DUMUX_RICHARDS_MODEL_HH
#define DUMUX_RICHARDS_MODEL_HH

#include <dumux/boxmodels/common/boxmodel.hh>

#include "richardslocalresidual.hh"
#include "richardsproblem.hh"

namespace Dumux
{
/*!
 * \ingroup BoxModels
 * \defgroup RichardsModel Richards box model
 */

/*!
 * \ingroup RichardsModel
 * \brief Implements the Richards model for quasi twophase flow.
 *
 * In the unsaturated zone, Richards' equation can be used. Fundamentally, the Richards-equation
 * is equivalent to the twophase model, i.e.
 \f[
 \frac{\partial\;\phi S_\alpha \rho_\alpha}{\partial t}
 -
 \mathbf{div} \left\{
 \frac{k_{r\alpha}}{\mu_\alpha}\;K
 \mathbf{grad}\left[
 p_alpha - g\rho_\alpha
 \right]
 \right\}
 =
 q_\alpha,
 \f]
 * where \f$\alpha \in \{w, n\}\f$ is the fluid phase, \f$\rho_\alpha\f$ is the fluid 
 * density, \f$S_\alpha\f$ is the fluid saturation, \f$\phi\f$ is the porosity, 
 * \f$k_{r\alpha}\f$ is the relative permeability of the fluid, \f$\mu_\alpha\f$ is
 * the fluid's dynamic viscosity, \f$K\f$ is the intrinsic permeability, \f$p_\alpha\f$ 
 * is the fluid pressure and \f$g\f$ is the potential of the gravity.
 * 
 * However, the Richards model assumes that the non-wetting is gas
 * which typically exhibits a much low viscosity than the liquid
 * wetting phase. (For example air has about \f$1\%\f$ of the viscosity of
 * liquid water.) As a consequence, the
 * \f$\frac{k_{r\alpha}}{\mu_\alpha}\f$ term is typically much larger for
 * the non-wetting phase than for the wetting phase. In the Richards
 * model it is now assumed that \f$\frac{k_{rn}}{\mu_n} \to \infty\f$
 * which means that the pressure of the non-wetting phase is
 * equivalent to hydrostatic pressure of the gas or can be externally
 * specified.  Therefore, in Richards' equation mass conservation only
 * needs to be considered for the wetting phase. The model thus choses
 * \f$p_w\f$ as its only primary variable and calculates the wetting phase
 * saturation using the inverse of the capilary pressure, i.e.
 \f[ S_w = p_c^{-1}(p_n - p_w)\, \f]
 * where \f$p_n\f$ is an externally given reference pressure.
 */
template<class TypeTag >
class RichardsModel : public BoxModel<TypeTag>
{
    typedef RichardsModel<TypeTag> ThisType;
    typedef BoxModel<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementBoundaryTypes)) ElementBoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexMapper)) VertexMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementMapper)) ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(RichardsIndices)) Indices;
    enum {
        dim = GridView::dimension,
        nPhaseIdx = Indices::nPhaseIdx,
        wPhaseIdx = Indices::wPhaseIdx,
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

public:
    /*!
     * \brief Returns the relative weight of a primary variable for
     *        calculating relative errors.
     *
     * \param vertIdx The global index of the vertex in question
     * \param pvIdx The index of the primary variable
     */
    Scalar primaryVarWeight(int vertIdx, int pvIdx) const
    {
        if (Indices::pwIdx == pvIdx)
            return 1e-6;
        return 1;
    }

    /*!
     * \brief All relevant primary and secondary of a given
     *        solution to an ouput writer.
     *
     * \param sol The current solution which ought to be written to disk
     * \param writer The Dumux::VtkMultiWriter which is be used to write the data
     */
    template <class MultiWriter>
    void addOutputVtkFields(const SolutionVector &sol, MultiWriter &writer)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->problem_().gridView().size(dim);
        ScalarField *pW = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *pN = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *pC = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *Sw = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *Sn = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *rhoW = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *rhoN = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *mobW = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *mobN = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *poro = writer.template createField<Scalar, 1> (numVertices);
        ScalarField *Te = writer.template createField<Scalar, 1> (numVertices);

        unsigned numElements = this->gridView_().size(0);
        ScalarField *rank =
                writer.template createField<Scalar, 1> (numElements);

        FVElementGeometry fvElemGeom;
        VolumeVariables volVars;

        ElementIterator elemIt = this->gridView_().template begin<0>();
        ElementIterator elemEndIt = this->gridView_().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt)
        {
            int idx = this->problem_().model().elementMapper().map(*elemIt);
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

                (*pW)[globalIdx] = volVars.pressure(wPhaseIdx);
                (*pN)[globalIdx] = volVars.pressure(nPhaseIdx);
                (*pC)[globalIdx] = volVars.capillaryPressure();
                (*Sw)[globalIdx] = volVars.saturation(wPhaseIdx);
                (*Sn)[globalIdx] = volVars.saturation(nPhaseIdx);
                (*rhoW)[globalIdx] = volVars.density(wPhaseIdx);
                (*rhoN)[globalIdx] = volVars.density(nPhaseIdx);
                (*mobW)[globalIdx] = volVars.mobility(wPhaseIdx);
                (*mobN)[globalIdx] = volVars.mobility(nPhaseIdx);
                (*poro)[globalIdx] = volVars.porosity();
                (*Te)[globalIdx] = volVars.temperature();
            };
        }

        writer.addVertexData(Sn, "Sn");
        writer.addVertexData(Sw, "Sw");
        writer.addVertexData(pN, "pn");
        writer.addVertexData(pW, "pw");
        writer.addVertexData(pC, "pc");
        writer.addVertexData(rhoW, "rhoW");
        writer.addVertexData(rhoN, "rhoN");
        writer.addVertexData(mobW, "mobW");
        writer.addVertexData(mobN, "mobN");
        writer.addVertexData(poro, "porosity");
        writer.addVertexData(Te, "temperature");
        writer.addCellData(rank, "process rank");
    }
};
}

#include "richardspropertydefaults.hh"

#endif
