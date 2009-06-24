// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
#ifndef DUMUX_ONEP_TWOC_BOX_JACOBIAN_HH
#define DUMUX_ONEP_TWOC_BOX_JACOBIAN_HH

#include <dumux/boxmodels/boxscheme/boxscheme.hh>

#include <dumux/boxmodels/1p2c/1p2cproperties.hh>

#include <dumux/boxmodels/1p2c/1p2cvertexdata.hh>
#include <dumux/boxmodels/1p2c/1p2celementdata.hh>
#include <dumux/boxmodels/1p2c/1p2cfluxdata.hh>

#include <dune/common/collectivecommunication.hh>
#include <vector>
#include <iostream>

namespace Dune
{
/*!
 * \brief Calculate the local Jacobian for the single-phase,
 *        two-component model in the BOX scheme.
 */
template<class TypeTag>
class OnePTwoCBoxJacobian : public BoxJacobian<TypeTag, OnePTwoCBoxJacobian<TypeTag> >
{
protected:
    typedef OnePTwoCBoxJacobian<TypeTag>     ThisType;
    typedef BoxJacobian<TypeTag, ThisType>   ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))      Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePTwoCIndices))  Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))       Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))     GridView;

    enum {
        dim            = GridView::dimension,
        dimWorld       = GridView::dimensionworld,

        numEq          = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        numPhases      = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
        numComponents  = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)),

        konti          = Indices::konti,
        transport      = Indices::transport,
    };

    typedef typename GridView::template Codim<0>::Entity   Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef FieldVector<Scalar, dim>       LocalPosition;
    typedef FieldVector<Scalar, dimWorld>  GlobalPosition;

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVarVector        PrimaryVarVector;
    typedef typename SolutionTypes::SolutionFunction        SolutionFunction;
    typedef typename SolutionTypes::SolutionOnElement       SolutionOnElement;

    typedef Dune::FieldVector<Scalar, numPhases> PhasesVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexData))   VertexData;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementData))  ElementData;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxData))     FluxData;

    typedef std::vector<VertexData>        VertexDataArray;
    typedef FieldMatrix<Scalar, dim, dim>  Tensor;

    static const Scalar upwindAlpha = GET_PROP_VALUE(TypeTag, PTAG(UpwindAlpha));

public:
    OnePTwoCBoxJacobian(Problem &problem)
        : ParentType(problem)
    {
    };

    /*!
     * \brief Evaluate the amount all conservation quantites
     *        (e.g. phase mass) within a finite volume.
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

        // storage term of continuity equation
        result[konti] = 0;

        // storage term of the transport equation
        result[transport] = vertDat.density * vertDat.porosity * vertDat.molefraction;
    }

    /*!
     * \brief Evaluates the mass flux over a face of a subcontrol
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

        // data attached to upstream and the downstream vertices
        // of the current phase
        const VertexData &up = this->curElemDat_[vars.upstreamIdx];
        const VertexData &dn = this->curElemDat_[vars.downstreamIdx];

        flux[konti] = vars.vDarcyNormal / vars.viscosityAtIP;

        // advective flux
        flux[transport] +=
            vars.vDarcyNormal *
            (  upwindAlpha*
               (  up.density * up.molefraction/up.viscosity )
               +
               (1 - upwindAlpha)*
               (  dn.density * dn.molefraction/dn.viscosity ) );

        // diffusive flux
        flux[transport] +=
            vars.densityAtIP * vars.diffCoeffPM *
            (vars.concentrationGrad*vars.face->normal);
    }

    /*!
     * \brief Calculate the source term of the equation
     */
    void computeSource(PrimaryVarVector &q, int localVertexIdx)
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
    template <class PrimaryVarVector>
    Scalar temperature(const PrimaryVarVector &sol)
    {
        return this->problem_.temperature(); /* constant temperature */
    }

    /*!
     * \brief Add the mass fraction of air in water to VTK output of
     *        the current timestep.
     */
    template <class MultiWriter>
    void addVtkFields(MultiWriter &writer, const SolutionFunction &globalSol)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->gridView_.size(dim);

        ScalarField *pressure = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *molefraction = writer.template createField<Scalar, 1>(numVertices);

        SolutionOnElement tmpSol;
        VertexDataArray   elemDat;

        ElementIterator elementIt = this->gridView_.template begin<0>();
        const ElementIterator &endit = this->gridView_.template end<0>();
        for (; elementIt != endit; ++elementIt)
        {
            int numLocalVerts = elementIt->template count<dim>();
            tmpSol.resize(numLocalVerts);

            setCurrentElement(*elementIt);
            this->restrictToElement(tmpSol, globalSol);
            updateElementData_(elemDat, tmpSol, false);

            for (int i = 0; i < numLocalVerts; ++i)
            {
                int globalIdx = this->problem_.model().vertexMapper().map(*elementIt, i, dim);

                (*pressure)[globalIdx] = elemDat[i].pressure;
                (*molefraction)[globalIdx] = elemDat[i].molefraction;

            };
        }

        writer.addVertexData(pressure, "pressure");
        writer.addVertexData(molefraction, "molefraction");

    }
};

}

#endif
