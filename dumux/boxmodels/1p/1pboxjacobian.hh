// $Id: 1pboxjacobian.hh 3738 2010-06-15 14:01:09Z lauser $
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
#ifndef DUMUX_1P_BOX_JACOBIAN_HH
#define DUMUX_1P_BOX_JACOBIAN_HH

#include <dumux/boxmodels/boxscheme/boxjacobian.hh>

#include "1pvertexdata.hh"
#include "1pelementdata.hh"
#include "1pfluxdata.hh"

namespace Dumux
{
/*!
 * \ingroup OnePBoxModel
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the one-phase box model.
 */
template<class TypeTag>
class OnePBoxJacobian : public BoxJacobian<TypeTag>
{
    typedef OnePBoxJacobian<TypeTag>         ThisType;
    typedef BoxJacobian<TypeTag>   ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))   Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))  GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))    Scalar;

    typedef typename GridView::template Codim<0>::Entity   Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::SolutionVector        SolutionVector;
    typedef typename SolutionTypes::SolutionOnElement       SolutionOnElement;
    typedef typename SolutionTypes::PrimaryVarVector        PrimaryVarVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePIndices)) Indices;

    enum {
        dim              = GridView::dimension,
        dimWorld         = GridView::dimensionworld,

        pressureIdx      = Indices::pressureIdx,
    };


    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexData))   VertexData;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluxData))     FluxData;
    typedef std::vector<VertexData> VertexDataArray;

    typedef Dune::FieldVector<Scalar, dim>       LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;

public:
    OnePBoxJacobian(Problem &problem)
        : ParentType(problem)
    {};

    /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub control
     *        volume of a finite volume element for the OneP
     *        model.
     *
     * This function should not include the source and sink terms.
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

        // partial time derivative of the wetting phase mass
        result[pressureIdx] =  vertDat.density * vertDat.porosity;
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

        flux[pressureIdx] = vars.densityAtIP * vars.vDarcyNormal / vars.viscosityAtIP;
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
    { return this->problem_.temperature(); /* constant temperature */ }

    /*!
     * \brief All relevant primary and secondary of a given
     *        solution to an ouput writer.
     */
    template <class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer, const SolutionVector &globalSol)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->gridView_.size(dim);
        ScalarField *p = writer.template createField<Scalar, 1>(numVertices);

        unsigned numElements = this->gridView_.size(0);
        ScalarField *rank = writer.template createField<Scalar, 1>(numElements);

        SolutionOnElement tmpSol;
        ElementIterator elementIt = this->gridView_.template begin<0>();
        const ElementIterator &endit = this->gridView_.template end<0>();;
        for (; elementIt != endit; ++elementIt)
        {
            int idx = this->problem_.model().elementMapper().map(*elementIt);
            (*rank)[idx] = this->gridView_.comm().rank();

            setCurrentElement(*elementIt);
            this->restrictToElement(tmpSol, globalSol);
            this->setCurrentSolution(tmpSol);

            int numVerts = elementIt->template count<dim>();
            for (int i = 0; i < numVerts; ++i)
            {
                int globalIdx = this->problem_.model().vertexMapper().map(*elementIt, i, dim);
                (*p)[globalIdx] = this->curElemDat_[i].pressure;
            };
        }

        writer.addVertexData(p, "p");
        writer.addCellData(rank, "process rank");
    }

private:
    ThisType &asImp_()
    { return *static_cast<ThisType *>(this); }

    const ThisType &asImp_() const
    { return *static_cast<const ThisType *>(this); }
};

};

#endif
