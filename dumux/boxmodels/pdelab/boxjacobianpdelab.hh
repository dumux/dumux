// $Id: boxjacobianpdelab.hh 3736 2010-06-15 09:52:10Z lauser $
/*****************************************************************************
 *   Copyright (C) 2009-2010 by Bernd Flemisch                               *
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
#ifdef HAVE_DUNE_PDELAB

#ifndef DUMUX_BOXJACOBIANPDELAB_HH
#define DUMUX_BOXJACOBIANPDELAB_HH

#include<vector>
#include<dune/common/fvector.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/quadraturerules.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>

template<class TypeTag>
class BoxJacobianPDELab
    :
//    : public Dune::PDELab::NumericalJacobianApplyVolume<BoxJacobianPDELab<TypeTag> >,
//public Dune::PDELab::NumericalJacobianVolume<BoxJacobianPDELab<TypeTag> >,
    public Dune::PDELab::FullVolumePattern,
    public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
    // pattern assembly flags
    enum { doPatternVolume = true };

    // residual assembly flags
    enum { doAlphaVolume = true };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model))    Model;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))  Problem;
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::SolutionOnElement       SolutionOnElement;
    enum{numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq))};

    BoxJacobianPDELab (Model &model)
        : model_(model)
    {}

    // volume integral depending on test and ansatz functions
    template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x,
            const LFSV& lfsv, R& r) const
    {
        typedef typename LFSU::Traits::SizeType size_type;

        model_.localJacobian().setCurrentElement(eg.entity());

        int numVertices = x.size()/numEq;

        SolutionOnElement localU(numVertices);
        model_.localJacobian().restrictToElement(localU, model_.curSol());
        model_.localJacobian().setCurrentSolution(localU);

        SolutionOnElement localUOld(numVertices);
        model_.localJacobian().restrictToElement(localUOld, model_.prevSol());
        model_.localJacobian().setPreviousSolution(localUOld);

        SolutionOnElement localResidual(numVertices);
        model_.localJacobian().evalLocalResidual(localResidual);
        for (size_type comp = 0; comp < r.size(); comp++)
            r[comp] = localResidual[comp%numVertices][comp/numVertices];
    }

    // jacobian of volume term
    template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
    void jacobian_volume (const EG& eg,
                          const LFSU& lfsu,
                          const X& x,
                          const LFSV& lfsv,
                          Dune::PDELab::LocalMatrix<R>& mat) const
    {
        typedef typename LFSU::Traits::SizeType size_type;

        model_.localJacobian().assemble(eg.entity());

        int numVertices = x.size()/numEq;
        for (size_type j=0; j<lfsu.size(); j++)
          for (size_type i=0; i<lfsu.size(); i++)
              mat(i,j) = (model_.localJacobian().mat(i%numVertices,j%numVertices))[i/numVertices][j/numVertices];
    }

private:
    Model& model_;
};

#endif

#endif
