// $Id$
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
#ifndef DUMUX_PDELAB_BOX_LOCAL_OPERATOR_HH
#define DUMUX_PDELAB_BOX_LOCAL_OPERATOR_HH

#include<vector>
#include<dune/common/fvector.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/quadraturerules.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspace.hh>
#include<dune/pdelab/gridoperatorspace/gridoperatorspaceutilities.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>

namespace Dumux {

namespace PDELab {

template<class TypeTag>
class BoxLocalOperator
    :
//    : public Dune::PDELab::NumericalJacobianApplyVolume<BoxLocalOperatorPDELab<TypeTag> >,
//public Dune::PDELab::NumericalJacobianVolume<BoxLocalOperatorPDELab<TypeTag> >,
    public Dune::PDELab::FullVolumePattern,
    public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
    // pattern assembly flags
    enum { doPatternVolume = true };

    // residual assembly flags
    enum { doAlphaVolume = true };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Model;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementSolutionVector)) ElementSolutionVector;
    enum{numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq))};

    BoxLocalOperator(Model &model)
        : model_(model)
    {}

    // volume integral depending on test and ansatz functions
    template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x,
                       const LFSV& lfsv, R& r) const
    {
        typedef typename LFSU::Traits::SizeType size_type;

        model_.localResidual().eval(eg.entity());
        
        int numVertices = x.size()/numEq;
        for (size_type comp = 0; comp < r.size(); comp++)
            r[comp] = model_.localResidual().residual(comp%numVertices)[comp/numVertices];
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

} // namespace PDELab
} // namespace Dumux

#endif
