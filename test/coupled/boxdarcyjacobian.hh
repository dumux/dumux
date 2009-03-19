#ifndef DUNE_BOXDARCYJACOBIAN_HH
#define DUNE_BOXDARCYJACOBIAN_HH

#include<map>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<sstream>

#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/quadraturerules.hh>
#include <dune/grid/utility/intersectiongetter.hh>

#include<dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dumux/operators/boxjacobian.hh>

namespace Dune
{
//! A class for computing local jacobian matrices
/*! A class for computing local jacobian matrix for the
  diffusion equation

  div j = q; j = -K grad u; in Omega

  u = g on Gamma1; j*n = J on Gamma2.

  Uses the box method.

  Template parameters are:

  - Grid     a DUNE grid type
  - Scalar    type used for return values
*/
template<class Grid, class Scalar, class BoxFunction = LeafP1Function<Grid, Scalar, 1> >
class BoxDarcyJacobian
    : public BoxJacobian<BoxDarcyJacobian<Grid,Scalar,BoxFunction>,Grid,Scalar,1,BoxFunction>
{
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Element::Geometry Geometry;
    typedef BoxDarcyJacobian<Grid,Scalar,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,Grid,Scalar,1>::VBlockType SolutionVector;
    typedef BoxJacobian<ThisType,Grid,Scalar,1,BoxFunction> BoxJacobianType;
    typedef Dune::FVElementGeometry<Grid> FVElementGeometry;

public:
    enum {dim=Grid::dimension};

    //! Constructor
    BoxDarcyJacobian (CoupledPorousMediaProblem<Grid,Scalar>& params,
                          bool levelBoundaryAsDirichlet_, const Grid& grid,
                          BoxFunction& sol,
                          bool procBoundaryAsDirichlet_=true)
        : BoxJacobian<ThisType,Grid,Scalar,1,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
          problem(params)
    {
        this->analytic = false;
    }

    void clearVisited ()
    {
        return;
    }

    template<class TypeTag>
    void localDefect(const Element& element, const SolutionVector* sol, bool withBC = true) {
        BoxJacobianType::template localDefect<TypeTag>(element, sol, withBC);

        this->template assembleBC<TypeTag>(element);

        Dune::GeometryType gt = element.geometry().type();
        const typename ReferenceElementContainer<Scalar,dim>::value_type& referenceElement = ReferenceElements<Scalar, dim>::general(gt);

        for (int vert=0; vert < this->fvGeom.numVertices; vert++) // begin loop over vertices / sub control volumes
            if (!this->fvGeom.subContVol[vert].inner)
            {
                typedef typename IntersectionIteratorGetter<Grid,TypeTag>::IntersectionIterator IntersectionIterator;

                FieldVector<Scalar,dim> averagedNormal(0);
                int faces = 0;
                IntersectionIterator endit = IntersectionIteratorGetter<Grid, TypeTag>::end(element);
                for (IntersectionIterator it = IntersectionIteratorGetter<Grid, TypeTag>::begin(element); it!=endit; ++it)
                {
                    if (it->boundary()) {
                        int faceIdx = it->numberInSelf();
                        int numVerticesOfFace = referenceElement.size(faceIdx, 1, dim);
                        for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++) {
                            int nodeInElement = referenceElement.subEntity(faceIdx, 1, nodeInFace, dim);
                            if (nodeInElement != vert)
                                continue;

                            faces++;
                        }
                    }
                }

                if (faces == 2 && this->fvGeom.numVertices == 4)
                {
                    this->def[vert] = sol[0] + sol[3] - sol[1] - sol[2];
                }
            }

        return;
    }

    SolutionVector computeM (const Element& e, const SolutionVector* sol, int node, bool old = false)
    {
        SolutionVector result(0);
        return result;
    }

    SolutionVector computeA (const Element& e, const SolutionVector* sol, int face)
    {
        FieldVector<Scalar, dim> gradP(0);
        for (int k = 0; k < this->fvGeom.numVertices; k++) {
            FieldVector<Scalar,dim> grad(this->fvGeom.subContVolFace[face].grad[k]);
            grad *= sol[k];
            gradP += grad;
        }

        FieldVector<Scalar,dim> KGradP(0);
        elData.K.umv(gradP, KGradP);
        SolutionVector flux = KGradP*this->fvGeom.subContVolFace[face].normal;

        return flux;
    }

    SolutionVector computeQ (const Element& e, const SolutionVector* sol, const int& node)
    {
        return problem.q(this->fvGeom.subContVol[node].global, e, this->fvGeom.subContVol[node].local);
    }

    void computeElementData (const Element& e)
    {
        // ASSUMING element-wise constant permeability, evaluate K at the cell center
        elData.K = problem.K(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);
    };

    virtual void updateVariableData(const Element& e, const SolutionVector* sol, int i, bool old = false)
    {
        return;
    }

    void updateVariableData(const Element& e, const SolutionVector* sol, bool old = false)
    {
        return;
    }

    virtual void updateStaticData (const Element& e, const SolutionVector* sol)
    {
        return;
    }

    struct ElementData {
        FieldMatrix<Scalar,dim,dim> K;
    };

    ElementData elData;
    CoupledPorousMediaProblem<Grid,Scalar>& problem;
};
}
#endif
