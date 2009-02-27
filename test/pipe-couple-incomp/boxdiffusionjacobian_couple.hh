#ifndef DUNE_BOXDIFFUSIONJACOBIAN_HH
#define DUNE_BOXDIFFUSIONJACOBIAN_HH

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
#include "diffusionparameters_couple.hh"

namespace Dune
{
//! A class for computing local jacobian matrices
/*! A class for computing local jacobian matrix for the
    diffusion equation

        div j = q; j = -K grad u; in Omega

        u = g on Gamma1; j*n = J on Gamma2.

    Uses the box method.

    Template parameters are:

    - G     a DUNE grid type
    - RT    type used for return values
 */
template<class G, class RT, class GlobalToPipeMapper, class VertexMapper, class VertexVectorOnLineType, class BoxFunction = LeafP1Function<G, RT, 1> >
class BoxDiffusionJacobian
: public BoxJacobian<BoxDiffusionJacobian<G,RT,GlobalToPipeMapper,VertexMapper,VertexVectorOnLineType,BoxFunction>,G,RT,1,BoxFunction>
{
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef BoxDiffusionJacobian<G,RT, GlobalToPipeMapper,VertexMapper,VertexVectorOnLineType,BoxFunction> ThisType;
    typedef typename LocalJacobian<ThisType,G,RT,1>::VBlockType VBlockType;
    typedef Dune::FVElementGeometry<G> FVElementGeometry;

public:
    enum {n=G::dimension};

    //! Constructor
    BoxDiffusionJacobian (DiffusionParameters<G,RT, GlobalToPipeMapper, VertexMapper, VertexVectorOnLineType>& params,
            bool levelBoundaryAsDirichlet_, const G& grid,
            BoxFunction& sol,
            bool procBoundaryAsDirichlet_=true)
            : BoxJacobian<ThisType,G,RT,1,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
            problem(params)
            {
        this->analytic = false;
            }

    void clearVisited ()
    {
        return;
    }

    VBlockType computeM (const Entity& e, const VBlockType* sol, int node, bool old = false)
    {
          int globalId = problem.vertexMapper.template map<n>(e, node);
          for (unsigned k = 0; k < problem.vertexVectorOnLine.size(); k++)
          {
              if (globalId == problem.vertexVectorOnLine[k].globalId)
              {
                  VBlockType result(problem.alphaEx*sol[node]);
                  return result;
              }
          }

          VBlockType result(0);
          return result;
    }

    VBlockType computeQ (const Entity& e, const VBlockType* sol, const int& node)
    {
        return problem.q(this->fvGeom.subContVol[node].global, e, this->fvGeom.subContVol[node].local, node);
    }

    VBlockType computeA (const Entity& e, const VBlockType* sol, int face)
    {
        FieldVector<RT, n> gradP(0);
        for (int k = 0; k < this->fvGeom.numVertices; k++) {
            FieldVector<RT,n> grad(this->fvGeom.subContVolFace[face].grad[k]);
            grad *= sol[k];
            gradP += grad;
        }

        FieldVector<RT,n> gravity = problem.gravity();
        double densityFace =  problem.materialLaw().density(problem.Temp());

        gravity *= (-densityFace);
        gradP += gravity;

        FieldVector<RT,n> KGradP(0);
        elData.K.umv(gradP, KGradP);

        VBlockType flux = KGradP*this->fvGeom.subContVolFace[face].normal;

        return flux;
    }

    void computeElementData (const Entity& e)
    {
        // ASSUMING element-wise constant permeability, evaluate K at the cell center
        elData.K = problem.K(this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal);
    };

     virtual void updateVariableData(const Entity& e, const VBlockType* sol, int i, bool old = false)
    {
         return;
    }

    void updateVariableData(const Entity& e, const VBlockType* sol, bool old = false)
    {
        return;
    }

    void updateStaticData (const Entity& e, const VBlockType* sol)
    {
        return;
    }

    struct ElementData {
        FieldMatrix<DT,n,n> K;
    };

    ElementData elData;
    DiffusionParameters<G,RT, GlobalToPipeMapper, VertexMapper, VertexVectorOnLineType>& problem;
};
}
#endif
