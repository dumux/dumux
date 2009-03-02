// $Id$

#ifndef DUNE_STOKESJACOBIAN_HH
#define DUNE_STOKESJACOBIAN_HH

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

#include"dumux/functions/mixedfunction.hh"
#include"dumux/operators/localjacobian.hh"
#include"dumux/stokes/stokesproblem.hh"

namespace Dune {

/** \todo Please doc me! */

template<class G, class RT, class MixedFunction = LeafMixedFunction<G, RT, 1> >
class StokesJacobian :
        public LocalJacobian<StokesJacobian<G,RT,MixedFunction>,G,RT,1> {
    // mapper: one data element per vertex
    template<int dim> struct P1Layout {
        bool contains(Dune::GeometryType gt) {
            return gt.dim() == 0;
        }
    };

    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename Entity::Geometry Geometry;
    typedef LocalJacobian<StokesJacobian<G,RT,MixedFunction>,G,RT,1> ThisType;
    typedef typename ThisType::VBlockType VBlockType;
    typedef typename ThisType::MBlockType MBlockType;
    typedef MultipleCodimMultipleGeomTypeMapper<G, typename G::Traits::LeafIndexSet, P1Layout> VertexMapper;

public:
    enum {dim=G::dimension};
    enum {SIZE=LagrangeShapeFunctionSetContainer<DT,RT,dim>::maxsize};

    //! Constructor
    StokesJacobian(StokesProblem<G,RT>& prob, bool dummy1, const G& grid, MixedFunction& sol, bool dummy2 = false) :
        vertexMapper(grid, grid.leafIndexSet()), problem(prob),
        currentSolution(sol), oldSolution(grid), dt(1) {
    }

    template<class TypeTag>
    void localDefect(const Entity& e, const VBlockType* sol, bool withBC = true) {
        for (int i=0; i < e.geometry().corners() + 1; i++)
            this->def[i] = 0.0;

        return;
    }

    void setLocalSolution(const Entity& e)
    {
        return;
    }

    void localToGlobal(const Entity& e, const VBlockType* sol)
    {
        return;
    }


    void setDt(double d) {
        dt = d;

        return;
    }

    double getDt() {
        return dt;
    }

    void setOldSolution(MixedFunction& uOld) {
        *oldSolution = *uOld;
    }

    virtual void updateStaticData(const Entity& e, VBlockType* sol) {
        return;
    }

    template<class TypeTag> void assembleBC(const Entity& e) {
        return;
    }

    // parameters given in constructor
    VertexMapper vertexMapper;
    StokesProblem<G,RT>& problem;
    MixedFunction& currentSolution;
    MixedFunction oldSolution;

public:
    double dt;
};
}
#endif
