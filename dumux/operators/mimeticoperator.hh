// $Id$

#ifndef DUNE_MIMETICOPERATOR_HH
#define DUNE_MIMETICOPERATOR_HH

#include<iostream>
#include<vector>
#include<set>
#include<map>
#include<stdio.h>
#include<stdlib.h>

#include<dune/common/timer.hh>
#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/mcmgmapper.hh>

#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/bcrsmatrix.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include"dumux/functions/CRfunction.hh"
#include"dumux/shapefunctions/CRshapefunctions.hh"
#include"localstiffnessext.hh"
#include"CRoperator.hh"

namespace Dune
{
/*! @brief Levelwise assembler

  This class serves as a base class for local assemblers. It provides
  space and access to the local stiffness matrix. The actual assembling is done
  in a derived class via the virtual assemble method.

  The template parameters are:

  - G    A grid type
  - RT   The field type used in the elements of the stiffness matrix
  - m    number of degrees of freedom per node (system size)
*/
template<class G, class RT, int m=1>
class MimeticOperatorAssembler : public LevelCROperatorAssembler<G, RT, m>
{
    template<int dim>
    struct ElementLayout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };

    enum {n=G::dimension};
    typedef typename G::ctype DT;
    typedef typename G::LevelGridView GV;
    typedef typename GV::IndexSet IS;
    typedef LevelCommunicate<G> LC;
    typedef LevelP0Function<G,RT,2*n> VType;
    typedef BlockVector< FieldVector<RT,1> > PType;
    typedef typename GV::template Codim<0>::Iterator Iterator;
    typedef MultipleCodimMultipleGeomTypeMapper<GV,ElementLayout> EM;

public:

    MimeticOperatorAssembler (const G& grid, int level)
        : LevelCROperatorAssembler<G, RT, m>(grid, level), elementmapper(grid.levelView(level))
    {}

    template<class I>
    void calculatePressure (LocalStiffnessExt<I,G,RT,m>& loc, CRFunction<G,RT,GV,LC,m>& u,
                            VType& velocity, PType& pressure)
    {
        // run over all level elements
        Iterator eendit = this->grid.template lend<0>(0);
        for (Iterator it = this->grid.template lbegin<0>(0); it!=eendit; ++it)
        {
            // get access to shape functions for CR elements
            Dune::GeometryType gt = it->geometry().type();
            const typename Dune::CRShapeFunctionSetContainer<DT,RT,n>::value_type&
                sfs = Dune::CRShapeFunctions<DT,RT,n>::general(gt,1);

            int elemId = elementmapper.map(*it);

            // get local to global id map and pressure traces
            Dune::FieldVector<DT,2*n> pressTrace(0);
            for (int k = 0; k < sfs.size(); k++)
            {
                pressTrace[k] = (*u)[this->facemapper.template map<1>(*it, k)];
            }

            // The notation is borrowed from Aarnes/Krogstadt/Lie 2006, Section 3.4.
            // The matrix W developed here corresponds to one element-associated
            // block of the matrix B^{-1} there.
            Dune::FieldVector<DT,2*n> faceVol(0);
            Dune::FieldMatrix<DT,2*n,2*n> W(0);
            Dune::FieldVector<DT,2*n> c(0);
            Dune::FieldMatrix<DT,2*n,2*n> Pi(0);
            Dune::FieldVector<RT,2*n> F(0);
            RT dinv = 0;
            RT qmean = 0;
            loc.assembleElementMatrices(*it, faceVol, W, c, Pi, dinv, F, qmean);

            pressure[elemId] = dinv*(qmean + (F*pressTrace));

            Dune::FieldVector<RT,2*n> v(0);
            for (int i = 0; i < 2*n; i++)
                for (int j = 0; j < 2*n; j++)
                    v[i] += W[i][j]*faceVol[j]*(pressure[elemId] - pressTrace[j]);

            (*velocity)[elemId] = v;
        }
    }

private:
    EM elementmapper;
};
}
#endif
