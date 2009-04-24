// $Id$

#ifndef DUNE_EXSOLUTION_HH
#define DUNE_EXSOLUTION_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>

#include<dumux/material/twophaserelations_deprecated.hh>
#include<dumux/material/linearlaw_deprecated.hh>

/**
 * @file
 * @brief  Basic class for embedding of an exact solution
 * @author Markus Wolff
 */

namespace Dune
{
/** \todo Please doc me! */

template<class G, class RT> class ExSolution
{
    typedef typename G::ctype DT;
    enum
        {    n=G::dimension,m=2};
    typedef BlockVector<FieldVector<RT, 1> > BV;
    typedef BlockVector<FieldVector<RT, m> > BVu;
    typedef BlockVector<FieldVector<RT, (m+1)> > BVuEx;

    typedef typename G::LevelGridView GV;
    typedef typename GV::IndexSet IS;
    typedef typename GV::template Codim<0>::Iterator Iterator;
    typedef typename G::Traits::template Codim<0>::Entity Entity;

    template<int dim> struct ElementLayout
    {
        bool contains(Dune::GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV,ElementLayout>
    ElementMapper;

    //private accsess functions
    void uExInit()
    {
        uEx.resize(mapper.size());
        uEx=0;
        return;
    }
protected:
    void uExIn(BV &uExIn)
    {
        uEx=uExIn;
        return;
    }

    void uExInVertex(RT &value, int &ElementIndex, int VariableIndex)
    {
        uEx[ElementIndex][VariableIndex]=value;
        return;
    }

public:

    ExSolution(const G &g, int lev = 0) :
        grid(g), mapper(grid.levelView(lev)), uEx(0), error(0),
        elementvolume(0)
    {
        uExInit();
    }

    virtual ~ExSolution()
    {
    }

    //public access function
    BVu& uExOut()
    {
        return uEx;
    }

    RT uExOutVertex(int &ElementIndex, int VariableIndex) const
    {
        return uEx[ElementIndex][VariableIndex];
    }

    void calcSatError(BV &Approx)
    {
        int size=mapper.size();
        error.resize(size);
        elementvolume.resize(size);

        error=0;
        elementvolume=0;

        const GV& gridview(grid.levelView(0));

        Iterator eendit = gridview.template end<0>();
        for (Iterator it = gridview.template begin<0>(); it
                 != eendit; ++it)
        {
            // get entity
            const Entity& entity = *it;

            int index = mapper.map(*it);

            elementvolume[index]= entity.geometry().volume();
            //                        std::cout<<"elementvolume = "<<elementvolume[index]<<std::endl;
        }

        double globalvolume = elementvolume.one_norm();
        //                std::cout<<"globalvolume = "<<globalvolume<<std::endl;

        for (int i=0; i<size; i++)
        {
            error[i]=uEx[i][0]-Approx[i];
            //            std::cout<<"error = "<<error[i]<<std::endl;
            //            std::cout<<"uEx = "<<uEx[i]<<std::endl;
            //            std::cout<<"Approx = "<<Approx[i]<<std::endl;
        }
        //        std::cout<<"error = "<<error<<std::endl;

        double diffNorm = error.two_norm();
        std::cout<<"diffNorm = "<<diffNorm<<std::endl;

        for (int i=0; i<size; i++)
            uEx[i][1] = diffNorm * pow((elementvolume[i]/globalvolume), 0.5);

        return;
    }

protected:
    const G &grid;
    ElementMapper mapper;
    BVu uEx;
    BV error;
    BV elementvolume;
};
}
#endif
