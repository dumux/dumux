#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <fstream>
#include<dune/grid/sgrid.hh> // load sgrid definition
#include<dune/grid/common/gridinfo.hh> // definition of gridinfo /*@\label{gs:inc1}@*/
#include<dune/grid/common/mcmgmapper.hh>
#include<dune/istl/bcrsmatrix.hh>
#include<dune/istl/io.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/grid/common/quadraturerules.hh>
#include<dune/istl/ilu.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/grid/utility/hierarchicsearch.hh>
#include "dumux/pardiso/pardiso.hh"
#include "isolatepipe/isolate.hh"
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

template<class G, class V>
void vtkout (const G& grid, const V& Sat,const V& press,const V& perm, char* name, int k);

//void inline TimeloopOptsPipe( double& tstart, double& tend, double& max_dt, double& first_dt, double& CFL_factor, int& flag,
//              int& n_iter, double& max_def, int& modulo, int& stages );

void inline TimeloopOptsPorous(double& tstart, double& tend, double& max_dt, double& max_def, int& modulo);

//#include"vtkout.hh" // Output
#include"bc_ic.hh"
//#include"dumux/material/properties.hh"
#include"pipeflow_incomp_newcouple.hh"
#include"pipephysic.hh"

//! Parameter for mapper class
template<int dim>
struct P0Layout
{
    bool contains (Dune::GeometryType gt)
    {
        if (gt.dim()==dim) return true;
        return false;
    }
};

//! Parameter for mapper class
template<int dim>
struct PDimLayout
{
    bool contains (Dune::GeometryType gt)
    {
        if (gt.dim()==0) return true;
        return false;
    }
};


