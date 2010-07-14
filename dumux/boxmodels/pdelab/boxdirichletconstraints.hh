// $Id: boxdirichletconstraints.hh 3834 2010-07-14 12:50:32Z bernd $
// -*- tab-width: 4; indent-tabs-mode: nil -*-
/*****************************************************************************
 *   Copyright (C) 2009-2010 by Andreas Lauser                               *
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
#ifndef DUMUX_BOXDIRICHLETCONSTRAINTS_HH
#define DUMUX_BOXDIRICHLETCONSTRAINTS_HH

#include <cstddef>

#include<dune/common/exceptions.hh>
#include<dune/grid/common/genericreferenceelements.hh>
#include<dune/grid/common/grid.hh>
#include<dune/common/geometrytype.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/finiteelementmap/conformingconstraints.hh>

#include<dumux/common/boundarytypes.hh>
#include<dumux/common/boundaryconditions.hh>

namespace Dumux {
//! Constraints construction
// works in any dimension and on all element types
template <class TypeTag>
class BoxDirichletConstraints // : public Dune::PDELab::ConformingDirichletConstraints
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))  Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))   Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;

    enum {numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq))};
    enum {dim = GridView::dimension};

    typedef Dumux::BoundaryTypes<numEq> BoundaryTypeVector;
    Problem &problem_;

public:
    enum { doBoundary = true };
    enum { doProcessor = false };
    enum { doSkeleton = false };
    enum { doVolume = false };

    BoxDirichletConstraints(Problem& problem)
        : problem_(problem)
    {}

    //! boundary constraints
    /**
     * \tparam F   grid function returning boundary condition type
     * \tparam IG  intersection geometry
     * \tparam LFS local function space
     * \tparam T   TransformationType
     */
    template<typename F, typename I, typename LFS,
             typename T>
    void boundary (const F& f,
                   const Dune::PDELab::IntersectionGeometry<I>& ig,
                   const LFS& lfs,
                   T& trafo) const
    {
        const ElementPointer& elementPointer = ig.inside();
        FVElementGeometry fvElemGeom(problem_.gridView());

        fvElemGeom.update(*elementPointer);
        BoundaryTypeVector bcTypes;

        //problem_.boundaryTypes(values, element, fvElemGeom, ig.intersection(), scvIdx, boundaryFaceIdx);
        //typename F::Traits::RangeType bctype;

        const int face = ig.indexInInside();

        // find all local indices of this face
        Dune::GeometryType gt = ig.inside()->type();
        typedef typename Dune::PDELab::IntersectionGeometry<I>::ctype DT;
        const int dim = Dune::PDELab::IntersectionGeometry<I>::Entity::Geometry::dimension;
        const Dune::GenericReferenceElement<DT,dim>& refelem = Dune::GenericReferenceElements<DT,dim>::general(gt);

        // empty map means Dirichlet constraint
        typename T::RowType empty;

        for (int faceVertIdx = 0; faceVertIdx < refelem.size(face, 1, dim); faceVertIdx ++){
            int elemVertIdx = refelem.subEntity(face, 1, faceVertIdx, dim);
            int boundaryFaceIdx = fvElemGeom.boundaryFaceIndex(face, faceVertIdx);

            bcTypes.reset();
            problem_.boundaryTypes(bcTypes, *elementPointer, fvElemGeom, ig.intersection(), elemVertIdx, boundaryFaceIdx);
            bcTypes.checkWellPosed();

            for (std::size_t i = 0; i < lfs.localFiniteElement().localCoefficients().size(); i++) {
                // The codim to which this dof is attached to
                unsigned int codim = lfs.localFiniteElement().localCoefficients().localKey(i).codim();

                if (codim!=dim)
                    continue;
                if (lfs.localFiniteElement().localCoefficients().localKey(i).subEntity() != elemVertIdx)
                    continue;

                if (bcTypes.isDirichlet(F::eqIdx)) {
                    trafo[i] = empty;
                }
            }
        }
    }
};

// extend constraints class by processor boundary
template <class TypeTag>
class NonoverlappingBoxDirichletConstraints : public BoxDirichletConstraints<TypeTag>
{
public:
  enum { doVolume = true };
  typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))  Problem;

  NonoverlappingBoxDirichletConstraints(Problem& problem)
  : BoxDirichletConstraints<TypeTag>(problem)
  {}


  template<typename E, typename LFS, typename T>
  void volume (const Dune::PDELab::ElementGeometry<E>& eg, const LFS& lfs, T& trafo) const
  {
    // nothing to do for interior entities
    if (eg.entity().partitionType()==Dune::InteriorEntity)
      return;

    // empty map means Dirichlet constraint
    typename T::RowType empty;

    typedef typename LFS::Traits::GridFunctionSpaceType::Traits::BackendType B;

    // loop over all degrees of freedom and check if it is not owned by this processor
    for (size_t i=0; i<lfs.localFiniteElement().localCoefficients().size(); i++)
      {
        if (ghost_[lfs.globalIndex(i)]!=0)
          {
            trafo[i] = empty;
          }
      }
  }

  template<class GFS>
  void compute_ghosts (const GFS& gfs)
  {
    typedef typename GFS::template VectorContainer<int>::Type V;
    V ighost(gfs);
    Dune::PDELab::GhostDataHandle<GFS,V> gdh(gfs,ighost);
    if (gfs.gridview().comm().size()>1)
      gfs.gridview().communicate(gdh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
    ighost.std_copy_to(ghost_);
    rank = gfs.gridview().comm().rank();
  }

  void print ()
  {
    std::cout << "/" << rank << "/ " << "ghost size="
              << ghost_.size() << std::endl;
    for (std::size_t i=0; i<ghost_.size(); i++)
      std::cout << "/" << rank << "/ " << "ghost[" << i << "]="
                << ghost_[i] << std::endl;
  }

private:
  int rank;
  std::vector<int> ghost_;
};

// extend constraints class by processor boundary
template <class TypeTag>
class OverlappingBoxDirichletConstraints : public BoxDirichletConstraints<TypeTag>
{
public:
  enum { doProcessor = true };
  typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem))  Problem;

  OverlappingBoxDirichletConstraints(Problem& problem)
  : BoxDirichletConstraints<TypeTag>(problem)
  {}

  // boundary constraints
  // IG : intersection geometry
  // LFS : local function space
  // T : TransformationType
  template<typename I, typename LFS, typename T>
  void processor (const Dune::PDELab::IntersectionGeometry<I>& ig,
                  const LFS& lfs, T& trafo) const
  {
    // determine face
    const int face = ig.indexInInside();

    // find all local indices of this face
    Dune::GeometryType gt = ig.inside()->type();
    typedef typename Dune::PDELab::IntersectionGeometry<I>::ctype DT;
    const int dim = Dune::PDELab::IntersectionGeometry<I>::Entity::Geometry::dimension;


    const Dune::GenericReferenceElement<DT,dim>& refelem = Dune::GenericReferenceElements<DT,dim>::general(gt);

    // empty map means Dirichlet constraint
    typename T::RowType empty;

    // loop over all degrees of freedom and check if it is on given face
    for (size_t i=0; i<lfs.localFiniteElement().localCoefficients().size(); i++)
      {
        // The codim to which this dof is attached to
        unsigned int codim = lfs.localFiniteElement().localCoefficients().localKey(i).codim();

        if (codim==0) continue;

        for (int j=0; j<refelem.size(face,1,codim); j++)
          if (lfs.localFiniteElement().localCoefficients().localKey(i).subEntity()==refelem.subEntity(face,1,j,codim))
            trafo[i] = empty;
      }
  }
};




}

#endif
