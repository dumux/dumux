// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
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
#ifndef DUNE_VARIABLECLASS2P_NEW_HH
#define DUNE_VARIABLECLASS2P_NEW_HH

#include <dune/istl/bvector.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

/**
 * @file
 * @brief  Class including the variables and data of discretized data of the constitutive relations
 * @author Markus Wolff
 */

namespace Dune
{
/*!
 * \ingroup fracflow
 * \ingroup diffusion
 * \ingroup transport
 */
 //! Class including the variables and data of discretized data of the constitutive relations.
 /*! The variables of two-phase flow, which are one pressure and one saturation are stored in this class.
 * Additionally, a velocity needed in the transport part of the decoupled two-phase flow is stored, as well as discretized data of constitutive relationships like
 * mobilities, fractional flow functions and capillary pressure. Thus, they have to be callculated just once in every time step or every iteration step.
 *
 * Template parameters are:

 - GridView      a DUNE gridview type
 - Scalar        type used for scalar quantities
 */
template<class GridView, class Scalar>
class VariableClass
{
private:
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

typedef    typename GridView::Grid Grid;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

public:
    typedef Dune::BlockVector< Dune::FieldVector<Scalar,1> > ScalarVectorType;//!<type for vector of scalars
    typedef Dune::BlockVector< FieldVector<FieldVector<Scalar, 1>, 2*dim> > PotType;//!<type for vector of vectors (of size 2 x dimension) of scalars
    typedef Dune::BlockVector< FieldVector<FieldVector<Scalar, dim>, 2*dim> > VelType;//!<type for vector of vectors (of size 2 x dimension) of vector (of size dimension) of scalars

private:
    GridView& gridViewDiffusion_;
    GridView& gridViewTransport_;
    const IndexSet& indexSetDiffusion_;
    const IndexSet& indexSetTransport_;
    const int gridSizeDiffusion_;
    const int gridSizeTransport_;

    ScalarVectorType saturation_;
    ScalarVectorType pressure_;
    ScalarVectorType mobilityWetting_;//store lambda for efficiency reasons
    ScalarVectorType mobilityNonWetting_;
    ScalarVectorType fracFlowFuncWetting_;
    ScalarVectorType fracFlowFuncNonWetting_;
    ScalarVectorType capillaryPressure_;
    VelType velocity_;
    PotType potentialWetting_;
    PotType potentialNonWetting_;

    bool multiscale_;
    const int codim_;
public:
    //! Constructs a VariableClass object
    /**
     *  @param gridViewDiff a DUNE gridview object corresponding to the diffusion equation
     *  @param gridViewTrans a DUNE gridview object corresponding to the transport equation
     *  @param initialSat initial value for the saturation (only necessary if only diffusion part is solved)
     *  @param initialVel initial value for the velocity (only necessary if only transport part is solved)
     */
    VariableClass(GridView& gridViewDiff, GridView& gridViewTrans, Scalar& initialSat = *(new Scalar(0)), Dune::FieldVector<Scalar, dim>& initialVel = *(new Dune::FieldVector<Scalar, dim> (0)))
    : gridViewDiffusion_(gridViewDiff), gridViewTransport_(gridViewTrans),
    indexSetDiffusion_(gridViewDiff.indexSet()),indexSetTransport_(gridViewTrans.indexSet()),
    gridSizeDiffusion_(indexSetDiffusion_.size(0)),gridSizeTransport_(indexSetTransport_.size(0)), multiscale_(true), codim_(0)
    {
        //resize to grid size
        pressure_.resize(gridSizeDiffusion_);
        velocity_.resize(gridSizeDiffusion_);
        potentialWetting_.resize(gridSizeDiffusion_);
        potentialNonWetting_.resize(gridSizeDiffusion_);
        saturation_.resize(gridSizeTransport_);
        mobilityWetting_.resize(gridSizeTransport_);//lambda is dependent on saturation! ->choose same size
        mobilityNonWetting_.resize(gridSizeTransport_);
        fracFlowFuncWetting_.resize(gridSizeTransport_);
        fracFlowFuncNonWetting_.resize(gridSizeTransport_);
        capillaryPressure_.resize(gridSizeTransport_);

        //initialise variables
        pressure_ = 0;
        velocity_ = initialVel;
        saturation_ = initialSat;
        mobilityWetting_ = 0;
        mobilityNonWetting_ = 0;
        fracFlowFuncWetting_ = 0;
        fracFlowFuncNonWetting_ = 0;
        capillaryPressure_ = 0;
        initializePotentials(initialVel);
    }

    //! Constructs a VariableClass object
    /**
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     *  @param initialSat initial value for the saturation (only necessary if only diffusion part is solved)
     *  @param initialVel initial value for the velocity (only necessary if only transport part is solved)
     */
    VariableClass(GridView& gridView, Scalar& initialSat = *(new Scalar(0)), Dune::FieldVector<Scalar, dim>& initialVel = *(new Dune::FieldVector<Scalar, dim> (0)))
    : gridViewDiffusion_(gridView), gridViewTransport_(gridView),
    indexSetDiffusion_(gridView.indexSet()),indexSetTransport_(gridView.indexSet()),
    gridSizeDiffusion_(indexSetDiffusion_.size(0)),gridSizeTransport_(indexSetTransport_.size(0)), multiscale_(false), codim_(0)
    {
        //resize to grid size
        pressure_.resize(gridSizeDiffusion_);
        velocity_.resize(gridSizeDiffusion_);
        potentialWetting_.resize(gridSizeDiffusion_);
        potentialNonWetting_.resize(gridSizeDiffusion_);
        saturation_.resize(gridSizeTransport_);
        mobilityWetting_.resize(gridSizeTransport_);//lambda is dependent on saturation! ->choose same size
        mobilityNonWetting_.resize(gridSizeTransport_);
        fracFlowFuncWetting_.resize(gridSizeTransport_);
        fracFlowFuncNonWetting_.resize(gridSizeTransport_);
        capillaryPressure_.resize(gridSizeTransport_);

        //initialise variables
        pressure_ = 0;
        velocity_ = initialVel;
        saturation_ = initialSat;
        mobilityWetting_ = 0;
        mobilityNonWetting_ = 0;
        fracFlowFuncWetting_ = 0;
        fracFlowFuncNonWetting_ = 0;
        capillaryPressure_ = 0;
        initializePotentials(initialVel);
    }

    //! Constructs a VariableClass object
    /**
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     *  @param codim codimension of the entity of which data has to be strored
     *  @param initialSat initial value for the saturation (only necessary if only diffusion part is solved)
     *  @param initialVel initial value for the velocity (only necessary if only transport part is solved)
     */
    VariableClass(GridView& gridView, int codim, Scalar& initialSat = *(new Scalar(0)), Dune::FieldVector<Scalar, dim>& initialVel = *(new Dune::FieldVector<Scalar, dim> (0)))
    : gridViewDiffusion_(gridView), gridViewTransport_(gridView),
    indexSetDiffusion_(gridView.indexSet()),indexSetTransport_(gridView.indexSet()),
    gridSizeDiffusion_(indexSetDiffusion_.size(codim)),gridSizeTransport_(indexSetTransport_.size(codim)), multiscale_(false), codim_(codim)
    {
        //resize to grid size
        pressure_.resize(gridSizeDiffusion_);
        velocity_.resize(gridSizeDiffusion_);
        potentialWetting_.resize(gridSizeDiffusion_);
        potentialNonWetting_.resize(gridSizeDiffusion_);
        saturation_.resize(gridSizeTransport_);
        mobilityWetting_.resize(gridSizeTransport_);//lambda is dependent on saturation! ->choose same size
        mobilityNonWetting_.resize(gridSizeTransport_);
        fracFlowFuncWetting_.resize(gridSizeTransport_);
        fracFlowFuncNonWetting_.resize(gridSizeTransport_);
        capillaryPressure_.resize(gridSizeTransport_);

        //initialise variables
        pressure_ = 0;
        velocity_ = initialVel;
        saturation_ = initialSat;
        mobilityWetting_ = 0;
        mobilityNonWetting_ = 0;
        fracFlowFuncWetting_ = 0;
        fracFlowFuncNonWetting_ = 0;
        capillaryPressure_ = 0;
        initializePotentials(initialVel);
    }
private:
    void initializePotentials (Dune::FieldVector<Scalar, dim>& initialVel)
    {
        if (initialVel.two_norm())
        {
            // compute update vector
            ElementIterator eItEnd = gridViewTransport_.template end<0>();
            for (ElementIterator eIt = gridViewTransport_.template begin<0>(); eIt != eItEnd; ++eIt)
            {
                // cell index
                int globalIdxI = indexSetTransport_.index(*eIt);

                // run through all intersections with neighbors and boundary
                IntersectionIterator
                isItEnd = gridViewTransport_.template iend(*eIt);
                for (IntersectionIterator
                        isIt = gridViewTransport_.template ibegin(*eIt); isIt
                        !=isItEnd; ++isIt)
                {
                    // local number of facet
                    int indexInInside = isIt->indexInInside();

                    // get geometry type of face
                    Dune::GeometryType faceGT = isIt->geometryInInside().type();

                    // center in face's reference element
                    const Dune::FieldVector<Scalar,dim-1>&
                    faceLocal = Dune::ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

                    Dune::FieldVector<Scalar,dimWorld> unitOuterNormal = isIt->unitOuterNormal(faceLocal);

                    potentialWetting_[globalIdxI][indexInInside] = initialVel*unitOuterNormal;
                    potentialNonWetting_[globalIdxI][indexInInside] = initialVel*unitOuterNormal;
                }
            }
        }
        else
        {
            potentialWetting_ = 0;
            potentialNonWetting_ = 0;
        }
        return;
    }
public:
    //! Return saturation vector
    ScalarVectorType& saturation()
    {
        return saturation_;
    }

    //! Return pressure vector
    ScalarVectorType& pressure()
    {
        return pressure_;
    }

    //! Return velocity vector
    VelType& velocity()
    {
        return velocity_;
    }

    //! Return vector of wetting phase potential gradients
    PotType& potentialWetting()
    {
        return potentialWetting_;
    }

    //! Return vector of non-wetting phase potential gradients
    PotType& potentialNonWetting()
    {
        return potentialNonWetting_;
    }

    //! Return vector of wetting phase mobilities
    ScalarVectorType& mobilityWetting()
    {
        return mobilityWetting_;
    }

    //! Return vector of non-wetting phase mobilities
    ScalarVectorType& mobilityNonWetting()
    {
        return mobilityNonWetting_;
    }

    //! Return vector of wetting phase fractional flow functions
    ScalarVectorType& fracFlowFuncWetting()
    {
        return fracFlowFuncWetting_;
    }

    //! Return vector of non-wetting phase fractional flow functions
    ScalarVectorType& fracFlowFuncNonWetting()
    {
        return fracFlowFuncNonWetting_;
    }

    //! Return capillary pressure vector
    ScalarVectorType& capillaryPressure()
    {
        return capillaryPressure_;
    }

    //! Get index of element (codim 0 entity) corresponding to the grid of the discretized diffusion equation.
    /*! Get index of element (codim 0 entity) corresponding to the grid of the discretized diffusion equation.
     * @param element codim 0 entity
     * \return element index
     */
    int indexDiffusion(const Element& element)
    {
        return indexSetDiffusion_.index(element);
    }

    //! Get index of element (codim 0 entity) corresponding to the grid of the discretized transport equation.
    /*! Get index of element (codim 0 entity) corresponding to the grid of the discretized transport equation.
     * @param element codim 0 entity
     * \return element index
     */
    int indexTransport(const Element& element)
    {
        return indexSetTransport_.index(element);
    }

    //!Return the number of data elements of the discretized diffusion equation
    int gridSizeDiffusion()
    {
        return gridSizeDiffusion_;
    }

    //!Return the number of data elements of the discretized transport equation
    int gridSizeTransport()
    {
        return gridSizeTransport_;
    }

    //!Return gridView on the grid of the discretized diffusion equation
    GridView& gridViewDiffusion()
    {
        return gridViewDiffusion_;
    }

    //!Return gridView on the grid of the discretized transport equation
    GridView& gridViewTransport()
    {
        return gridViewTransport_;
    }

    //! Get saturation
    /*! evaluate saturation at given element
     @param  element      entity of codim 0
     \return     value of saturation
     */
    const Dune::FieldVector<Scalar,1>& satElement(const Element& element) const
    {
        return saturation_[indexSetTransport_.index(element)];;
    }

    //! Get pressure
    /*! evaluate pressure at given element
     @param  element      entity of codim 0
     \return     value of pressure
     */
    const Dune::FieldVector<Scalar,1>& pressElement(const Element& element) const
    {
        return pressure_[indexSetDiffusion_.index(element)];
    }

    //! Get velocity at given element face
    /*! evaluate velocity at given location
     @param  element      entity of codim 0
     @param  indexInInside     index in reference element
     \return     vector of velocity
     */
    const Dune::FieldVector<Scalar,dim>& vTotalElementFace(const Element& element,
            const int indexInInside) const
    {
        int elemId = indexSetTransport_.index(element);

        return (velocity_[elemId][indexInInside]);
    }

    //! \brief Write data files
    /*!
     *  \param name file name
     *  \param k format parameter
     */
    void vtkout(const char* name, int k)
    {
        vtkoutMultiLevel(name, k);
    }

    //! \brief Write saturation and pressure into file
    /*!
     *  \param name file name
     *  \param k format parameter
     */
    void vtkoutMultiLevel(const char* name, int k) const
    {
        if (codim_ == 0)
        {
            if (multiscale_)
            {
                Dune::VTKWriter<GridView> vtkwriterpressure(gridViewDiffusion_);
                char fname[128];
                sprintf(fname, "%s-press%05d", name, k);
                vtkwriterpressure.addCellData(pressure_, "pressure");
                vtkwriterpressure.write(fname, Dune::VTKOptions::ascii);

                Dune::VTKWriter<GridView> vtkwritersaturation(gridViewTransport_);
                sprintf(fname, "%s-%05d", name, k);
                vtkwritersaturation.addCellData(saturation_, "saturation");
                vtkwritersaturation.write(fname, VTKOptions::ascii);
            }
            else
            {
                VTKWriter<GridView> vtkwriter(gridViewDiffusion_);
                char fname[128];
                sprintf(fname, "%s-%05d", name, k);
                vtkwriter.addCellData(saturation_, "saturation");
                vtkwriter.addCellData(pressure_, "pressure");
                vtkwriter.write(fname, VTKOptions::ascii);
            }
        }
        if (codim_ == dim)
        {
            if (multiscale_)
            {
                Dune::VTKWriter<GridView> vtkwriterpressure(gridViewDiffusion_);
                char fname[128];
                sprintf(fname, "%s-press%05d", name, k);
                vtkwriterpressure.addVertexData(pressure_, "pressure");
                vtkwriterpressure.write(fname, Dune::VTKOptions::ascii);

                Dune::VTKWriter<GridView> vtkwritersaturation(gridViewTransport_);
                sprintf(fname, "%s-%05d", name, k);
                vtkwritersaturation.addVertexData(saturation_, "saturation");
                vtkwritersaturation.write(fname, VTKOptions::ascii);
            }
            else
            {
                VTKWriter<GridView> vtkwriter(gridViewDiffusion_);
                char fname[128];
                sprintf(fname, "%s-%05d", name, k);
                vtkwriter.addVertexData(saturation_, "saturation");
                vtkwriter.addVertexData(pressure_, "pressure");
                vtkwriter.write(fname, VTKOptions::ascii);
            }
        }

        return;
    }
};
}
#endif
