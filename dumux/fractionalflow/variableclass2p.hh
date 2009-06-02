// $Id$
#ifndef DUNE_VARIABLECLASS2P_NEW_HH
#define DUNE_VARIABLECLASS2P_NEW_HH

#include <dune/istl/bvector.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

/**
 * @file
 * @brief  class including the variables
 * @author Markus Wolff
 */

namespace Dune
{
/** \todo Please doc me! */

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
    typedef Dune::BlockVector< Dune::FieldVector<Scalar,1> > ScalarVectorType;
    typedef Dune::BlockVector< FieldVector<FieldVector<Scalar, 1>, 2*dim> > PotType;
    typedef Dune::BlockVector< FieldVector<FieldVector<Scalar, dim>, 2*dim> > VelType;

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
public:

    VariableClass(GridView& gridViewDiff, GridView& gridViewTrans, Scalar& initialSat = *(new Scalar(0)), Dune::FieldVector<Scalar, dim>& initialVel = *(new Dune::FieldVector<Scalar, dim> (0)))
    : gridViewDiffusion_(gridViewDiff), gridViewTransport_(gridViewTrans),
    indexSetDiffusion_(gridViewDiff.indexSet()),indexSetTransport_(gridViewTrans.indexSet()),
    gridSizeDiffusion_(indexSetDiffusion_.size(0)),gridSizeTransport_(indexSetTransport_.size(0)), multiscale_(true)
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
    VariableClass(GridView& gridView, Scalar& initialSat = *(new Scalar(0)), Dune::FieldVector<Scalar, dim>& initialVel = *(new Dune::FieldVector<Scalar, dim> (0)))
    : gridViewDiffusion_(gridView), gridViewTransport_(gridView),
    indexSetDiffusion_(gridView.indexSet()),indexSetTransport_(gridView.indexSet()),
    gridSizeDiffusion_(indexSetDiffusion_.size(0)),gridSizeTransport_(indexSetTransport_.size(0)), multiscale_(false)
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

    ScalarVectorType& saturation()
    {
        return saturation_;
    }
    ScalarVectorType& pressure()
    {
        return pressure_;
    }
    VelType& velocity()
    {
        return velocity_;
    }
    PotType& potentialWetting()
    {
        return potentialWetting_;
    }
    PotType& potentialNonWetting()
    {
        return potentialNonWetting_;
    }
    ScalarVectorType& mobilityWetting()
    {
        return mobilityWetting_;
    }
    ScalarVectorType& mobilityNonWetting()
    {
        return mobilityNonWetting_;
    }
    ScalarVectorType& fracFlowFuncWetting()
    {
        return fracFlowFuncWetting_;
    }
    ScalarVectorType& fracFlowFuncNonWetting()
    {
        return fracFlowFuncNonWetting_;
    }
    ScalarVectorType& capillaryPressure()
    {
        return capillaryPressure_;
    }

    int indexDiffusion(const Element& element)
    {
        return indexSetDiffusion_.index(element);
    }
    int indexTransport(const Element& element)
    {
        return indexSetTransport_.index(element);
    }
    int gridSizeDiffusion()
    {
        return gridSizeDiffusion_;
    }
    int gridSizeTransport()
    {
        return gridSizeTransport_;
    }
    GridView& gridViewDiffusion()
    {
        return gridViewDiffusion_;
    }
    GridView& gridViewTransport()
    {
        return gridViewTransport_;
    }

    const Dune::FieldVector<Scalar,1>& satElement(const GlobalPosition globalPos,
            const Element& element, const LocalPosition localPos) const
    {
        return saturation_[indexSetTransport_.index(element)];;
    }

    const Dune::FieldVector<Scalar,1>& pressElement(const GlobalPosition globalPos,
            const Element& element, const LocalPosition localPos) const
    {
        return pressure_[indexSetDiffusion_.index(element)];
    }

    const Dune::FieldVector<Scalar,dim>& vTotalElementFace(const Element& element,
            const int numberInSelf) const
    {
        int elemId = indexSetTransport_.index(element);

        return (velocity_[elemId][numberInSelf]);
    }

    /*! @brief prints the saturation to a VTK file
     *
     *  The file name is "<name>-<k>.vtu" where k is an integer number.
     *  @param name specifies the name of the VTK file
     *  @param k specifies a number
     */

    void vtkout(const char* name, int k)
    {
        vtkoutMultiLevel(name, k);
    }

    void vtkoutMultiLevel(const char* name, int k) const
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

        return;
    }
};
}
#endif
