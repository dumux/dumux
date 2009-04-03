// $Id: variableclass.hh 1329 2009-03-02 15:14:45Z lauser $
#ifndef DUNE_VARIABLECLASS_HH
#define DUNE_VARIABLECLASS_HH

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

template<class Grid, class Scalar>
class VariableClass
{
    enum
        {
            dim = Grid::dimension, dimWorld = Grid::dimensionworld
        };

    typedef typename Grid::LevelGridView GridView;
    typedef    typename Grid::Traits::template Codim<0>::Entity Element;

    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

    typedef typename GridView::IndexSet IndexSet;

public:
    typedef Dune::BlockVector< Dune::FieldVector<Scalar,1> > ScalarVectorType;
    typedef Dune::BlockVector< FieldVector<FieldVector<Scalar, dim>, 2*dim> > VelType;

    Grid& grid;
    const int levelDiffusion;
    const int levelTransport;
    const IndexSet& indexSetDiffusion;
    const IndexSet& indexSetTransport;
    const int gridSizeDiffusion;
    const int gridSizeTransport;

    ScalarVectorType saturation;
    ScalarVectorType pressure;
    ScalarVectorType mobilityWetting;//store lambda for efficiency reasons
    ScalarVectorType mobilityNonWetting;
    ScalarVectorType fracFlowFuncWetting;
    ScalarVectorType fracFlowFuncNonWetting;
    ScalarVectorType capillaryPressure;
    VelType velocity;

    VariableClass(Grid& grid, Scalar& initialSat = *(new Scalar(0)), Scalar& initialPress = *(new Scalar(0)), Dune::FieldVector<Scalar, dim>& initialVel = *(new Dune::FieldVector<Scalar, dim> (0)),int transLev = -1, int diffLev = -1)
        : grid(grid),
          levelDiffusion((diffLev >= 0) ? diffLev : grid.maxLevel()), levelTransport((transLev >= 0) ? transLev : grid.maxLevel()),
          indexSetDiffusion(grid.levelView(levelDiffusion).indexSet()),indexSetTransport(grid.levelView(levelTransport).indexSet()),
          gridSizeDiffusion(indexSetDiffusion.size(0)),gridSizeTransport(indexSetTransport.size(0))
    {
        //resize to grid size
        pressure.resize(gridSizeDiffusion);
        velocity.resize(gridSizeDiffusion);
        saturation.resize(gridSizeTransport);
        mobilityWetting.resize(gridSizeTransport);//lambda is dependent on saturation! ->choose same size
        mobilityNonWetting.resize(gridSizeTransport);
        fracFlowFuncWetting.resize(gridSizeTransport);
        fracFlowFuncNonWetting.resize(gridSizeTransport);
        capillaryPressure.resize(gridSizeTransport);

        //initialise variables
        pressure = initialPress;
        velocity = initialVel;
        saturation = initialSat;
        mobilityWetting = 0;
        mobilityNonWetting = 0;
        fracFlowFuncWetting = 0;
        fracFlowFuncNonWetting = 0;
        capillaryPressure = 0;
    }

    const Dune::FieldVector<Scalar,1>& sat(const GlobalPosition globalPos,
                                           const Element& element, const LocalPosition localPos) const
    {
        return saturation[indexSetTransport.index(element)];;
    }

    const Dune::FieldVector<Scalar,1>& press(const GlobalPosition globalPos,
                                             const Element& element, const LocalPosition localPos) const
    {
        return pressure[indexSetDiffusion.index(element)];
    }

    const Dune::FieldVector<Scalar,dim>& vTotal(const Element& element,
                                                const int numberInSelf) const
    {
        int elemId = indexSetTransport.index(element);

        return (velocity[elemId][numberInSelf]);
    }

    /*! @brief prints the saturation to a VTK file
     *
     *  The file name is "<name>-<k>.vtu" where k is an integer number.
     *  @param name specifies the name of the VTK file
     *  @param k specifies a number
     */

    void vtkout (const char* name, int k)
    {
        vtkoutMultiLevel(name,k,levelDiffusion,levelTransport);
    }

    void vtkoutMultiLevel(const char* name, int k, int pressureLevel = 0, int satLevel = 0) const
    {
        if (pressureLevel == satLevel)
        {
            VTKWriter<typename Grid::LeafGridView> vtkwriter(grid.leafView());
            char fname[128];
            sprintf(fname, "%s-%05d", name, k);
            vtkwriter.addCellData(saturation, "saturation");
            vtkwriter.addCellData(pressure, "total pressure p~");
            vtkwriter.write(fname, VTKOptions::ascii);
        }
        else
        {
            Dune::VTKWriter<typename Grid::LevelGridView>
                vtkwriterpressure(grid.levelView(pressureLevel));
            char fname[128];
            sprintf(fname, "%s-press%05d", name, k);
            vtkwriterpressure.addCellData(pressure, "total pressure p~");
            vtkwriterpressure.write(fname, Dune::VTKOptions::ascii);

            Dune::VTKWriter<typename Grid::LevelGridView>
                vtkwritersaturation(grid.levelView(satLevel));
            sprintf(fname, "%s-%05d", name, k);
            vtkwritersaturation.addCellData(saturation, "saturation");
            vtkwritersaturation.write(fname, VTKOptions::ascii);
        }
        return;
    }
};
}
#endif
