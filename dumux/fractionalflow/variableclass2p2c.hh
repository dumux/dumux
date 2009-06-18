// $Id$

/*****************************************************************************
* Copyright (C) 2009 by Jochen Fritz                                         *
* Institute of Hydraulic Engineering                                         *
* University of Stuttgart, Germany                                           *
* email: <givenname>.<name>@iws.uni-stuttgart.de                             *
*                                                                            *
* This program is free software; you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by       *
* the Free Software Foundation; either version 2 of the License, or          *
* (at your option) any later version, as long as this copyright notice       *
* is included in its original form.                                          *
*                                                                            *
* This program is distributed WITHOUT ANY WARRANTY.                          *
*****************************************************************************/

#ifndef DUNE_VARIABLECLASS2P2C_HH
#define DUNE_VARIABLECLASS2P2C_HH

namespace Dune {

/** \ingroup decoupled2p2c
 *  \brief container for the variables needed for decoupled 2p2c computations.
 *
 *  All required variables for compositional two-phase computations are provided and
 *  can be accessed via this class.
 *  Template arguments: GridView is a View of one of the Dune Grid implementations.
 *  Scalar is the desired type to represent Scalar values.
 *
 */
template<class GridView, class Scalar> class VariableClass2p2c
{
private:
    enum {n=GridView::dimension};
    typedef typename GridView::Grid::ctype ct;
    typedef Dune::BlockVector< Dune::FieldVector<Scalar,1> > ScalarType;
    typedef typename GridView::template Codim<0>::Entity Entity;

    GridView& gridview;
    const typename GridView::IndexSet& indexset;
    int size;

public:
    ScalarType saturation;
    ScalarType pressure;
    ScalarType totalConcentration;
    ScalarType wet_X1, nonwet_X1;
    ScalarType density_wet, density_nonwet;
    ScalarType mobility_wet, mobility_nonwet;
    ScalarType volErr;

    template<int dim> struct ElementLayout
    {
        bool contains(Dune::GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };

    //! Constructor
    /**
     * \param gv an object of template argument GridView
     */
    VariableClass2p2c(GridView& gv) :
        gridview(gv), indexset(gv.indexSet()), size(gv.size(0))
    {
        saturation.resize(size);
        pressure.resize(size);
        totalConcentration.resize(2*size);
        wet_X1.resize(size);
        nonwet_X1.resize(size);
        density_wet.resize(size);
        density_nonwet.resize(size);
        mobility_wet.resize(size);
        mobility_nonwet.resize(size);
        volErr.resize(size);

        volErr = 0;
    }

    //! returns a reference to the saturation vector
    ScalarType& sat() const
    {
        return saturation;
    }

    //! returns a reference to the pressure vector
    ScalarType& press() const
    {
        return pressure;
    }

    //! returns a reference to the saturation at a distinct location
    const Dune::FieldVector<Scalar,1>& sat(const Dune::FieldVector<ct,n>& x, const Entity& e,
                                       const Dune::FieldVector<ct,n>& xi) const
    {
        return saturation[indexset.index(e)];;
    }

    //! returns a reference to the pressure at a distinct location
    const Dune::FieldVector<Scalar,1>& press(const Dune::FieldVector<ct,n>& x, const Entity& e,
                                         const Dune::FieldVector<ct,n>& xi) const
    {
        return pressure[indexset.index(e)];
    }

    /*! @brief writes all variables to a VTK File
     *
     *  The file name is "<name>-<k>.vtu" where k is an integer number.
     *  @param name specifies the name of the VTK file
     *  @param k specifies a number
     */
    void vtkout(const char* name, int k) const
    {
        ScalarType C1, C2;
        C1.resize(size); C2.resize(size);
        for (int i = 0; i < size; i++)
        {
            C1[i] = totalConcentration[i];
            C2[i] = totalConcentration[i + size];
        }
        VTKWriter<GridView> vtkwriter(gridview);
        char fname[128];
        sprintf(fname, "%s-%05d", name, k);
        vtkwriter.addCellData(saturation, "saturation [-]");
        vtkwriter.addCellData(pressure, "pressure [Pa]");
        vtkwriter.addCellData(C1, "total concentration 1 [kg/m^3]");
        vtkwriter.addCellData(C2, "total concentration 2 [kg/m^3]");
        vtkwriter.addCellData(volErr, "volume error [%]");
        vtkwriter.addCellData(wet_X1, "Mass fraction 1 in wetting phase [-]");
        vtkwriter.addCellData(nonwet_X1, "Mass fraction 1 in non-wetting phase [-]");
        vtkwriter.addCellData(density_wet, "wetting phase density [kg/m^3]");
        vtkwriter.addCellData(density_nonwet, "non-wetting phase density [kg/m^3]");
        vtkwriter.addCellData(mobility_wet, "wetting phase mobility [m*s/kg]");
        vtkwriter.addCellData(mobility_nonwet, "non-wetting phase mobility [m*s/kg]");
        vtkwriter.write(fname, VTKOptions::ascii);
        return;
    }

    /*! @brief writes only the pressure to a VTK File
     *
     *  The file name is "<name>-<k>.vtu" where k is an integer number.
     *  @param name specifies the name of the VTK file
     *  @param k specifies a number
     */
    void vtkoutpressure(const char* name, int k) const {
        VTKWriter<GridView> vtkwriter(gridview);
        char fname[128];
        sprintf(fname, "%s-press%05d", name, k);
        vtkwriter.addCellData(pressure, "total pressure p~");
        vtkwriter.write(fname, VTKOptions::ascii);
    }
};
}
#endif
