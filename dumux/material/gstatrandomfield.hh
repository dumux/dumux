// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Material
 * \brief Creating random fields using gstat
 */
#ifndef DUMUX_MATERIAL_GSTAT_RANDOM_FIELD_HH
#define DUMUX_MATERIAL_GSTAT_RANDOM_FIELD_HH

#include <iostream>
#include <algorithm>
#include <string>
#include <cmath>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk.hh>

namespace Dumux {

/*!
 * \ingroup Material
 * \brief Creating random fields using gstat
 *
 * gstat is an open source software tool which can (among other things) generate
 * geostatistical random fields (see <a href="http://www.gstat.org">http://www.gstat.org</a>
 * or \cite Pebesma1998a).
 *
 * To use this class, execute the installexternal.py from your DuMuX root
 * directory or download, unpack and install the tarball from the gstat-website.
 * Then rerun cmake (in the second case set GSTAT_ROOT in your input file to the
 * path where gstat is installed).
 */
template<class GridView, class Scalar>
class GstatRandomField
{
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    using DataVector = std::vector<Scalar>;
    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;

public:
    // Add field types if you want to implement e.g. tensor permeabilities.
    enum FieldType { scalar, log10 };

    /*!
     * \brief Constructor
     *
     * \param gridView the used gridView
     * \param elementMapper	Maps elements of the given grid view
     */
    GstatRandomField(const GridView& gridView, const ElementMapper& elementMapper)
    : gridView_(gridView)
    , elementMapper_(elementMapper)
    , data_(gridView.size(0))
    {}

    /*!
     * \brief Creates a new field with random variables, if desired.
     * Otherwise creates a data field from already available data.
     * For the random field generation three files are necessary.
     *
     * A \a gstatControlFile in which all commands and in/output files for gstat are specified.
     * A \a gstatInputFile contains all coordinates (cell centers) of the grid, so that
     * gstat can perform its random realization. The filename must be same as in the gstatControlFile.
     * A \a gstatOutputFile in which gstat writes the random values to this file.
     * The filename must be the same as in the gstatControlFile.
     * \param fieldType
     * \param gstatControlFile name of control file for gstat
     * \param gstatInputFile name of input file for gstat
     * \param gstatOutputFile name of the gstat output file
     * \param createNew set true to create a new field
     */
    void create(const std::string& gstatControlFile,
                const std::string& gstatInputFile = "gstatInput.txt",
                const std::string& gstatOutputFile = "permeab.dat",
                FieldType fieldType = FieldType::scalar,
                bool createNew = true)
    {
        fieldType_ = fieldType;
        if (createNew)
        {
#if !DUMUX_HAVE_GSTAT
            DUNE_THROW(Dune::InvalidStateException, "Requested data field generation with gstat"
              << " but gstat was not found on your system. Set GSTAT_ROOT to the path where gstat "
              << " is installed and pass it to CMake, e.g. through an opts file.");
#else
            std::ofstream gstatInput(gstatInputFile);
            for (const auto& element : elements(gridView_))
                // get global coordinates of cell centers
                gstatInput << element.geometry().center() << std::endl;

            gstatInput.close();

            std::string syscom;
            syscom = GSTAT_EXECUTABLE;
            syscom += " ";
            syscom += gstatControlFile;

            if (!gstatInput.good())
            {
                DUNE_THROW(Dune::IOError, "Reading the gstat control file: "
                             << gstatControlFile << " failed." << std::endl);
            }

            if (system(syscom.c_str()))
            {
                DUNE_THROW(Dune::IOError, "Executing gstat failed.");
            }
#endif
        }

        std::ifstream gstatOutput(gstatOutputFile);
        if (!gstatOutput.good())
        {
            DUNE_THROW(Dune::IOError, "Reading from file: "
                         << gstatOutputFile << " failed." << std::endl);
        }

        std::string line;
        std::getline(gstatOutput, line);

        Scalar trash, dataValue;
        for (const auto& element : elements(gridView_))
        {
            std::getline(gstatOutput, line);
            std::istringstream curLine(line);
            if (dim == 1)
                curLine >> trash >> dataValue;
            else if (dim == 2)
                curLine >> trash >> trash >> dataValue;
            else if (dim == 3)
                curLine >> trash >> trash >> trash >> dataValue;
            else
                DUNE_THROW(Dune::InvalidStateException, "Invalid dimension " << dim);

            data_[elementMapper_.index(element)] = dataValue;
        }
        gstatOutput.close();

        // post processing
        using std::pow;
        if (fieldType_ == FieldType::log10)
          std::for_each(data_.begin(), data_.end(), [](Scalar& s){ s = pow(10.0, s); });
    }

    //! \brief Return an entry of the data vector
    Scalar data(const Element& e) const
    {
        return data_[elementMapper_.index(e)];
    }

    //! \brief Return the data vector for analysis or external vtk output
    const DataVector& data() const
    {
        return data_;
    }

    //! \brief Write the data to a vtk file
    void writeVtk(const std::string& vtkName,
                  const std::string& dataName = "data") const
    {
        Dune::VTKWriter<GridView> vtkwriter(gridView_);
        vtkwriter.addCellData(data_, dataName);

        DataVector logPerm;
        if (fieldType_ == FieldType::log10)
        {
            logPerm = data_;
            using std::log10;
            std::for_each(logPerm.begin(), logPerm.end(), [](Scalar& s){ s = log10(s); });
            vtkwriter.addCellData(logPerm, "log10 of " + dataName);
        }
        vtkwriter.write(vtkName, Dune::VTK::OutputType::ascii);
    }

private:
    const GridView gridView_;
    const ElementMapper& elementMapper_;
    DataVector data_;
    FieldType fieldType_;
};

} // end namespace Dumux

#endif
