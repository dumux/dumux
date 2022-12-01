// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup InputOutput
 * \brief  A simple writer class for raster images.
 */
#ifndef DUMUX_RASTER_IMAGE_WRITER_HH
#define DUMUX_RASTER_IMAGE_WRITER_HH

#include <cassert>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <iterator>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dumux/common/stringutilities.hh>
#include <dumux/io/rasterimagedata.hh>

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief A simple reader class for the Netpbm format (https://en.wikipedia.org/wiki/Netpbm_format).
 *        So far, only black and white (*.pbm) and grayscale (*pgm) images are supported.
 */
class NetPBMWriter
{
    template<typename T>
    using Result = Detail::RasterImageData::Result<T>;

    using HeaderData = Detail::RasterImageData::HeaderData;

public:

    template<class ValueType>
    static void write(const std::string& writeFileName,
                      Result<ValueType>& img,
                      const bool useDuneGridOrdering = true)
    {
        // Pass the image to the writer
        writeRasterImageFile_(writeFileName, img, useDuneGridOrdering);
    }

    template<class ValueType>
    static void write(const std::string& writeFileName,
                      const std::size_t& nCols,
                      const std::size_t& nRows,
                      const std::string& magicNumber,
                      const std::string& type,
                      const std::string& encoding,
                      const std::vector<ValueType>& img,
                      const bool useDuneGridOrdering = true)
    {
        // Fill header data and collect img data
        HeaderData headerData;
        headerData.nCols = nCols;
        headerData.nRows = nRows;
        headerData.format.magicNumber = magicNumber;
        Result<ValueType> result(std::move(img), std::move(headerData));

        writeRasterImageFile_(writeFileName, result, useDuneGridOrdering);
    }


    /*!
     * \brief Change the ordering of the pixels according
     *        to Dune's convention, shifting the origin from upper left to lower left.
     *
     * \param result The image's pixel values ordered from top to bottom.
     */
    template<class T>
    static void applyDuneGridOrdering(Result<T>& result)
    {
        auto tmp = result;
        for (std::size_t i = 0; i < result.size(); i += result.header().nCols)
            std::swap_ranges((result.begin() + i), (result.begin() + i + result.header().nCols), (tmp.end() - i - result.header().nCols));
    }

private:

    template <class T>
    static void writeRasterImageFile_(const std::string& writeFileName,
                                      Result<T>& result,
                                      const bool useDuneGridOrdering = true)
    {
            // This will reverse any reordering
        if (useDuneGridOrdering)
            applyDuneGridOrdering(result);

        // Write the corrected image and header information to a file
        std::ofstream outfile(writeFileName, std::ios::trunc);
        outfile << result.header().format.magicNumber << "\n";
        outfile << result.header().nCols << " " << result.header().nRows << "\n";
        if ((result.header().format.magicNumber == "P2") || (result.header().format.magicNumber == "P5"))
        {
            for (int i = 0; i < result.size(); i++)
                outfile << result[i] << "\n";
        }
        else
        {
            for (int i = 0; i < result.size(); i++)
            {
                if (i % result.header().nCols == 0) // wrap once per row
                    outfile << "\n";
                outfile << result[i];
            }
        }
    }

};

} // namespace Dumux

#endif
