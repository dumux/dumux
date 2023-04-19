// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief  A data class for raster image information
 */
#ifndef DUMUX_RASTER_IMAGE_DATA_HH
#define DUMUX_RASTER_IMAGE_DATA_HH

#include <string>
#include <vector>
#include <fstream>

namespace Dumux::Detail::RasterImageData {

/*!
    * \brief A struct that holds all information of the image format.
    */
struct Format
{
    std::string magicNumber;
    std::string type;
    std::string encoding;
};

/*!
    * \brief A struct that contains all header data of the image.
    */
struct HeaderData
{
    Format format;
    std::size_t nCols;
    std::size_t nRows;
    std::size_t maxValue = 1;
};

/*!
    * \brief The return type of the reading functions.
    *        Holds the actual pixel values and the header data.
    */
template<class T>
class Result : private std::vector<T>
{
    using Parent = std::vector<T>;
public:
    Result() = delete;

    /*!
        * \brief Construct from data and header by copy
        */
    Result(const std::vector<T>& data, const HeaderData& header)
    : Parent(data)
    , header_(header)
    {}

    /*!
        * \brief Construct from data and header by move
        */
    Result(std::vector<T>&& data, HeaderData&& header)
    : Parent(std::move(data))
    , header_(std::move(header))
    {}

    //! Returns the header data.
    const HeaderData& header() const { return header_; }

    // expose some methods of std::vector
    using Parent::operator[];
    using Parent::begin;
    using Parent::end;
    using Parent::size;

private:
    HeaderData header_;
};

} // end namespace Dumux::Detail::RasterImageData

#endif
