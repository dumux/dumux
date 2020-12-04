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
 * \brief  A simple reader class for raster images.
 */
#ifndef DUMUX_RASTER_IMAGE_READER_HH
#define DUMUX_RASTER_IMAGE_READER_HH

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

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief A simple reader class for the Netpbm format (https://en.wikipedia.org/wiki/Netpbm_format).
 *        So far, only black and white (*.pbm) and grayscale (*pgm) images are supported.
 */
class NetPBMReader
{
public:
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
         * \brief Contruct from data and header by copy
         */
        Result(const std::vector<T>& data, const HeaderData& header)
        : Parent(data)
        , header_(header)
        {}

        /*!
         * \brief Contruct from data and header by move
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

    /*!
     * \brief A helper function to retrieve the format from tokens of the file's first line.
     *
     * \param firstLineTokes The tokens extracted from the first line of the file
     */
    static Format getFormat(const std::vector<std::string_view>& firstLineTokes)
    {
        static const auto format = []{
            std::map<std::string, Format> format;
            format["P1"] = Format{"P1", "Portable BitMap", "ASCII"};
            format["P2"] = Format{"P2", "Portable GrayMap", "ASCII"};
            format["P3"] = Format{"P3", "Portable PixMap", "ASCII"};
            format["P4"] = Format{"P4", "Portable BitMap", "binary"};
            format["P5"] = Format{"P5", "Portable GrayMap", "binary"};
            format["P6"] = Format{"P6", "Portable PixMap", "binary"};
            return format;
        }();

        const std::string& magicNumber = std::string(firstLineTokes[0]);

        if (!format.count(magicNumber))
            DUNE_THROW(Dune::IOError, magicNumber << " is not a valid magic number for the Netpbm format");

        return format.at(magicNumber);
    }

    /*!
     * \brief Reads a *pbm (black and white) in ASCII or binary encoding.
     *        Returns a struct that contains both the pixel values and the header data.
     *
     * \param fileName The file name (*.pbm).
     * \param useDuneGridOrdering If set true, the ordering of the pixels will be changed according
     *                            to Dune's convention, shifting the origin from upper left to lower left.
     */
    static Result<bool> readPBM(const std::string& fileName, const bool useDuneGridOrdering = true)
    {
        std::ifstream infile(fileName, std::ios::binary);

        if (!infile)
            DUNE_THROW(Dune::IOError, "Raster data file not found or corrupt");

        HeaderData headerData = readHeader(infile);
        std::vector<bool> values;

        if (headerData.format.magicNumber == "P1")
            values = readPBMDataASCII_(infile, headerData);
        else if (headerData.format.magicNumber == "P4")
            values = readPBMDataBinary_(infile, headerData);
        else
            DUNE_THROW(Dune::IOError, headerData.format.magicNumber << " not supported. Use readPBM for P1 and P4 or readPGM for P2 and P5");

        Result<bool> result(std::move(values), std::move(headerData));
        printInfo(result);

        if (useDuneGridOrdering)
            applyDuneGridOrdering(result);

        return result;
    }

    /*!
     * \brief Reads a *.pgm (grayscale) in ASCII or binary encoding.
     *        Returns a struct that contains both the pixel values and the header data.
     *
     * \tparam ValueType The value type representing the pixel data. By default, std::uint8_t (0-255) is used.
     *                   Since this type is often defined as unsigned char, some conversion issues during I/O may occur.
     *                   Hence the type may be changed, e.g., to std::size_t.
     *
     * \param fileName The file name (*.pbm).
     * \param useDuneGridOrdering If set true, the ordering of the pixels will be changed according
     *                            to Dune's convention, shifting the origin from upper left to lower left.
     */
    template<class ValueType = std::uint8_t>
    static Result<ValueType> readPGM(const std::string& fileName, const bool useDuneGridOrdering = true)
    {
        std::ifstream infile(fileName, std::ios::binary);

        if (!infile)
            DUNE_THROW(Dune::IOError, "Raster data file not found or corrupt");

        HeaderData headerData = readHeader(infile);
        std::vector<ValueType> values;

        if (headerData.format.magicNumber == "P2")
            values = NetPBMReader::template readPGMDataASCII_<ValueType>(infile, headerData);
        else if (headerData.format.magicNumber == "P5")
            values = NetPBMReader::template readPGMDataBinary_<ValueType>(infile, headerData);
        else
            DUNE_THROW(Dune::IOError, headerData.format.magicNumber << " not supported. Use readPBM for P1 and P4 or readPGM for P2 and P5");

        Result<ValueType> result(std::move(values), std::move(headerData));
        printInfo(result);

        if (useDuneGridOrdering)
            applyDuneGridOrdering(result);

        return result;
    }

    /*!
     * \brief Returns the header data of the image file.
     *
     * \param infile The input stream used to process the image file.
     */
    static HeaderData readHeader(std::ifstream& infile)
    {
        HeaderData headerData;
        std::string inputLine;
        std::size_t lineNumber = 0;

        // First line : get format.
        std::getline(infile, inputLine);
        ++lineNumber;

        const auto firstLineTokens = tokenize(inputLine, " ");
        headerData.format = getFormat(firstLineTokens);
        const auto magicNumber = headerData.format.magicNumber;

        // dimensions could be given right after magic number (in same line)
        if (firstLineTokens.size() > 2)
        {
            if (isBlackAndWhite_(magicNumber) && firstLineTokens.size() != 3)
                DUNE_THROW(Dune::IOError, "Could not read first line for B/W image");

            headerData.nCols = std::stoi(std::string(firstLineTokens[1]));
            headerData.nRows = std::stoi(std::string(firstLineTokens[2]));

            if (isGrayScale_(magicNumber))
            {
                if (firstLineTokens.size() == 4)
                    headerData.maxValue = std::stoi(std::string(firstLineTokens[3]));
                if (firstLineTokens.size() > 4)
                    DUNE_THROW(Dune::IOError, "Could not read first line for grayscale image");
            }
        }
        else
        {
            // Read dimensions and maximum value (for non-b/w images).
            while (!infile.eof())
            {
                std::getline(infile, inputLine);
                ++lineNumber;

                // Skip comments.
                if (isComment_(inputLine))
                    continue;

                const auto tokens = tokenize(inputLine, " ");

                // The first line after the comments contains the dimensions.
                if (tokens.size() != 2)
                    DUNE_THROW(Dune::IOError, "Expecting " << [](auto size){ return size < 2 ? "both" : "only"; }(tokens.size()) << " dimensions (2 numbers) in line " << lineNumber);

                headerData.nCols = std::stoi(std::string(tokens[0]));
                headerData.nRows = std::stoi(std::string(tokens[1]));

                // Grayscale images additionaly contain a maxium value in the header.
                if (isGrayScale_(magicNumber))
                {
                    std::getline(infile, inputLine);
                    ++lineNumber;

                    const auto token = tokenize(inputLine, " ");
                    if (token.size() != 1)
                        DUNE_THROW(Dune::IOError, "Expecting" << [](auto size){ return size == 0 ? "" : " only"; }(token.size()) << " intensity (one number) in line " << lineNumber);

                    headerData.maxValue = std::stoi(std::string(token[0]));
                }
                break;
            }
        }

        return headerData;
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

    /*!
     * \brief Print the data contained in the header.
     *
     * \param result The object storing both the image's pixel values and the header data.
     */
    template<class T>
    static void printInfo(const Result<T>& result)
    {
        const Format& format = result.header().format;
        std::cout << "Reading " << format.type << " file (" << format.encoding << ")" << std::endl;
        std::cout << "Dimensions : " << result.header().nCols << " " << result.header().nRows << std::endl;
        std::cout << "Maximum value : " << result.header().maxValue << std::endl;
    }

    /*!
     * \brief Fill a pre-defined 2D image object, e.g. std::vector<std::vector<bool>>,
     *        with the pixel values stored in a 1D container.
     *
     * \param image The 2D image to be filled with values. Needs to have to correct dimensions (nCols x nRows).
     * \param result The image's pixel values stored in a 1D container.
     */
    template<class Image, class T>
    static void fillImage(Image& image, const Result<T>& result)
    {
        const auto nCols = result.header().nCols;
        const auto nRows = result.header().nRows;
        using RowType = std::decay_t<decltype(image[0])>;
        image.resize(nRows, RowType(nCols));

        std::size_t rowIdx = 0;
        std::size_t colIdx = 0;
        for (const auto& val : result)
        {
            image[rowIdx][colIdx] = val;

            // start a new line after nCols entries
            if (++colIdx == nCols)
            {
                colIdx = 0;
                ++rowIdx;
            }
        }
    }

    /*!
     * \brief Flattens a 2D image object to a 1D container.
     *
     * \param image The 2D image to be flattened.
     */
    template<class Image>
    static auto flattenImageToVector(const Image& image)
    {
        using OutputValueType = std::decay_t<decltype(image[0][0])>;

        std::vector<OutputValueType> data;
        data.reserve(image.size()*image[0].size());
        for (const auto& row: image)
            for (const auto& col : row)
                data.push_back(col);

        return data;
    }

private:

    static bool isBlackAndWhite_(const std::string& magicNumber)
    {
        return magicNumber == "P1" || magicNumber == "P4";
    }

    static bool isGrayScale_(const std::string& magicNumber)
    {
        return magicNumber == "P2" || magicNumber == "P5";
    }

    /*!
     * \brief Reads the data block of a *.pbm (black and white) file in ASCII encoding.
     *        Returns a vector that contains the pixel values.
     * \note readHeader(infile) must be called prior to calling this function. This assumes all values to be either 0 or 1.
     *
     * \param infile The input stream used to process the image file.
     * \param headerData The data contained in the file's header.
     */
    static std::vector<bool> readPBMDataASCII_(std::ifstream& infile,
                                               const HeaderData& headerData)
    {
        std::string inputLine;
        std::vector<bool> data;
        data.reserve(numPixel_(headerData));

        while (!infile.eof())
        {
            std::getline(infile, inputLine);
            if (!isComment_(inputLine))
            {
                inputLine.erase(std::remove_if(inputLine.begin(), inputLine.end(), [](unsigned char c){ return std::isspace(c); }), inputLine.end());
                if (!inputLine.empty())
                {
                    for (const auto& value : inputLine)
                    {
                        assert(value == '0' || value == '1');
                        data.push_back(value - '0');  // convert char to int
                    }
                }
            }
        }

        return data;
    }

    /*!
     * \brief Reads the data block of a *.pbm (black and white) file in binary encoding.
     *        Returns a vector that contains the pixel values.
     * \note readHeader(infile) must be called prior to calling this function.
     *
     * \param infile The input stream used to process the image file.
     * \param headerData The data contained in the file's header.
     */
    static std::vector<bool> readPBMDataBinary_(std::ifstream& infile,
                                                const HeaderData& headerData)
    {
        std::vector<bool> data(numPixel_(headerData));

        // Skip potentially remaining comments in header section
        // before reading binary content. We detect a comment by
        // reading a line with std::getline and checking the resulting string.
        // We continue reading new lines until no more comments are found. Then we
        // need to set infile's current position to one line before the actual binary
        // content, otherwise the following steps will fail.
        std::string inputLine;
        while (!infile.eof())
        {
            // store position before calling std::getline
            const auto lastPos = infile.tellg();
            std::getline(infile, inputLine);

            // stop the loop if no more comment is found and go back one line
            if (!isComment_(inputLine))
            {
                infile.seekg(lastPos);
                break;
            }
        }

        // read actual binary content
        std::size_t nBytes = 0;
        std::size_t bitIndex = 0;
        using Bit = std::uint8_t;
        Bit b = 0;
        for (std::size_t j = 0; j < headerData.nRows; j++)
        {
            for (std::size_t i = 0; i < headerData.nCols; i++)
            {
                if (i%8 == 0)
                {
                    char tmp;
                    infile.read(&tmp, 1);
                    b = static_cast<Bit>(tmp);
                    if (infile.eof())
                        DUNE_THROW(Dune::IOError, "Failed reading byte " << nBytes);

                    ++nBytes;
                }

                const Bit k = 7 - (i % 8);
                data[bitIndex++] = static_cast<bool>((b >> k) & 1);
            }
        }

        return data;
    }

    /*!
     * \brief Reads the data block of a *.pgm (grayscale) file in ASCII encoding.
     *        Returns a vector that contains the pixel values.
     * \note readHeader(infile) must be called prior to calling this function.
     *
     * \tparam ValueType The value type representing the pixel data. By default, std::uint8_t (0-255) is used.
     *                   Since this type is often defined as unsigned char, some conversion issues during I/O may occur.
     *                   Hence the type may be changed, e.g., to std::size_t.
     *
     * \param infile The input stream used to process the image file.
     * \param headerData The data contained in the file's header.
     */
    template<class ValueType = std::uint8_t>
    static std::vector<ValueType> readPGMDataASCII_(std::ifstream& infile,
                                                    const HeaderData& headerData)
    {
        std::string inputLine;

        std::vector<ValueType> data;
        data.reserve(numPixel_(headerData));

        while (!infile.eof())
        {
            std::getline(infile, inputLine);
            if (inputLine.empty())
                continue;

            // if the line contains multiple comma-separated values, store them individually in a vector
            if (inputLine.find(" ") != std::string::npos)
            {
                std::istringstream iss(inputLine);
                std::vector<std::string> tokens;
                std::copy(std::istream_iterator<std::string>(iss),
                          std::istream_iterator<std::string>(),
                          std::back_inserter(tokens));

                for (const auto& t : tokens)
                    data.push_back(std::stoi(t)); // convert string to integer type
            }
            else
                data.push_back(std::stoi(inputLine));
        }

        return data;
    }

    /*!
     * \brief Reads the data block of a *.pgm (grayscale) file in binary encoding.
     *        Returns a vector that contains the pixel values.
     * \note readHeader(infile) must be called prior to calling this function.
     *
     * \tparam ValueType The value type representing the pixel data. By default, std::uint8_t (0-255) is used.
     *                   Since this type is often defined as unsigned char, some conversion issues during I/O may occur.
     *                   Hence the type may be changed, e.g., to std::size_t.
     *
     * \param infile The input stream used to process the image file.
     * \param headerData The data contained in the file's header.
     */
    template<class ValueType = std::uint8_t>
    static std::vector<ValueType> readPGMDataBinary_(std::ifstream& infile,
                                                     const HeaderData& headerData)
    {
        std::string inputLine;

        // check the size of the binary part of the file
        const auto curPos = infile.tellg();
        infile.seekg(0, std::ios::end);
        const auto endPos = infile.tellg();
        const auto size = endPos - curPos;
        if (size != numPixel_(headerData))
            DUNE_THROW(Dune::IOError, "Binary file size does not match with raster image size");

        // reset to the current position
        infile.seekg(curPos, std::ios::beg);

        // extract the binary data
        std::vector<std::uint8_t> data(size);
        infile.read(reinterpret_cast<char*>(&data[0]), size);

        // convert std::uint8_t to ValueType
        return std::vector<ValueType>(data.begin(), data.end());
    }

    /*!
     * \brief Returns the image's number of pixels.
     *
     * \param headerData The data contained in the file's header.
     */
    static std::size_t numPixel_(const HeaderData& headerData)
    {
        return headerData.nRows*headerData.nCols;
    }

    /*!
     * \brief Returns true if a given line is a comment (starting with #)
     */
    static bool isComment_(const std::string_view line)
    {
        return line[0] == '#';
    }
};

} // namespace Dumux

#endif
