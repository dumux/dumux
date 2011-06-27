/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Simplifies handling of buffers to be used in conjunction with MPI
 */
#ifndef DUMUX_MPI_BUFFER_HH
#define DUMUX_MPI_BUFFER_HH

#if HAVE_MPI
#include <mpi.h>
#endif

#include <type_traits>

namespace Dumux
{
/*!
 * \brief Simplifies handling of buffers to be used in conjunction with MPI
 */
template <class DataType>
class MpiBuffer
{
public:
    MpiBuffer(int size)
    {
        data_ = new DataType[size];
        dataSize_ = size;

#if HAVE_MPI
        mpiDataSize_ = size;

        // set the MPI data type
        if (std::is_same<DataType, char>::value)
            mpiDataType_ = MPI_CHAR;
        else if (std::is_same<DataType, unsigned char>::value)
            mpiDataType_ = MPI_UNSIGNED_CHAR;
        else if (std::is_same<DataType, short>::value  || std::is_same<DataType, unsigned short>::value)
            mpiDataType_ = MPI_SHORT;
        else if (std::is_same<DataType, int>::value || std::is_same<DataType, unsigned>::value)
            mpiDataType_ = MPI_INT;
        else if (std::is_same<DataType, long>::value || std::is_same<DataType, unsigned long>::value)
            mpiDataType_ = MPI_LONG;
        else if (std::is_same<DataType, float>::value)
            mpiDataType_ = MPI_FLOAT;
        else if (std::is_same<DataType, double>::value)
            mpiDataType_ = MPI_DOUBLE;
        else if (std::is_same<DataType, long double>::value)
            mpiDataType_ = MPI_LONG_DOUBLE;
        else {
            mpiDataType_ = MPI_BYTE;
            mpiDataSize_ *= sizeof(DataType);
        }
#endif // HAVE_MPI
    };

    ~MpiBuffer()
    {
        delete[] data_;
    }

    /*!
     * \brief Send the buffer asyncronously to a peer process.
     */
    void send(int peerRank, bool setNoAccess=true)
    {
#if HAVE_MPI
        MPI_Isend(data_,
                  mpiDataSize_,
                  mpiDataType_,
                  peerRank,
                  0, // tag
                  MPI_COMM_WORLD,
                  &mpiRequest_);
#endif
    };

    /*!
     * \brief Wait until the buffer was send to the peer completely.
     */
    void wait()
    {
#if HAVE_MPI
        MPI_Wait(&mpiRequest_, &mpiStatus_);
#endif // HAVE_MPI
    };

    /*!
     * \brief Receive the buffer syncronously from a peer rank
     */
    void receive(int peerRank)
    {
#if HAVE_MPI
        MPI_Recv(data_,
                 mpiDataSize_,
                 mpiDataType_,
                 peerRank,
                 0, // tag
                 MPI_COMM_WORLD,
                 &mpiStatus_);
        assert(! mpiStatus_.MPI_ERROR);
#endif // HAVE_MPI
    };

#if HAVE_MPI
    /*!
     * \brief Returns the current MPI_Request object.
     *
     * This object is only well defined after the send() method.
     */
    MPI_Request &request()
    { return mpiRequest_; }
    /*!
     * \brief Returns the current MPI_Request object.
     *
     * This object is only well defined after the send() method.
     */
    const MPI_Request &request() const
    { return mpiRequest_; }

    /*!
     * \brief Returns the current MPI_Status object.
     *
     * This object is only well defined after the receive() and wait() methods.
     */
    MPI_Status &status()
    { return mpiStatus_; }
    /*!
     * \brief Returns the current MPI_Status object.
     *
     * This object is only well defined after the receive() and wait() methods.
     */
    const MPI_Status &status() const
    { return mpiStatus_; }
#endif // HAVE_MPI


    /*!
     * \brief Returns the number of data objects in the buffer
     */
    int size() const
    { return dataSize_; }

    /*!
     * \brief Provide access to the buffer data.
     */
    DataType &operator[](int i)
    {
        assert(0 <= i && i < dataSize_);
        return data_[i];
    }

    /*!
     * \brief Provide access to the buffer data.
     */
    const DataType &operator[](int i) const
    {
        assert(0 <= i && i < dataSize_);
        return data_[i];
    }

private:
    DataType *data_;
    int dataSize_;
#if HAVE_MPI
    int mpiDataSize_;
    MPI_Datatype mpiDataType_;
    MPI_Request mpiRequest_;
    MPI_Status mpiStatus_;
#endif // HAVE_MPI
};

}

#endif
