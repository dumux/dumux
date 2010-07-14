// $Id: boundarytypespdelab.hh 3736 2010-06-15 09:52:10Z lauser $
/*****************************************************************************
 *   Copyright (C) 2009-2010 by Bernd Flemisch, Andreas Lauser               *
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
#if HAVE_DUNE_PDELAB && !defined DUMUX_BOUNDARYTYPESPDELAB_HH
#define DUMUX_BOUNDARYTYPESPDELAB_HH

#include<dune/common/fvector.hh>
#include<dune/pdelab/common/function.hh>

template <class TypeTag>
class BoundaryIndexHelperPDELab
{
    enum { numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)) };
public:
    template <int i>
    struct Child {
        struct Type {
            enum { isLeaf = true };
            enum { eqIdx = i };
        };
    };

    template <int i>
    const typename Child<i>::Type &getChild() const
    {
        static typename Child<i>::Type dummy;
        return dummy;
    };

    enum { isLeaf = false };
    enum { CHILDREN = numEq };

    BoundaryIndexHelperPDELab()
    {}
};

#endif
