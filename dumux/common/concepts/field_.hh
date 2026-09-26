// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Concepts
 * \brief Concepts for fields that can be bound to an element and evaluated on its faces
 */
#ifndef DUMUX_CONCEPTS_FIELD__HH
#define DUMUX_CONCEPTS_FIELD__HH

#include <concepts>
#include <utility>

namespace Dumux::Concept {

/*!
 * \ingroup Concepts
 * \brief A field that can be bound to a grid element.
 *
 * Binding returns an object rather than mutating the field, so that a bound object
 * outlives neither more nor less than what it owns and a single field can be bound
 * concurrently to several elements.
 */
template<class F, class Element>
concept BindableField = requires(const F& field, const Element& element)
{
    field.bind(element);
};

/*!
 * \ingroup Concepts
 * \brief A bound field that can be evaluated on a face.
 *
 * The face type comes from the discretization the field was bound against, so this is
 * checked where a field and a discretization meet rather than on the field alone.
 */
template<class B, class Face>
concept FaceEvaluatable = requires(const B& bound, const Face& face)
{
    bound(face);
};

/*!
 * \ingroup Concepts
 * \brief A field bindable to an element and evaluatable on that element's faces.
 */
template<class F, class Element, class Face>
concept FaceField = BindableField<F, Element>
    && FaceEvaluatable<decltype(std::declval<const F&>().bind(std::declval<const Element&>())), Face>;

} // end namespace Dumux::Concept

#endif
