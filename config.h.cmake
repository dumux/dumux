#ifdef CONFIG_H
#  error "config.h included more than once!"
#endif
#define CONFIG_H

#define DUNE_MINIMAL_DEBUG_LEVEL 4
#cmakedefine HAVE_BOOST 1
#cmakedefine HAVE_DUNE 1
#cmakedefine HAVE_DUNE_GRID 1
#cmakedefine HAVE_DUNE_DISC 1
#cmakedefine HAVE_DUNE_FEM 1
#cmakedefine HAVE_DUNE_ISTL 1
#cmakedefine HAVE_DUNE_LOCALFUNCTIONS 1
#cmakedefine HAVE_DUNE_PDELAB 1

/* If this is set, the member 'size' of FieldVector is a method rather than an
   enum */
#define DUNE_COMMON_FIELDVECTOR_SIZE_IS_METHOD 1

/* Define to the version of dune-common */
#cmakedefine DUNE_COMMON_VERSION "${DUNE_COMMON_VERSION}"

/* Define to the major version of dune-common */
#cmakedefine DUNE_COMMON_VERSION_MAJOR ${DUNE_COMMON_VERSION_MAJOR}

/* Define to the minor version of dune-common */
#define DUNE_COMMON_VERSION_MINOR ${DUNE_COMMON_VERSION_MINOR}

/* Define to the revision of dune-common */
#define DUNE_COMMON_VERSION_REVISION ${DUNE_COMMON_VERSION_REVISION}

#cmakedefine HAVE_MPI 1

#cmakedefine HAVE_UG ENABLE_UG
#if HAVE_MPI && HAVE_UG
// use parallel UG if both UG and MPI are available
#   define ModelP
#endif

#ifdef ENABLE_ALUGRID
#cmakedefine HAVE_ALUGRID 1
#endif

#ifdef ENABLE_METIS
#cmakedefine HAVE_METIS 1
#endif

#ifdef ENABLE_ALBERTA
#cmakedefine HAVE_ALBERTA 1
#endif

#cmakedefine PROJECT_NAME             "${PROJECT_NAME}"
#cmakedefine PROJECT_VERSION          "${PROJECT_VERSION}"
#cmakedefine PROJECT_MAINTAINER       "${PROJECT_MAINTAINER}"
#cmakedefine PROJECT_MAINTAINER_EMAIL "${PROJECT_MAINTAINER_EMAIL}"

#cmakedefine HAVE_SUPERLU ENABLE_SUPERLU
#ifdef HAVE_SUPERLU
#define SUPERLU_POST_2005_VERSION
#cmakedefine SUPERLU_MIN_VERSION_4_3
#endif

/* tr1/array. */
/*#cmakedefine HAVE_TR1_ARRAY 1*/

/* Define to 1 if the <array> C++0x is available and support array::fill */
#cmakedefine HAVE_ARRAY 1

/* Define to 1 if you have the <memory> header file. */
#cmakedefine HAVE_MEMORY 1

/* The namespace in which SHARED_PTR can be found */
#cmakedefine SHARED_PTR_NAMESPACE ${SHARED_PTR_NAMESPACE}

/* The header in which SHARED_PTR can be found */
#cmakedefine SHARED_PTR_HEADER ${SHARED_PTR_HEADER}

/* Define to 1 if SHARED_PTR_NAMESPACE::make_shared is usable */
#cmakedefine HAVE_MAKE_SHARED 1

/* Define to 1 if you have <boost/make_shared.hpp> */
#cmakedefine HAVE_BOOST_MAKE_SHARED_HPP 1

#cmakedefine HAVE_NULLPTR 1
#cmakedefine HAVE_STATIC_ASSERT 1
#cmakedefine HAVE_VARIADIC_TEMPLATES 1
#cmakedefine HAVE_VARIADIC_CONSTRUCTOR_SFINAE 1
#cmakedefine HAVE_RVALUE_REFERENCES 1
#cmakedefine HAVE_MALLOC_H 1
#cmakedefine HAVE_VALGRIND 1

#include <dune/common/deprecated.hh>

#if DUNE_COMMON_VERSION_MAJOR >= 2 && DUNE_COMMON_VERSION_MINOR >= 2
#include <dune/common/unused.hh>
#endif

#cmakedefine HAVE_ATTRIBUTE_DEPRECATED 1
#if HAVE_ATTRIBUTE_DEPRECATED
#  define DUMUX_DEPRECATED __attribute__((deprecated))
#else
#  define DUMUX_DEPRECATED
#endif

#cmakedefine HAVE_ATTRIBUTE_DEPRECATED_MSG 1
#if HAVE_ATTRIBUTE_DEPRECATED_MSG
#  define DUMUX_DEPRECATED_MSG(msg) __attribute__((deprecated(msg)))
#else
#  define DUMUX_DEPRECATED_MSG DUMUX_DEPRECATED
#endif

#cmakedefine HAVE_ATTRIBUTE_ALWAYS_INLINE 1
#if HAVE_ATTRIBUTE_ALWAYS_INLINE
#  define DUMUX_ALWAYS_INLINE __attribute__((always_inline))
#else
#  define DUMUX_ALWAYS_INLINE
#endif

#cmakedefine HAVE_ATTRIBUTE_UNUSED 1
#if HAVE_ATTRIBUTE_UNUSED
#  define DUMUX_UNUSED __attribute__((unused))
#else
#  define DUMUX_UNUSED
#endif
