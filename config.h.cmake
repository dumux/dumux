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

#cmakedefine HAVE_MPI 1

#cmakedefine HAVE_UG 1
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
#define HAVE_MEM_USAGE_T_EXPANSIONS 1
#endif

/* tr1/array. */
//#cmakedefine HAVE_TR1_ARRAY 1
#cmakedefine HAVE_NULLPTR 1
#cmakedefine HAVE_STATIC_ASSERT 1
#cmakedefine HAVE_VARIADIC_TEMPLATES 1
#cmakedefine HAVE_VARIADIC_CONSTRUCTOR_SFINAE 1
#cmakedefine HAVE_RVALUE_REFERENCES 1
#cmakedefine HAVE_MALLOC_H 1
#cmakedefine HAVE_VALGRIND 1

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
