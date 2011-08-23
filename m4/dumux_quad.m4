AC_DEFUN([DUMUX_CHECK_QUAD],
[
AC_ARG_WITH([quad],
  AS_HELP_STRING([--with-quad@<:@=yes|no@:>@], [enable quad-precision floating point math (default is no)]),
  [
        if test "$withval" = "no"; then
	   want_quad="no"
        elif test "$withval" = "yes"; then
           want_quad="yes"
        else
           want_quad="no"
        fi
  ],
  [want_quad="no"])

  # store old values
  ac_save_CPPFLAGS="$CPPFLAGS"
  ac_save_LIBS="$LIBS"

  if test "x$want_quad" = "xyes"; then
    QUADMATH_LIB_PATH=""
    QUADMATH_INCLUDE_PATH=""

    QUAD_LIBS="-lquadmath"
    QUAD_CPPFLAGS=""
        
    dnl overwrite ld flags if we have required special directory with
    dnl --with-quad-libdir parameter
    dnl if test "$ac_quad_lib_path" != ""; then
    dnl    QUAD_LDFLAGS="-L$ac_quad_lib_path -lquadmath"
    dnl fi
    LIBS="$LIBS -L/usr/lib64"
    AC_LANG_PUSH(C++)
    AC_SEARCH_LIBS([sqrtq], [quadmath],
    [
      HAVE_QUAD="1"
      QUAD_LIBS="$QUAD_LIBS"
    ],[
      HAVE_QUAD="0"
      AC_MSG_ERROR([libquadmath not found!])
    ])
    AC_LANG_POP([C++])

    AC_SUBST(QUAD_LIBS, $QUAD_LIBS)
    AC_SUBST(QUAD_LDFLAGS, $QUAD_LDFLAGS)
    AC_SUBST(QUAD_CPPFLAGS, $QUAD_CPPFLAGS)
    AC_MSG_RESULT(ok)

    # add to global list
    DUNE_PKG_LIBS="$QUAD_LIBS $DUNE_PKG_LIBS"
    DUNE_PKG_LDFLAGS="$DUNE_PKG_LDFLAGS $QUAD_LDFLAGS"
    DUNE_PKG_CPPFLAGS="$DUNE_PKG_CPPFLAGS $QUAD_CPPFLAGS"

    AC_DEFINE([HAVE_QUAD],1,[Are quad-precision floating point values usable?])
    with_quad="yes"
  else
    with_quad="no"
  fi

  # tell automake
  AM_CONDITIONAL(QUAD, test "x$HAVE_QUAD" = "x1")

  # restore variables
  LDFLAGS="$ac_save_LDFLAGS"
  CPPFLAGS="$ac_save_CPPFLAGS"
  LIBS="$ac_save_LIBS"

  DUNE_ADD_SUMMARY_ENTRY([quadruple precision math],[$with_quad])
])
