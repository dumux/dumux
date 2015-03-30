dnl This macro introduces a configure flag --disable-dumux-deprecated-autotools
dnl to suppress the deprecration warning that Autotools will be removed in 2.8

AC_DEFUN([DUMUX_DEPRECATED_AUTOTOOLS],[
  AC_ARG_ENABLE(dumux-deprecated-autotools,
    AS_HELP_STRING([--disable-dumux-deprecated-autotools],[Suppresses DuMuX' warning about deprecated Autotools support.]))

  AS_IF([test "x$enable_dumux_deprecated_autotools" != "xno"],
    AC_DEFINE(ENABLE_DUMUX_DEPRECATED_AUTOTOOLS, 1, [Warn about deprecated Autotools support]))

  AH_BOTTOM([
#if ENABLE_DUMUX_DEPRECATED_AUTOTOOLS
#warning Support for Autotools in DuMuX is deprecated and will be removed after 2.7.
#endif
])
])
