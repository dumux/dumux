AC_DEFUN([CONSTEXPR_CHECK],[
  AC_REQUIRE([AC_PROG_CXX])
  AC_LANG_PUSH([C++])
  AC_MSG_CHECKING([whether constexpr is supported])
  AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM(
      [],
      [[constexpr int peng = 123;]])
    ],
    [HAVE_CONSTEXPR=yes
      AC_MSG_RESULT(yes)],
    [HAVE_CONSTEXPR=no
      AC_MSG_RESULT(no)])
  AC_LANG_POP([C++])
  if test "x$HAVE_CONSTEXPR" = xyes; then
    AC_DEFINE(HAVE_CONSTEXPR, 1, [Define to 1 if constexpr is supported])
  fi
])
