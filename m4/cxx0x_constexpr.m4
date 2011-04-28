AC_DEFUN([CONSTEXPR_CHECK],[
  AC_REQUIRE([AC_PROG_CXX])
  AC_REQUIRE([GXX0X])
  AC_LANG_PUSH([C++])
  AC_MSG_CHECKING([whether constexpr is supported])
  AC_TRY_COMPILE([],[constexpr int peng = 123;], [
    HAVE_CONSTEXPR=yes
    AC_MSG_RESULT(yes)], [
    HAVE_CONSTEXPR=no
    AC_MSG_RESULT(no)])
  if test "x$HAVE_CONSTEXPR" = xyes; then
    AC_DEFINE(HAVE_CONSTEXPR, 1, [Define to 1 if constexpr is supported])
  fi
  AC_LANG_POP
])
