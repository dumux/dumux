AC_DEFUN([CHECK_PATCHED_PDELAB],[
  AC_REQUIRE([AC_PROG_CXX])
  AC_LANG_PUSH([C++])
  ac_save_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $DUNE_COMMON_CPPFLAGS $DUNE_ISTL_CPPFLAGS $DUNE_PDELAB_CPPFLAGS $DUNE_PKG_CPPFLAGS"
  AC_MSG_CHECKING([whether dune-pdelab is patched to be usable by DuMuX])
  AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM(
      [#include <dune/pdelab/backend/istlvectorbackend.hh>],
      [[
        Dune::PDELab::ISTLBlockVectorContainer
          <std::vector<double>, double, 3> blockVectorContainer;
        return blockVectorContainer.size();
      ]])
    ],
    [DUNE_PDELAB_IS_PATCHED_FOR_DUMUX=yes
      AC_MSG_RESULT(yes)],
    [DUNE_PDELAB_IS_PATCHED_FOR_DUMUX=no
      AC_MSG_RESULT(no)])
  CPPFLAGS="$ac_save_CPPFLAGS"
  AC_LANG_POP([C++])

  # add variable to config.h
  if test "x$DUNE_PDELAB_IS_PATCHED_FOR_DUMUX" = xyes; then
    AC_DEFINE(DUNE_PDELAB_IS_PATCHED_FOR_DUMUX, 1,
      [Define to 1 if dune-pdelab is patched to be usable by DuMuX])
  fi
])
