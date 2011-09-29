# Check for the right way to create the deprecation warning
#
# this file is a copy from dune-common which was modified to add a
# version of the macro with an argument. TODO: please remove this file
# as soon as we depend on a DUNE version which provides
# DUNE_DEPRECATED_MSG.

AC_DEFUN([DUMUX_CHECKDEPRECATED],[
	AC_MSG_CHECKING([for deprecated attribute])
        AC_LANG_PUSH([C++])
        AC_TRY_COMPILE([#define DEP __attribute__((deprecated("msg") ))
                    class bar { bar() DEP; };
                    class peng { } DEP;
                    template <class T>
                    class t_bar { t_bar() DEP; };
                    template <class T>
                    class t_peng { t_peng() {}; } DEP;
                    void foo() DEP;
                    void foo() {};],[],
                                  [DUMUX_DEPRECATED="__attribute__((deprecated))"
                                   DUMUX_DEPRECATED_MSG="__attribute__((deprecated(msg) ))"
                     AC_MSG_RESULT(yes)],
                                  [DUMUX_DEPRECATED=""
                                   DUMUX_DEPRECATED_MSG=""
                    AC_MSG_RESULT(no)])
        AC_LANG_POP([C++])

    AC_DEFINE_UNQUOTED(DUMUX_DEPRECATED, $DUMUX_DEPRECATED,
                      [how to create a deprecation warning])
    AC_DEFINE_UNQUOTED(DUMUX_DEPRECATED_MSG(msg), $DUMUX_DEPRECATED_MSG,
                      [how to create a deprecation warning with an additional message])

])
