# Check for presence of pdfLaTeX
AC_DEFUN([AC_PROG_PDFLATEX],[
  AC_CHECK_PROG(pdflatex, pdflatex, pdflatex)
  if test -z "$pdflatex"; then
    AC_MSG_WARN([pdflatex not found, no handbook will be built.])
  fi
  export pdflatex;

  AC_SUBST(pdflatex)
])