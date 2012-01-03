dnl @synopsis AC_PROG_DVIPDF
dnl
dnl This macro test if dvipdf is installed. If dvipdf is installed, it
dnl set $dvipdf to the right value
dnl
dnl @category InstalledPackages
dnl @category LaTeX
dnl @author Mathieu Boretti <boretti@bss-network.com>
dnl @version 2005-01-21
dnl @license GPLWithACException

AC_DEFUN([AC_PROG_DVIPDF],[
AC_CHECK_PROGS(dvipdf,dvipdf,no)
export dvipdf;
if test $dvipdf = "no" ;
then
	AC_MSG_ERROR([Unable to find a dvipdf application]);
fi;
AC_SUBST(dvipdf)
])
