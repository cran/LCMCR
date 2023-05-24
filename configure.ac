AC_INIT([LCMCR],[0.4.13],[dmanriqu@indiana.edu])
AC_CONFIG_MACRO_DIRS([.])
AC_CONFIG_MACRO_DIRS([tools])
AX_PATH_GSL([2.5],,
    AC_MSG_ERROR([could not find required version of GSL]))
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)

AC_C_INLINE
if test "$ac_cv_c_inline" != no ; then
    INL=-DHAVE_INLINE
    AC_SUBST(INL)
fi

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
