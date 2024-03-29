#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(altai, 0.0.1, xbz5174@psu.edu, shaomingfu@gmail.com)
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/main.cc])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB

# Check home directories for external libraries
AC_ARG_WITH(boost, AS_HELP_STRING([--with-boost], [home directory for boost]), BOOST_HOME=$withval, BOOST_HOME=)
AC_ARG_WITH(htslib, AS_HELP_STRING([--with-htslib], [home directory for htslib]), HTSLIB_HOME=$withval, HTSLIB_HOME=)
AC_ARG_WITH(zlib, AS_HELP_STRING([--with-zlib], [home directory for zlib]), ZLIB_HOME=$withval, ZLIB_HOME=)
AC_ARG_WITH(lapack, AS_HELP_STRING([--with-lapack], [home directory for lapack]), LAPACK_HOME=$withval, LAPACK_HOME=)

# Assign CXXFLAGS
CXXFLAGS='-Wall -std=c++11'
AS_IF([test "x$BOOST_HOME" != "x"], [AC_SUBST([CXXFLAGS], ["-I$BOOST_HOME $CXXFLAGS"])], [])
AS_IF([test "x$HTSLIB_HOME" != "x"], [AC_SUBST([CXXFLAGS], ["-I$HTSLIB_HOME/include $CXXFLAGS"])], [])

# Check whether static is enabled
AC_ARG_ENABLE([static], [AS_HELP_STRING([--enable-static], [enable static linking])])
AS_IF([test "x$enable_useclp" = "xyes"], [AC_SUBST([CXXFLAGS], ["-DUSECLP $CXXFLAGS"])])

# Check whether enable USECLP
AC_ARG_ENABLE([useclp], [AS_HELP_STRING([--enable-useclp], [use LP to decompose unsplitable vertices])])

# If statis is enabled, directly assign .a files to LIBS
AS_IF([test "x$enable_static" = "xyes"], [AC_SUBST([LIBS], ["-lz $LIBS"])])
AS_IF([test "x$enable_static" = "xyes"], AS_IF([test "x$HTSLIB_HOME" != "x"], [AC_SUBST([LIBS], ["$HTSLIB_HOME/lib/libhts.a $LIBS"])]))
AS_IF([test "x$enable_static" = "xyes"], AS_IF([test "x$enable_useclp" = "xyes"], AS_IF([test "x$CLP_HOME" != "x"], [AC_SUBST([LIBS], ["$CLP_HOME/lib/libClp.a $LIBS"])])))
AS_IF([test "x$enable_static" = "xyes"], AS_IF([test "x$enable_useclp" = "xyes"], AS_IF([test "x$CLP_HOME" != "x"], [AC_SUBST([LIBS], ["$CLP_HOME/lib/libCoinUtils.a $LIBS"])])))

# If static is NOT enabled, assign LDFLAGS and LIBS
AS_IF([test "x$enable_static" != "xyes"], AS_IF([test "x$HTSLIB_HOME" != "x"], [AC_SUBST([LDFLAGS], ["-L$HTSLIB_HOME/lib $LDFLAGS"])]))
AS_IF([test "x$enable_static" != "xyes"], AC_CHECK_LIB([hts], [hts_open], [], [AC_MSG_ERROR([could not find htslib])]))
AS_IF([test "x$enable_static" != "xyes"], AS_IF([test "x$enable_useclp" = "xyes"], AS_IF([test "x$CLP_HOME" != "x"], AC_SUBST([LDFLAGS], ["-L$CLP_HOME/lib $LDFLAGS"]))))
AS_IF([test "x$enable_static" != "xyes"], AS_IF([test "x$enable_useclp" = "xyes"], AC_CHECK_LIB([Clp], [Clp_dual], [], [AC_MSG_ERROR([could not find libClp])])))
AS_IF([test "x$enable_static" != "xyes"], AS_IF([test "x$enable_useclp" = "xyes"], AC_CHECK_LIB([CoinUtils], [_Z9CoinIsnand], [], [AC_MSG_ERROR([could not find libCoinUtils])])))

# Checks for header files.
AC_CHECK_HEADERS([stdint.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT32_T

# Checks for library functions.

AC_CONFIG_FILES([Makefile
				 lib/Makefile
				 lib/util/Makefile
				 lib/gtf/Makefile
				 lib/graph/Makefile
                 src/Makefile])
AC_OUTPUT
