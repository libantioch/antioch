 dnl                                               -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(antioch, 0.0.0, pbauman@ices.utexas.edu)
AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_HEADER(antioch_config.h.tmp)
dnl Prefixes everything in the temporary file with FINS_
AX_PREFIX_CONFIG_H(antioch_config.h,ANTIOCH,antioch_config.h.tmp)

AC_CONFIG_AUX_DIR([build-aux])
AC_CANONICAL_TARGET()
AM_INIT_AUTOMAKE([dist-bzip2 tar-ustar 1.11 no-define color-tests parallel-tests])
AM_SILENT_RULES(yes)  # use silent rules where available - automake 1.11

LT_INIT

dnl Release versioning
AX_SPLIT_VERSION
GENERIC_MAJOR_VERSION=$AX_MAJOR_VERSION
GENERIC_MINOR_VERSION=$AX_MINOR_VERSION
GENERIC_MICRO_VERSION=$AX_POINT_VERSION
GENERIC_LIBRARY_NAME=antioch

dnl shared library versioning

PACKAGE=$GENERIC_LIBRARY_NAME
GENERIC_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION.$GENERIC_MICRO_VERSION
GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION

dnl AC_SUBST(GENERIC_LIBRARY_VERSION)
dnl AC_SUBST(GENERIC_LIBRARY_NAME)

AC_SUBST(GENERIC_VERSION)
AC_SUBST(GENERIC_RELEASE)
AC_SUBST(GENERIC_MAJOR_VERSION)
AC_SUBST(GENERIC_MINOR_VERSION)
AC_SUBST(GENERIC_MICRO_VERSION)

AC_DEFINE_UNQUOTED([MAJOR_VERSION],$GENERIC_MAJOR_VERSION,[Major version])
AC_DEFINE_UNQUOTED([MINOR_VERSION],$GENERIC_MINOR_VERSION,[Minor version])
AC_DEFINE_UNQUOTED([MICRO_VERSION],$GENERIC_MICRO_VERSION,[Micro version])



dnl------------------------------
dnl Checks for compilers
dnl------------------------------

AC_PROG_CC
AC_PROG_CXX

dnl Added for AX_CODE_COVERAGE macro to work correctly, even though there's no Fortran here.
AC_PROG_FC

dnl--------------------------
dnl Checks for code coverage
dnl--------------------------
AX_CODE_COVERAGE

dnl---------------------------------
dnl Query configuration environment
dnl Must be done *before* Doxygen to
dnl get some happy variables defined.
dnl---------------------------------
AX_SUMMARIZE_ENV

dnl------------------------------------------
dnl            Doxygen support
dnl------------------------------------------

dnl DX_HTML_FEATURE(ON)
dnl DX_CHM_FEATURE(OFF)
dnl DX_CHI_FEATURE(OFF)
dnl DX_MAN_FEATURE(OFF)
dnl DX_RTF_FEATURE(OFF)
dnl DX_XML_FEATURE(OFF)
dnl DX_PDF_FEATURE(ON)
dnl DX_PS_FEATURE(OFF)

dnl DX_INIT_DOXYGEN(ANTIOCH, doxygen/antioch.dox, docs/doxygen)

dnl AC_CONFIG_FILES(doxygen/txt_common/about.page)

dnl-----------------------------------------------
dnl Generate files for unit and regression testing
dnl-----------------------------------------------
dnl AC_CONFIG_FILES(test/test_axi_ns_poiseuille_flow.sh, [chmod +x test/test_axi_ns_poiseuille_flow.sh])

dnl-----------------------------------------------
dnl Generate header files
dnl-----------------------------------------------
AC_CONFIG_FILES(src/utilities/include/antioch/antioch_version.h)
            
dnl-----------------------------------------------
dnl Generate Makefiles
dnl-----------------------------------------------
AC_CONFIG_FILES([Makefile src/Makefile test/Makefile]) dnl  doxygen/Makefile examples/Makefile])

# Must still call AC_OUTPUT() after generating all the files
AC_OUTPUT()

dnl------------------------------
dnl Final summary
dnl------------------------------
AX_SUMMARIZE_CONFIG
