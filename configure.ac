 dnl                                               -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([Antioch], [0.2.0], [antioch-users@googlegroups.com], [antioch])
PACKAGE_DESCRIPTION="A New, Templated Implementation of Chemistry for Hydrodynamics"
AC_SUBST([PACKAGE_DESCRIPTION])
PACKAGE_URL="https://github.com/libantioch/antioch/"
AC_SUBST([PACKAGE_URL])
AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_HEADER(antioch_config.h.tmp)
dnl Prefixes everything in the temporary file with FINS_
AX_PREFIX_CONFIG_H(antioch_config.h,ANTIOCH,antioch_config.h.tmp)

AC_CONFIG_AUX_DIR([build-aux])
AC_CANONICAL_TARGET()
AM_INIT_AUTOMAKE([dist-bzip2 tar-ustar filename-length-max=99 1.11 no-define color-tests parallel-tests subdir-objects])
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

dnl -Wall warnings, -Wall the time.
AX_CXXFLAGS_WARN_ALL

dnl Allow users to make vectorized tests larger
AC_ARG_WITH([n_tuples],
	    AC_HELP_STRING([--with-n_tuples=<1...1000000000>],
			   [number of tuples used in each vector test]),
	    [n_tuples="$withval"],
	    [n_tuples=2])
AC_DEFINE_UNQUOTED(N_TUPLES, $n_tuples, [number of tuples used in each vector test])

dnl Allow users to make vectorized tests longer
AC_ARG_WITH([n_runs],
	    AC_HELP_STRING([--with-n_runs=<1...1000000000>],
			   [number of runs done in each vector test]),
	    [n_runs="$withval"],
	    [n_runs=1])
AC_DEFINE_UNQUOTED(N_RUNS, $n_runs, [number of runs done in each vector test])


dnl--------------------------
dnl Checks for third-party libraries
dnl--------------------------
antioch_optional_INCLUDES=""
antioch_optional_test_INCLUDES=""
antioch_optional_test_LDFLAGS=""
antioch_optional_test_LIBS=""

dnl Enable VexCL if found
AX_PATH_VEXCL(0.0.0, no)
if (test x$HAVE_VEXCL = x1); then
  antioch_optional_test_INCLUDES="$VEXCL_CPPFLAGS $antioch_optional_test_INCLUDES"
  antioch_optional_test_LDFLAGS="$VEXCL_LDFLAGS $antioch_optional_test_LDFLAGS"
  antioch_optional_test_LIBS="$VEXCL_LIBS $antioch_optional_test_LIBS"
fi

dnl Enable ViennaCL if found
AX_PATH_VIENNACL(0.0.0, no)
if (test x$HAVE_VIENNACL = x1); then
  antioch_optional_test_INCLUDES="$VIENNACL_CPPFLAGS $antioch_optional_test_INCLUDES"
  antioch_optional_test_LDFLAGS="$VIENNACL_LDFLAGS $antioch_optional_test_LDFLAGS"
  antioch_optional_test_LIBS="$VIENNACL_LIBS $antioch_optional_test_LIBS"
fi


dnl Enable Eigen if a new enough version is found
enableoptional=yes
CONFIGURE_EIGEN(3.1.2, no)
if (test x$enableeigen = xyes); then
  antioch_optional_test_INCLUDES="$EIGEN_INCLUDE $antioch_optional_test_INCLUDES"
fi
AM_CONDITIONAL(ANTIOCH_ENABLE_EIGEN, test x$enableeigen = xyes)

dnl Enable MetaPhysicL if a new enough version is found
AX_PATH_METAPHYSICL(0.0.0, no)
if (test x$HAVE_METAPHYSICL = x1); then
  antioch_optional_test_INCLUDES="$METAPHYSICL_CPPFLAGS $antioch_optional_test_INCLUDES"
fi

dnl Use GRVY timers if a new enough version is found
AX_PATH_GRVY_NEW([0.30],[no])
if (test x$HAVE_GRVY = x1); then
  antioch_optional_test_INCLUDES="$GRVY_CFLAGS $antioch_optional_test_INCLUDES"
  antioch_optional_test_LIBS="$GRVY_LIBS $antioch_optional_test_LIBS"
fi

AC_SUBST(antioch_optional_INCLUDES)
AC_SUBST(antioch_optional_test_INCLUDES)
AC_SUBST(antioch_optional_test_LDFLAGS)
AC_SUBST(antioch_optional_test_LIBS)


dnl--------------------------
dnl Even if no dependency needed it,
dnl C++11 makes our asserts nicer,
dnl so let's try and enable it.
dnl--------------------------
AC_ARG_ENABLE(cxx11,
              AC_HELP_STRING([--enable-cxx11],
                             [build with C++11 support]),
	      [case "${enableval}" in
	        yes)  enablecxx11=yes ;;
	         no)  enablecxx11=no ;;
	          *)  AC_MSG_ERROR(bad value ${enableval} for --enable-cxx11) ;;
	       esac],
	       [enablecxx11=optional])

if (test x$enablecxx11 = xyes); then
  AX_CXX_COMPILE_STDCXX_11(noext, mandatory)
  AX_CXX_AUTO_THIS(optional)
elif (test x$enablecxx11 = xoptional); then
  AX_CXX_COMPILE_STDCXX_11(noext, optional)
  AX_CXX_AUTO_THIS(optional)
fi


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

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)

DX_INIT_DOXYGEN(antioch, doxygen/antioch.dox, docs/doxygen)

AC_CONFIG_FILES(doxygen/txt_common/about_vpath.page)

dnl-----------------------------------------------
dnl Generate files for unit and regression testing
dnl-----------------------------------------------
AC_CONFIG_FILES(test/kinetics_unit_air_5sp.sh,           [chmod +x test/kinetics_unit_air_5sp.sh])
AC_CONFIG_FILES(test/kinetics_unit_thrust_air_5sp.sh,    [chmod +x test/kinetics_unit_thrust_air_5sp.sh])
AC_CONFIG_FILES(test/kinetics_vec_unit_air_5sp.sh,       [chmod +x test/kinetics_vec_unit_air_5sp.sh])
AC_CONFIG_FILES(test/kinetics_regression_air_5sp.sh,     [chmod +x test/kinetics_regression_air_5sp.sh])
AC_CONFIG_FILES(test/kinetics_regression_vec_air_5sp.sh, [chmod +x test/kinetics_regression_vec_air_5sp.sh])
AC_CONFIG_FILES(test/photochemical_rate_unit.sh,         [chmod +x test/photochemical_rate_unit.sh])
AC_CONFIG_FILES(test/sigma_bin_converter_unit.sh,        [chmod +x test/sigma_bin_converter_unit.sh])
AC_CONFIG_FILES(test/parsing_xml.sh,                     [chmod +x test/parsing_xml.sh])
AC_CONFIG_FILES(test/kinetics_reactive_scheme_unit.sh,   [chmod +x test/kinetics_reactive_scheme_unit.sh])

dnl-----------------------------------------------
dnl Generate header files
dnl-----------------------------------------------
AC_CONFIG_FILES(src/utilities/include/antioch/antioch_version.h)

dnl-----------------------------------------------
dnl Generate Makefiles and other templated things
dnl-----------------------------------------------
AC_CONFIG_FILES([
  antioch.pc
  doxygen/Makefile
  doxygen/antioch.dox
  Makefile
  src/Makefile
  test/Makefile
])

dnl-----------------------------------------------
dnl Collect files for "header stamping" here
dnl-----------------------------------------------

# PB: Although the documentation says that ac_top_srcdir should be defined,
#     it doesn't get activated until _AC_SRCDIRS gets called
#     (used to be _AC_SRCPATHS), according to this thread:
#     http://lists.gnu.org/archive/html/autoconf-patches/2003-02/msg00003.html
#     My own hacking revealed that we could also do the following:
#     TOP_SEARCH_DIR=$ac_pwd/$srcdir
#     I'm not sure which will be more "future proof", but this is the alternative
#     to using $(shell) calls in the Makefile.am which 1. requires
#     GNU Make and 2. Causes automake to issue a warning
_AC_SRCDIRS(.)
TOP_SEARCH_DIR=$ac_abs_top_srcdir

# We have to do this by subdirectory because otherwise distcheck
# breaks as we start picking up files in the directories
# that it uses.
ANTIOCH_STAMPED_FILES=$(find $TOP_SEARCH_DIR/{src,test} -name "*.h" -or -name "*.C" | tr "\n" " ")
AC_SUBST(ANTIOCH_STAMPED_FILES)

# defaults files, set here for tests
# changed in Makefile.am at install
# step (see install-data-hook)
# to make it $datadir/antioch_default_files
# then we define the default files
ANTIOCH_DEFAULT_FILES_PATH_TMP=$TOP_SEARCH_DIR/share/antioch_default_files/
AC_DEFINE_UNQUOTED(DEFAULT_FILES_PATH,"$ANTIOCH_DEFAULT_FILES_PATH_TMP",[Path to defaults files])
AC_DEFINE(DEFAULT_SPECIES_LIST,"antioch_default_chemical_species.dat",[Default list of species])
AC_DEFINE(DEFAULT_CHEMICAL_MIXTURE,"antioch_default_chemical_mixture.dat",[Mandatory default informations on species])
AC_DEFINE(DEFAULT_VIBRATIONAL_DATA,"antioch_default_vibrational_data.dat",[Vibrational default data])
AC_DEFINE(DEFAULT_ELECTRONIC_DATA,"antioch_default_electronic_data.dat",[Electronic default data])


# Must still call AC_OUTPUT() after generating all the files
AC_OUTPUT()

dnl------------------------------
dnl Final summary
dnl------------------------------
AX_SUMMARIZE_CONFIG
