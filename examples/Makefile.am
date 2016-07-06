EXTRA_DIST = input_files

check_PROGRAMS  =
check_PROGRAMS += antioch_init

AM_CPPFLAGS  =
AM_CPPFLAGS += -I$(top_srcdir)/src/core/include
AM_CPPFLAGS += -I$(top_srcdir)/src/units/include
AM_CPPFLAGS += -I$(top_srcdir)/src/particles_flux/include
AM_CPPFLAGS += -I$(top_srcdir)/src/kinetics/include
AM_CPPFLAGS += -I$(top_srcdir)/src/parsing/include
AM_CPPFLAGS += -I$(top_srcdir)/src/thermo/include
AM_CPPFLAGS += -I$(top_srcdir)/src/viscosity/include
AM_CPPFLAGS += -I$(top_srcdir)/src/diffusion/include
AM_CPPFLAGS += -I$(top_srcdir)/src/thermal_conduction/include
AM_CPPFLAGS += -I$(top_srcdir)/src/transport/include
AM_CPPFLAGS += -I$(top_srcdir)/src/utilities/include
AM_CPPFLAGS += -I$(top_builddir)/src/utilities/include #antioch_version.h
AM_CPPFLAGS += -I$(top_srcdir)/test/standard_unit
AM_CPPFLAGS +=  $(antioch_optional_test_INCLUDES)

AM_LDFLAGS = $(antioch_optional_test_LDFLAGS)

LIBS = $(antioch_optional_test_LIBS)

LDADD = $(top_builddir)/src/libantioch.la

pkginclude_HEADERS =

standard_unit_tests_SOURCES = antioch_init.C

#Define examples to actually be run
TESTS  =
XFAIL_TESTS  =
TESTS += antioch_init

CLEANFILES =
if CODE_COVERAGE_ENABLED
  CLEANFILES += *.gcda *.gcno
endif
