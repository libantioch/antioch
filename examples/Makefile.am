EXTRA_DIST = 

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
AM_CPPFLAGS += -I$(top_builddir)/src/utilities/include 
AM_CPPFLAGS += -I$(top_srcdir)/examples/antioch_test
AM_CPPFLAGS +=  $(antioch_optional_test_INCLUDES)

AM_LDFLAGS = $(antioch_optional_test_LDFLAGS)

LIBS = $(antioch_optional_test_LIBS)

LDADD = $(top_builddir)/src/libantioch.la

pkginclude_HEADERS =

#
# CPP TESTS
#
antioch_initdir = $(prefix)/share/examples/antioch_init
antioch_init_PROGRAMS = antioch_init
antioch_init_SOURCES = antioch_init.C
antioch_init_DATA    = ${antioch_init_SOURCES}

#Define examples to actually be run
TESTS  =
XFAIL_TESTS  =
TESTS += antioch_init

CLEANFILES =
if CODE_COVERAGE_ENABLED
  CLEANFILES += *.gcda *.gcno
endif
