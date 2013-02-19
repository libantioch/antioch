BUILT_SOURCES   = .license.stamp

check_PROGRAMS  =
check_PROGRAMS += chem_mixture_unit
check_PROGRAMS += chem_mixture_vec_unit
check_PROGRAMS += cea_thermo_unit
check_PROGRAMS += cea_thermo_vec_unit
check_PROGRAMS += arrhenius_rate_unit
check_PROGRAMS += arrhenius_rate_vec_unit
check_PROGRAMS += kinetics_unit
check_PROGRAMS += kinetics_regression
check_PROGRAMS += blottner_viscosity_unit
check_PROGRAMS += blottner_viscosity_vec_unit
check_PROGRAMS += sutherland_viscosity_unit
check_PROGRAMS += mixture_viscosity_unit

AM_CPPFLAGS  = 
AM_CPPFLAGS += -I$(top_srcdir)/src/core/include
AM_CPPFLAGS += -I$(top_srcdir)/src/kinetics/include
AM_CPPFLAGS += -I$(top_srcdir)/src/parsing/include
AM_CPPFLAGS += -I$(top_srcdir)/src/thermo/include
AM_CPPFLAGS += -I$(top_srcdir)/src/transport/include
AM_CPPFLAGS += -I$(top_srcdir)/src/utilities/include
AM_CPPFLAGS += -I$(top_builddir)/src/utilities/include #antioch_version.h

AM_LDFLAGS = $(top_builddir)/src/libantioch.la

# Sources for these tests
chem_mixture_unit_SOURCES = chem_mixture_unit.C
chem_mixture_vec_unit_SOURCES = chem_mixture_vec_unit.C
cea_thermo_unit_SOURCES = cea_thermo_unit.C
cea_thermo_vec_unit_SOURCES = cea_thermo_vec_unit.C
arrhenius_rate_unit_SOURCES = arrhenius_rate_unit.C
arrhenius_rate_vec_unit_SOURCES = arrhenius_rate_vec_unit.C
kinetics_unit_SOURCES = kinetics_unit.C
kinetics_regression_SOURCES = kinetics_regression.C
blottner_viscosity_unit_SOURCES = blottner_viscosity_unit.C
blottner_viscosity_vec_unit_SOURCES = blottner_viscosity_vec_unit.C
sutherland_viscosity_unit_SOURCES = sutherland_viscosity_unit.C
mixture_viscosity_unit_SOURCES = mixture_viscosity_unit.C

#Define tests to actually be run
TESTS  =
TESTS += chem_mixture_unit
TESTS += chem_mixture_vec_unit
TESTS += cea_thermo_unit
TESTS += cea_thermo_vec_unit
TESTS += arrhenius_rate_unit
TESTS += arrhenius_rate_vec_unit
TESTS += kinetics_unit_air_5sp.sh
TESTS += kinetics_regression_air_5sp.sh
TESTS += blottner_viscosity_unit
TESTS += blottner_viscosity_vec_unit
TESTS += sutherland_viscosity_unit
TESTS += mixture_viscosity_unit

if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE \
	$(top_srcdir)/test/*.C
	echo 'updated source license headers' >$@

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
