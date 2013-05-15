BUILT_SOURCES   = .license.stamp

EXTRA_DIST = input_files

check_PROGRAMS  =
check_PROGRAMS += chem_mixture_unit
check_PROGRAMS += chem_mixture_vec_unit
check_PROGRAMS += cea_evaluator_unit
check_PROGRAMS += cea_evaluator_vec_unit
check_PROGRAMS += cea_thermo_unit
check_PROGRAMS += cea_thermo_vec_unit
check_PROGRAMS += arrhenius_rate_unit
check_PROGRAMS += arrhenius_rate_vec_unit
check_PROGRAMS += kinetics_unit
check_PROGRAMS += kinetics_vec_unit
check_PROGRAMS += kinetics_regression
check_PROGRAMS += kinetics_regression_vec
check_PROGRAMS += blottner_viscosity_unit
check_PROGRAMS += blottner_viscosity_vec_unit
check_PROGRAMS += sutherland_viscosity_unit
check_PROGRAMS += sutherland_viscosity_vec_unit
check_PROGRAMS += mixture_viscosity_unit
check_PROGRAMS += mixture_viscosity_vec_unit
check_PROGRAMS += stat_mech_thermo_unit
check_PROGRAMS += eucken_thermal_cond_unit
check_PROGRAMS += wilke_transport_unit
check_PROGRAMS += wilke_transport_vec_unit
check_PROGRAMS += constant_lewis_unit

if ANTIOCH_ENABLE_EIGEN
check_PROGRAMS += stat_mech_thermo_unit_eigen
endif

AM_CPPFLAGS  = 
AM_CPPFLAGS += -I$(top_srcdir)/src/core/include
AM_CPPFLAGS += -I$(top_srcdir)/src/kinetics/include
AM_CPPFLAGS += -I$(top_srcdir)/src/parsing/include
AM_CPPFLAGS += -I$(top_srcdir)/src/thermo/include
AM_CPPFLAGS += -I$(top_srcdir)/src/transport/include
AM_CPPFLAGS += -I$(top_srcdir)/src/utilities/include
AM_CPPFLAGS += -I$(top_builddir)/src/utilities/include #antioch_version.h
AM_CPPFLAGS +=  $(antioch_optional_test_INCLUDES)

AM_LDFLAGS = $(top_builddir)/src/libantioch.la

# Sources for these tests
chem_mixture_unit_SOURCES = chem_mixture_unit.C
chem_mixture_vec_unit_SOURCES = chem_mixture_vec_unit.C
cea_evaluator_unit_SOURCES = cea_evaluator_unit.C
cea_evaluator_vec_unit_SOURCES = cea_evaluator_vec_unit.C
cea_thermo_unit_SOURCES = cea_thermo_unit.C
cea_thermo_vec_unit_SOURCES = cea_thermo_vec_unit.C
arrhenius_rate_unit_SOURCES = arrhenius_rate_unit.C
arrhenius_rate_vec_unit_SOURCES = arrhenius_rate_vec_unit.C
kinetics_unit_SOURCES = kinetics_unit.C
kinetics_vec_unit_SOURCES = kinetics_vec_unit.C
kinetics_regression_SOURCES = kinetics_regression.C
kinetics_regression_vec_SOURCES = kinetics_regression_vec.C
blottner_viscosity_unit_SOURCES = blottner_viscosity_unit.C
blottner_viscosity_vec_unit_SOURCES = blottner_viscosity_vec_unit.C
sutherland_viscosity_unit_SOURCES = sutherland_viscosity_unit.C
sutherland_viscosity_vec_unit_SOURCES = sutherland_viscosity_vec_unit.C
mixture_viscosity_unit_SOURCES = mixture_viscosity_unit.C
mixture_viscosity_vec_unit_SOURCES = mixture_viscosity_vec_unit.C
stat_mech_thermo_unit_SOURCES = stat_mech_thermo_unit.C
eucken_thermal_cond_unit_SOURCES = eucken_thermal_cond_unit.C
wilke_transport_unit_SOURCES = wilke_transport_unit.C
wilke_transport_vec_unit_SOURCES = wilke_transport_vec_unit.C
constant_lewis_unit_SOURCES = constant_lewis_unit.C

if ANTIOCH_ENABLE_EIGEN
stat_mech_thermo_unit_eigen_SOURCES = stat_mech_thermo_unit_eigen.C
endif

#Define tests to actually be run
TESTS  =
TESTS += chem_mixture_unit
TESTS += chem_mixture_vec_unit
TESTS += cea_evaluator_unit
TESTS += cea_evaluator_vec_unit
TESTS += cea_thermo_unit
TESTS += cea_thermo_vec_unit
TESTS += arrhenius_rate_unit
TESTS += arrhenius_rate_vec_unit
TESTS += kinetics_unit_air_5sp.sh
TESTS += kinetics_vec_unit_air_5sp.sh
TESTS += kinetics_regression_air_5sp.sh
TESTS += kinetics_regression_vec_air_5sp.sh
TESTS += blottner_viscosity_unit
TESTS += blottner_viscosity_vec_unit
TESTS += sutherland_viscosity_unit
TESTS += sutherland_viscosity_vec_unit
TESTS += mixture_viscosity_unit
TESTS += mixture_viscosity_vec_unit
TESTS += stat_mech_thermo_unit
TESTS += eucken_thermal_cond_unit
TESTS += wilke_transport_unit
TESTS += wilke_transport_vec_unit
TESTS += constant_lewis_unit

if ANTIOCH_ENABLE_EIGEN
TESTS += stat_mech_thermo_unit_eigen
endif


CLEANFILES =
if CODE_COVERAGE_ENABLED
  CLEANFILES += *.gcda *.gcno
endif

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE \
	$(top_srcdir)/test/*.C
	echo 'updated source license headers' >$@

CLEANFILES += .license.stamp

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
