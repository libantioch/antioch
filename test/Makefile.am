EXTRA_DIST = input_files

check_PROGRAMS  =
check_PROGRAMS += chem_mixture_unit
check_PROGRAMS += chem_mixture_vec_unit
check_PROGRAMS += cea_evaluator_unit
check_PROGRAMS += cea_evaluator_vec_unit
check_PROGRAMS += cea_thermo_unit
check_PROGRAMS += cea_thermo_vec_unit
check_PROGRAMS += constant_rate_unit
check_PROGRAMS += constant_rate_vec_unit
check_PROGRAMS += hercourtessen_rate_unit
check_PROGRAMS += hercourtessen_rate_vec_unit
check_PROGRAMS += berthelot_rate_unit
check_PROGRAMS += berthelot_rate_vec_unit
check_PROGRAMS += arrhenius_rate_deriv
check_PROGRAMS += arrhenius_rate_unit
check_PROGRAMS += arrhenius_rate_vec_unit
check_PROGRAMS += kooij_rate_unit
check_PROGRAMS += kooij_rate_vec_unit
check_PROGRAMS += berthelothercourtessen_rate_unit
check_PROGRAMS += berthelothercourtessen_rate_vec_unit
check_PROGRAMS += vanthoff_rate_unit
check_PROGRAMS += vanthoff_rate_vec_unit
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
check_PROGRAMS += elementary_process_unit
check_PROGRAMS += duplicate_process_unit
check_PROGRAMS += threebody_process_unit
check_PROGRAMS += lindemann_falloff_unit
check_PROGRAMS += troe_falloff_unit
check_PROGRAMS += photochemical_rate_unit
check_PROGRAMS += kinetics_reversibility_unit
check_PROGRAMS += parsing_xml
check_PROGRAMS += parsing_chemkin
check_PROGRAMS += units_unit
check_PROGRAMS += kinetics_reactive_scheme_unit
check_PROGRAMS += sigma_bin_converter_unit
check_PROGRAMS += sigma_bin_converter_vec_unit
check_PROGRAMS += kinetics_settings_unit

if ANTIOCH_ENABLE_EIGEN
check_PROGRAMS += stat_mech_thermo_unit_eigen
endif

AM_CPPFLAGS  = 
AM_CPPFLAGS += -I$(top_srcdir)/src/core/include
AM_CPPFLAGS += -I$(top_srcdir)/src/units/include
AM_CPPFLAGS += -I$(top_srcdir)/src/particles_flux/include
AM_CPPFLAGS += -I$(top_srcdir)/src/kinetics/include
AM_CPPFLAGS += -I$(top_srcdir)/src/parsing/include
AM_CPPFLAGS += -I$(top_srcdir)/src/thermo/include
AM_CPPFLAGS += -I$(top_srcdir)/src/transport/include
AM_CPPFLAGS += -I$(top_srcdir)/src/utilities/include
AM_CPPFLAGS += -I$(top_builddir)/src/utilities/include #antioch_version.h
AM_CPPFLAGS +=  $(antioch_optional_test_INCLUDES)

AM_LDFLAGS = $(antioch_optional_test_LDFLAGS)
LIBS = $(antioch_optional_test_LIBS)

# Sources for these tests
chem_mixture_unit_SOURCES = chem_mixture_unit.C
chem_mixture_vec_unit_SOURCES = chem_mixture_vec_unit.C
cea_evaluator_unit_SOURCES = cea_evaluator_unit.C
cea_evaluator_vec_unit_SOURCES = cea_evaluator_vec_unit.C
cea_thermo_unit_SOURCES = cea_thermo_unit.C
cea_thermo_vec_unit_SOURCES = cea_thermo_vec_unit.C
constant_rate_unit_SOURCES = constant_rate_unit.C
constant_rate_vec_unit_SOURCES = constant_rate_vec_unit.C
hercourtessen_rate_unit_SOURCES = hercourtessen_rate_unit.C
hercourtessen_rate_vec_unit_SOURCES = hercourtessen_rate_vec_unit.C
berthelot_rate_unit_SOURCES = berthelot_rate_unit.C
berthelot_rate_vec_unit_SOURCES = berthelot_rate_vec_unit.C
arrhenius_rate_deriv_SOURCES = arrhenius_rate_deriv.C
arrhenius_rate_unit_SOURCES = arrhenius_rate_unit.C
arrhenius_rate_vec_unit_SOURCES = arrhenius_rate_vec_unit.C
berthelothercourtessen_rate_unit_SOURCES = berthelothercourtessen_rate_unit.C
berthelothercourtessen_rate_vec_unit_SOURCES = berthelothercourtessen_rate_vec_unit.C
kooij_rate_unit_SOURCES = kooij_rate_unit.C
kooij_rate_vec_unit_SOURCES = kooij_rate_vec_unit.C
vanthoff_rate_unit_SOURCES = vanthoff_rate_unit.C
vanthoff_rate_vec_unit_SOURCES = vanthoff_rate_vec_unit.C
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
elementary_process_unit_SOURCES = elementary_process_unit.C
duplicate_process_unit_SOURCES = duplicate_process_unit.C
threebody_process_unit_SOURCES = threebody_process_unit.C
lindemann_falloff_unit_SOURCES = lindemann_falloff_unit.C
troe_falloff_unit_SOURCES = troe_falloff_unit.C
photochemical_rate_unit_SOURCES = photochemical_rate_unit.C
kinetics_reversibility_unit_SOURCES = kinetics_reversibility_unit.C
parsing_xml_SOURCES = parsing_xml.C
parsing_chemkin_SOURCES = parsing_chemkin.C
units_unit_SOURCES = units_unit.C
kinetics_reactive_scheme_unit_SOURCES = kinetics_reactive_scheme_unit.C
sigma_bin_converter_unit_SOURCES = sigma_bin_converter_unit.C
sigma_bin_converter_vec_unit_SOURCES = sigma_bin_converter_vec_unit.C
kinetics_settings_unit_SOURCES = kinetics_settings_unit.C

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
TESTS += constant_rate_unit
TESTS += constant_rate_vec_unit
TESTS += hercourtessen_rate_unit
TESTS += hercourtessen_rate_unit
TESTS += berthelot_rate_vec_unit
TESTS += berthelot_rate_vec_unit
TESTS += arrhenius_rate_deriv
TESTS += arrhenius_rate_unit
TESTS += arrhenius_rate_vec_unit
TESTS += berthelothercourtessen_rate_unit
TESTS += berthelothercourtessen_rate_vec_unit
TESTS += kooij_rate_unit
TESTS += kooij_rate_vec_unit
TESTS += vanthoff_rate_unit
TESTS += vanthoff_rate_vec_unit
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
TESTS += elementary_process_unit
TESTS += duplicate_process_unit
TESTS += threebody_process_unit
TESTS += lindemann_falloff_unit
TESTS += troe_falloff_unit
TESTS += photochemical_rate_unit.sh
TESTS += kinetics_reversibility_unit
TESTS += parsing_xml.sh
TESTS += parsing_chemkin.sh
TESTS += units_unit
TESTS += kinetics_reactive_scheme_unit.sh
TESTS += sigma_bin_converter_unit
TESTS += sigma_bin_converter_vec_unit
TESTS += sigma_bin_converter_unit.sh
TESTS += kinetics_settings_unit

if ANTIOCH_ENABLE_EIGEN
TESTS += stat_mech_thermo_unit_eigen
endif

CLEANFILES =
if CODE_COVERAGE_ENABLED
  CLEANFILES += *.gcda *.gcno
endif

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
