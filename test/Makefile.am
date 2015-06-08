EXTRA_DIST = input_files

check_PROGRAMS  =
check_PROGRAMS += chem_mixture_unit
check_PROGRAMS += chem_mixture_vec_unit
check_PROGRAMS += nasa_evaluator_unit
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
check_PROGRAMS += mixture_conductivity_regression
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
check_PROGRAMS += kinetics_theory_thermal_conductivity_unit
check_PROGRAMS += kinetics_theory_thermal_conductivity_vec_unit
check_PROGRAMS += rotational_relaxation_unit
check_PROGRAMS += rotational_relaxation_vec_unit
check_PROGRAMS += ideal_gas_micro_thermo_unit
check_PROGRAMS += ascii_parser_unit
check_PROGRAMS += lindemann_falloff_threebody_unit
check_PROGRAMS += troe_falloff_threebody_unit

#GSL Tests
check_PROGRAMS += molecular_binary_diffusion_unit
check_PROGRAMS += molecular_binary_diffusion_vec_unit
check_PROGRAMS += kinetics_theory_viscosity_unit
check_PROGRAMS += kinetics_theory_viscosity_vec_unit
check_PROGRAMS += gsl_spline_unit
check_PROGRAMS += gsl_spline_vec_unit
check_PROGRAMS += Stockmayer_unit

# Eigen Tests
check_PROGRAMS += stat_mech_thermo_unit_eigen

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
AM_CPPFLAGS +=  $(antioch_optional_test_INCLUDES)

AM_LDFLAGS = $(antioch_optional_test_LDFLAGS)
LIBS = $(antioch_optional_test_LIBS)

LDADD = $(top_builddir)/src/libantioch.la

# Sources for these tests
chem_mixture_unit_SOURCES = chem_mixture_unit.C
chem_mixture_vec_unit_SOURCES = chem_mixture_vec_unit.C
nasa_evaluator_unit_SOURCES = nasa_evaluator_unit.C
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
mixture_conductivity_regression_SOURCES = mixture_conductivity_regression.C
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
kinetics_theory_thermal_conductivity_unit_SOURCES = kinetics_theory_thermal_conductivity_unit.C
kinetics_theory_thermal_conductivity_vec_unit_SOURCES = kinetics_theory_thermal_conductivity_vec_unit.C
rotational_relaxation_unit_SOURCES = rotational_relaxation_unit.C
rotational_relaxation_vec_unit_SOURCES = rotational_relaxation_vec_unit.C
ascii_parser_unit_SOURCES = ascii_parser_unit.C
lindemann_falloff_threebody_unit_SOURCES = lindemann_falloff_threebody_unit.C
troe_falloff_threebody_unit_SOURCES = troe_falloff_threebody_unit.C

# GSL Tests
molecular_binary_diffusion_unit_SOURCES = molecular_binary_diffusion_unit.C
molecular_binary_diffusion_vec_unit_SOURCES = molecular_binary_diffusion_vec_unit.C
kinetics_theory_viscosity_unit_SOURCES = kinetics_theory_viscosity_unit.C
kinetics_theory_viscosity_vec_unit_SOURCES = kinetics_theory_viscosity_vec_unit.C
gsl_spline_unit_SOURCES = gsl_spline_unit.C
gsl_spline_vec_unit_SOURCES = gsl_spline_vec_unit.C
ideal_gas_micro_thermo_unit_SOURCES =  ideal_gas_micro_thermo_unit.C
Stockmayer_unit_SOURCES = Stockmayer_unit.C


# Eigen Tests
stat_mech_thermo_unit_eigen_SOURCES = stat_mech_thermo_unit_eigen.C


#Define tests to actually be run
TESTS  =
TESTS += chem_mixture_unit
TESTS += chem_mixture_vec_unit
TESTS += nasa_evaluator_unit.sh
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
TESTS += mixture_conductivity_regression
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
TESTS += kinetics_settings_unit
TESTS += kinetics_theory_thermal_conductivity_unit
TESTS += kinetics_theory_thermal_conductivity_vec_unit
TESTS += rotational_relaxation_unit
TESTS += rotational_relaxation_vec_unit
TESTS += ideal_gas_micro_thermo_unit
TESTS += ascii_parser_unit.sh
TESTS += lindemann_falloff_threebody_unit
TESTS += troe_falloff_threebody_unit

# GSL Tests
TESTS += molecular_binary_diffusion_unit
TESTS += molecular_binary_diffusion_vec_unit
TESTS += kinetics_theory_viscosity_unit
TESTS += kinetics_theory_viscosity_vec_unit
TESTS += gsl_spline_unit
TESTS += gsl_spline_vec_unit
TESTS += Stockmayer_unit


# Eigen Tests
TESTS += stat_mech_thermo_unit_eigen


CLEANFILES =
if CODE_COVERAGE_ENABLED
  CLEANFILES += *.gcda *.gcno
endif

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
