BUILT_SOURCES   = .license.stamp

#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

bin_PROGRAMS    = antioch_version 

#----------------------------
# INCLUDES we want distributed
#----------------------------
pkginclude_HEADERS  =

# core
pkginclude_HEADERS += core/include/antioch/chemical_species.h
pkginclude_HEADERS += core/include/antioch/chemical_mixture.h
pkginclude_HEADERS += core/include/antioch/species_enum.h

# unit
pkginclude_HEADERS += units/include/antioch/converter.h
pkginclude_HEADERS += units/include/antioch/insi.h
pkginclude_HEADERS += units/include/antioch/siprefix.h
pkginclude_HEADERS += units/include/antioch/unit_base.h
pkginclude_HEADERS += units/include/antioch/unit_defs.h
pkginclude_HEADERS += units/include/antioch/units.h
pkginclude_HEADERS += units/include/antioch/unit_store.h

# particles_flux
pkginclude_HEADERS += particles_flux/include/antioch/sigma_bin_converter.h
pkginclude_HEADERS += particles_flux/include/antioch/particle_flux.h

# kinetics-kinetics models
pkginclude_HEADERS += kinetics/include/antioch/kinetics_enum.h
pkginclude_HEADERS += kinetics/include/antioch/kinetics_type.h
pkginclude_HEADERS += kinetics/include/antioch/hercourtessen_rate.h
pkginclude_HEADERS += kinetics/include/antioch/berthelot_rate.h
pkginclude_HEADERS += kinetics/include/antioch/arrhenius_rate.h
pkginclude_HEADERS += kinetics/include/antioch/berthelothercourtessen_rate.h
pkginclude_HEADERS += kinetics/include/antioch/kooij_rate.h
pkginclude_HEADERS += kinetics/include/antioch/vanthoff_rate.h
pkginclude_HEADERS += kinetics/include/antioch/photochemical_rate.h
# kinetics-chemical process
pkginclude_HEADERS += kinetics/include/antioch/reaction_enum.h
pkginclude_HEADERS += kinetics/include/antioch/reaction.h
pkginclude_HEADERS += kinetics/include/antioch/elementary_reaction.h
pkginclude_HEADERS += kinetics/include/antioch/duplicate_reaction.h
pkginclude_HEADERS += kinetics/include/antioch/threebody_reaction.h
pkginclude_HEADERS += kinetics/include/antioch/falloff_reaction.h
pkginclude_HEADERS += kinetics/include/antioch/lindemann_falloff.h
pkginclude_HEADERS += kinetics/include/antioch/troe_falloff.h
# kinetics-other
pkginclude_HEADERS += kinetics/include/antioch/reaction_set.h
pkginclude_HEADERS += kinetics/include/antioch/reaction_parsing.h
pkginclude_HEADERS += kinetics/include/antioch/kinetics_parsing.h
pkginclude_HEADERS += kinetics/include/antioch/kinetics_evaluator.h

# parsing
pkginclude_HEADERS += parsing/include/antioch/tinyxml2.h
pkginclude_HEADERS += parsing/include/antioch/tinyxml2_imp.h
pkginclude_HEADERS += parsing/include/antioch/read_reaction_set_data_xml.h
pkginclude_HEADERS += parsing/include/antioch/blottner_parsing.h
pkginclude_HEADERS += parsing/include/antioch/sutherland_parsing.h
pkginclude_HEADERS += parsing/include/antioch/cea_thermo_ascii_parsing.h
pkginclude_HEADERS += parsing/include/antioch/cea_mixture_ascii_parsing.h
pkginclude_HEADERS += parsing/include/antioch/species_ascii_parsing.h

# thermo
pkginclude_HEADERS += thermo/include/antioch/cea_curve_fit.h
pkginclude_HEADERS += thermo/include/antioch/cea_mixture.h
pkginclude_HEADERS += thermo/include/antioch/temp_cache.h
pkginclude_HEADERS += thermo/include/antioch/cea_evaluator.h
pkginclude_HEADERS += thermo/include/antioch/cea_thermo.h
pkginclude_HEADERS += thermo/include/antioch/stat_mech_thermo.h

# transport
pkginclude_HEADERS += transport/include/antioch/blottner_viscosity.h
pkginclude_HEADERS += transport/include/antioch/sutherland_viscosity.h
pkginclude_HEADERS += transport/include/antioch/mixture_viscosity.h
pkginclude_HEADERS += transport/include/antioch/wilke_mixture.h
pkginclude_HEADERS += transport/include/antioch/eucken_thermal_conductivity.h
pkginclude_HEADERS += transport/include/antioch/wilke_evaluator.h
pkginclude_HEADERS += transport/include/antioch/constant_lewis_diffusivity.h

# utilities
pkginclude_HEADERS += utilities/include/antioch/antioch_exceptions.h
pkginclude_HEADERS += utilities/include/antioch/antioch_asserts.h
pkginclude_HEADERS += utilities/include/antioch/cmath_shims.h
pkginclude_HEADERS += utilities/include/antioch/eigen_utils.h
pkginclude_HEADERS += utilities/include/antioch/eigen_utils_decl.h
pkginclude_HEADERS += utilities/include/antioch/input_utils.h
pkginclude_HEADERS += utilities/include/antioch/math_constants.h
pkginclude_HEADERS += utilities/include/antioch/metaprogramming.h
pkginclude_HEADERS += utilities/include/antioch/metaprogramming_decl.h
pkginclude_HEADERS += utilities/include/antioch/metaphysicl_utils.h
pkginclude_HEADERS += utilities/include/antioch/metaphysicl_utils_decl.h
pkginclude_HEADERS += utilities/include/antioch/physical_constants.h
pkginclude_HEADERS += utilities/include/antioch/string_utils.h
pkginclude_HEADERS += utilities/include/antioch/valarray_utils.h
pkginclude_HEADERS += utilities/include/antioch/valarray_utils_decl.h
pkginclude_HEADERS += utilities/include/antioch/vector_utils.h
pkginclude_HEADERS += utilities/include/antioch/vector_utils_decl.h
pkginclude_HEADERS += utilities/include/antioch/vexcl_utils.h
pkginclude_HEADERS += utilities/include/antioch/vexcl_utils_decl.h

# Needs to be builddir since this is generated by configure
pkginclude_HEADERS += $(top_builddir)/src/utilities/include/antioch/antioch_version.h

# Version app
antioch_version_SOURCES = apps/version.C

#--------------------------------------
#Local Directories to include for build
#--------------------------------------
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
AM_CPPFLAGS +=  $(antioch_optional_INCLUDES)

CLEANFILES =

#--------------
# Code coverage
#--------------
if CODE_COVERAGE_ENABLED
  CLEANFILES += *.gcda *.gcno
endif

#---------------------------------
# Embedded license header support
#---------------------------------
STAMPED_FILES = $(pkginclude_HEADERS) $(antioch_version_SOURCES)

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl -S=$(top_srcdir)/src $(top_srcdir)/LICENSE $(STAMPED_FILES) 
	echo 'updated source license headers' >$@

CLEANFILES += .license.stamp

EXTRA_DIST = common

