BUILT_SOURCES   = .license.stamp

#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

bin_PROGRAMS    = antioch_version 

lib_LTLIBRARIES = libantioch.la

#----------------------------------------------
# List of source files to build dynamic library
#----------------------------------------------
libantioch_la_SOURCES =

# core
libantioch_la_SOURCES += core/src/chemical_species.C
libantioch_la_SOURCES += core/src/chemical_mixture.C

# preproc/core
libantioch_la_SOURCES += preproc/core/src/Error.cpp
libantioch_la_SOURCES += preproc/core/src/SIPrefix.cpp
libantioch_la_SOURCES += preproc/core/src/Converter.cpp
libantioch_la_SOURCES += preproc/core/src/InSI.cpp
libantioch_la_SOURCES += preproc/core/src/Units.cpp
libantioch_la_SOURCES += preproc/core/src/CoreUnc.cpp
libantioch_la_SOURCES += preproc/core/src/AtomicParameter.cpp
libantioch_la_SOURCES += preproc/core/src/PDF.cpp
libantioch_la_SOURCES += preproc/core/src/ParameterPhy.cpp
libantioch_la_SOURCES += preproc/core/src/phymath.cpp
libantioch_la_SOURCES += preproc/core/src/PhyFunct.cpp
libantioch_la_SOURCES += preproc/core/src/PhyHarmFunct.cpp
# preproc/chemistry
libantioch_la_SOURCES += preproc/chemistry/src/allKineticsProcesses.cpp
libantioch_la_SOURCES += preproc/chemistry/src/kineticsCrossRoad.cpp
libantioch_la_SOURCES += preproc/chemistry/src/SimpleChem.cpp
libantioch_la_SOURCES += preproc/chemistry/src/KineticsChem.cpp
libantioch_la_SOURCES += preproc/chemistry/src/ReactionChem.cpp
# preproc/kinetics/process
libantioch_la_SOURCES += preproc/kinetics/process/src/Process.cpp
libantioch_la_SOURCES += preproc/kinetics/process/src/ElementaryProcess.cpp
libantioch_la_SOURCES += preproc/kinetics/process/src/DuplicateProcess.cpp
libantioch_la_SOURCES += preproc/kinetics/process/src/ThreeBody.cpp
libantioch_la_SOURCES += preproc/kinetics/process/src/FalloffProcess.cpp
libantioch_la_SOURCES += preproc/kinetics/process/src/LindemannFalloff.cpp
libantioch_la_SOURCES += preproc/kinetics/process/src/TroeFalloff.cpp
# preproc/kinetics/model
libantioch_la_SOURCES += preproc/kinetics/model/src/KineticsModel.cpp
libantioch_la_SOURCES += preproc/kinetics/model/src/Arrhenius.cpp
libantioch_la_SOURCES += preproc/kinetics/model/src/HercourtEssen.cpp
libantioch_la_SOURCES += preproc/kinetics/model/src/Berthelot.cpp
libantioch_la_SOURCES += preproc/kinetics/model/src/BerthelotHercourtEssen.cpp
libantioch_la_SOURCES += preproc/kinetics/model/src/MultHercourtEssen.cpp
libantioch_la_SOURCES += preproc/kinetics/model/src/Kooij.cpp
libantioch_la_SOURCES += preproc/kinetics/model/src/VantHoff.cpp

# kinetics
libantioch_la_SOURCES += kinetics/src/reaction.C
libantioch_la_SOURCES += kinetics/src/reaction_set.C
libantioch_la_SOURCES += kinetics/src/kinetics_evaluator.C

# parsing
libantioch_la_SOURCES += parsing/src/tinyxml2.cpp
libantioch_la_SOURCES += parsing/src/read_reaction_set_data_xml.C
libantioch_la_SOURCES += parsing/src/blottner_parsing.C
libantioch_la_SOURCES += parsing/src/sutherland_parsing.C
libantioch_la_SOURCES += parsing/src/cea_thermo_ascii_parsing.C
libantioch_la_SOURCES += parsing/src/species_ascii_parsing.C

# thermo
libantioch_la_SOURCES += thermo/src/cea_curve_fit.C
libantioch_la_SOURCES += thermo/src/cea_thermo.C
libantioch_la_SOURCES += thermo/src/stat_mech_thermo.C

# transport
libantioch_la_SOURCES += transport/src/wilke_mixture.C

# utilities
libantioch_la_SOURCES += utilities/src/antioch_version.C
libantioch_la_SOURCES += utilities/src/input_utils.C


#----------------------------
# INCLUDES we want distributed
#----------------------------
includedir = $(prefix)/include/antioch
include_HEADERS =

# core
include_HEADERS += core/include/antioch/chemical_species.h
include_HEADERS += core/include/antioch/chemical_mixture.h
include_HEADERS += core/include/antioch/species_enum.h

# preproc/core
include_HEADERS += preproc/core/include/antioch/Error.hpp
include_HEADERS += preproc/core/include/antioch/SIPrefix.hpp
include_HEADERS += preproc/core/include/antioch/Converter.hpp
include_HEADERS += preproc/core/include/antioch/InSI.hpp
include_HEADERS += preproc/core/include/antioch/Units.hpp
include_HEADERS += preproc/core/include/antioch/unit_defs.hpp
include_HEADERS += preproc/core/include/antioch/CoreUnc.hpp
include_HEADERS += preproc/core/include/antioch/AtomicParameter.hpp
include_HEADERS += preproc/core/include/antioch/PDF.hpp
include_HEADERS += preproc/core/include/antioch/ParameterPhy.hpp
include_HEADERS += preproc/core/include/antioch/phymath.hpp
include_HEADERS += preproc/core/include/antioch/PhyFunct.hpp
include_HEADERS += preproc/core/include/antioch/PhyHarmFunct.hpp
# preproc/chemistry
include_HEADERS += preproc/chemistry/include/antioch/SimpleChem.hpp
include_HEADERS += preproc/chemistry/include/antioch/KineticsChem.hpp
include_HEADERS += preproc/chemistry/include/antioch/ReactionChem.hpp
include_HEADERS += preproc/chemistry/include/antioch/allKineticsProcesses.hpp
include_HEADERS += preproc/chemistry/include/antioch/kineticsCrossRoad.hpp
include_HEADERS += preproc/chemistry/include/antioch/kineticsCrossRoad.hpp
# preproc/kinetics/process
include_HEADERS += preproc/kinetics/process/include/antioch/Process.hpp
include_HEADERS += preproc/kinetics/process/include/antioch/ElementaryProcess.hpp
include_HEADERS += preproc/kinetics/process/include/antioch/DuplicateProcess.hpp
include_HEADERS += preproc/kinetics/process/include/antioch/ThreeBody.hpp
include_HEADERS += preproc/kinetics/process/include/antioch/FalloffProcess.hpp
include_HEADERS += preproc/kinetics/process/include/antioch/LindemannFalloff.hpp
include_HEADERS += preproc/kinetics/process/include/antioch/TroeFalloff.hpp
# preproc/kinetics/model
include_HEADERS += preproc/kinetics/model/include/antioch/KineticsModel.hpp
include_HEADERS += preproc/kinetics/model/include/antioch/HercourtEssen.hpp
include_HEADERS += preproc/kinetics/model/include/antioch/Arrhenius.hpp
include_HEADERS += preproc/kinetics/model/include/antioch/Berthelot.hpp
include_HEADERS += preproc/kinetics/model/include/antioch/BerthelotHercourtEssen.hpp
include_HEADERS += preproc/kinetics/model/include/antioch/Kooij.hpp
include_HEADERS += preproc/kinetics/model/include/antioch/VantHoff.hpp
# preproc/constants
include_HEADERS += preproc/constants/include/antioch/physical_constants.hpp

# kinetics
include_HEADERS += kinetics/include/antioch/arrhenius_rate.h
include_HEADERS += kinetics/include/antioch/hercourtessen_rate.h
include_HEADERS += kinetics/include/antioch/berthelot_rate.h
include_HEADERS += kinetics/include/antioch/kooij_rate.h
include_HEADERS += kinetics/include/antioch/berthelothercourtessen_rate.h
include_HEADERS += kinetics/include/antioch/vanthoff_rate.h
include_HEADERS += kinetics/include/antioch/reaction.h
include_HEADERS += kinetics/include/antioch/reaction_enum.h
include_HEADERS += kinetics/include/antioch/reaction_set.h
include_HEADERS += kinetics/include/antioch/kinetics_evaluator.h

# parsing
include_HEADERS += parsing/include/antioch/tinyxml2.h
include_HEADERS += parsing/include/antioch/read_reaction_set_data_xml.h
include_HEADERS += parsing/include/antioch/blottner_parsing.h
include_HEADERS += parsing/include/antioch/sutherland_parsing.h
include_HEADERS += parsing/include/antioch/cea_thermo_ascii_parsing.h
include_HEADERS += parsing/include/antioch/species_ascii_parsing.h

# thermo
include_HEADERS += thermo/include/antioch/cea_curve_fit.h
include_HEADERS += thermo/include/antioch/cea_thermo.h
include_HEADERS += thermo/include/antioch/stat_mech_thermo.h

# transport
include_HEADERS += transport/include/antioch/blottner_viscosity.h
include_HEADERS += transport/include/antioch/sutherland_viscosity.h
include_HEADERS += transport/include/antioch/mixture_viscosity.h
include_HEADERS += transport/include/antioch/wilke_mixture.h
include_HEADERS += transport/include/antioch/eucken_thermal_conductivity.h
include_HEADERS += transport/include/antioch/wilke_evaluator.h
include_HEADERS += transport/include/antioch/constant_lewis_diffusivity.h

# utilities
include_HEADERS += utilities/include/antioch/antioch_exceptions.h
include_HEADERS += utilities/include/antioch/antioch_asserts.h
include_HEADERS += utilities/include/antioch/eigen_utils.h
include_HEADERS += utilities/include/antioch/input_utils.h
include_HEADERS += utilities/include/antioch/metaprogramming.h
include_HEADERS += utilities/include/antioch/metaprogramming_decl.h
include_HEADERS += utilities/include/antioch/metaphysicl_utils.h
include_HEADERS += utilities/include/antioch/physical_constants.h
include_HEADERS += utilities/include/antioch/string_utils.h
include_HEADERS += utilities/include/antioch/valarray_utils.h
include_HEADERS += utilities/include/antioch/vector_utils.h

# Needs to be builddir since this is generated by configure
include_HEADERS += $(top_builddir)/src/utilities/include/antioch/antioch_version.h

# Version app
antioch_version_SOURCES = apps/version.C
antioch_version_LDADD = libantioch.la

#--------------------------------------
#Local Directories to include for build
#--------------------------------------
AM_CPPFLAGS  =
AM_CPPFLAGS += -I$(top_srcdir)/src/core/include
AM_CPPFLAGS += -I$(top_srcdir)/src/kinetics/include
AM_CPPFLAGS += -I$(top_srcdir)/src/parsing/include
AM_CPPFLAGS += -I$(top_srcdir)/src/thermo/include
AM_CPPFLAGS += -I$(top_srcdir)/src/transport/include
AM_CPPFLAGS += -I$(top_srcdir)/src/preproc/core/include
AM_CPPFLAGS += -I$(top_srcdir)/src/preproc/chemistry/include
AM_CPPFLAGS += -I$(top_srcdir)/src/preproc/kinetics/model/include
AM_CPPFLAGS += -I$(top_srcdir)/src/preproc/kinetics/process/include
AM_CPPFLAGS += -I$(top_srcdir)/src/preproc/constants/include
AM_CPPFLAGS += -I$(top_srcdir)/src/utilities/include
AM_CPPFLAGS += -I$(top_builddir)/src/utilities/include #antioch_version.h
AM_CPPFLAGS +=  $(antioch_optional_INCLUDES)

#--------------
# Code coverage
#--------------
if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

#---------------------------------
# Embedded license header support
#---------------------------------

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE \
        $(top_srcdir)/src/core/include/antioch/*.h $(top_srcdir)/src/core/src/*.C \
        $(top_srcdir)/src/kinetics/include/antioch/*.h $(top_srcdir)/src/kinetics/src/*.C \
        $(top_srcdir)/src/thermo/include/antioch/*.h $(top_srcdir)/src/thermo/src/*.C \
        $(top_srcdir)/src/preproc/core/include/antioch/*.hpp $(top_srcdir)/src/preproc/core/src/*.cpp \
        $(top_srcdir)/src/preproc/chemistry/include/antioch/*.hpp $(top_srcdir)/src/preproc/chemistry/src/*.cpp \
        $(top_srcdir)/src/preproc/kinetics/process/include/antioch/*.hpp $(top_srcdir)/src/preproc/kinetics/process/src/*.cpp \
        $(top_srcdir)/src/preproc/kinetics/model/include/antioch/*.hpp $(top_srcdir)/src/preproc/kinetics/model/src/*.cpp \
        $(top_srcdir)/src/transport/include/antioch/*.h \
	$(top_srcdir)/src/utilities/include/antioch/*.h \
        $(top_srcdir)/src/utilities/src/*.C \
        $(top_srcdir)/src/apps/*.C
	echo 'updated source license headers' >$@

EXTRA_DIST = common

