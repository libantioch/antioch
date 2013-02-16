BUILT_SOURCES   = .license.stamp

check_PROGRAMS  =
check_PROGRAMS += chem_mixture_unit
check_PROGRAMS += cea_thermo_unit
check_PROGRAMS += arrhenius_rate_unit

AM_CPPFLAGS  = 
AM_CPPFLAGS += -I$(top_srcdir)/src/core/include
AM_CPPFLAGS += -I$(top_srcdir)/src/kinetics/include
AM_CPPFLAGS += -I$(top_srcddir)/src/parsing/include
AM_CPPFLAGS += -I$(top_srcdir)/src/thermo/include
AM_CPPFLAGS += -I$(top_srcdir)/src/utilities/include
AM_CPPFLAGS += -I$(top_builddir)/src/utilities/include #antioch_version.h

AM_LDFLAGS = $(top_builddir)/src/libantioch.la

# Sources for these tests
chem_mixture_unit_SOURCES = chem_mixture_unit.C
cea_thermo_unit_SOURCES = cea_thermo_unit.C
arrhenius_rate_unit_SOURCES = arrhenius_rate_unit.C

#Define tests to actually be run
TESTS  =
TESTS += chem_mixture_unit
TESTS += cea_thermo_unit
TESTS += arrhenius_rate_unit

if CODE_COVERAGE_ENABLED
  CLEANFILES = *.gcda *.gcno
endif

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE \
	$(top_srcdir)/test/*.C
	echo 'updated source license headers' >$@

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
