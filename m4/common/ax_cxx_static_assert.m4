# ============================================================================
#  http://www.gnu.org/software/autoconf-archive/ax_cxx_auto_this.html
# ============================================================================
#
# SYNOPSIS
#
#   AX_CXX_STATIC_ASSERT([mandatory|optional])
#
# DESCRIPTION
#
#   This macro is based on AX_CXX_COMPILE_STDCXX_11, which can
#   usefully precede it to configure CXXFLAGS for C++11 support.
#
#   Check for coverage in the compiler for the presence of
#   static_assert.
#
#   The argument, if specified 'mandatory' or if left unspecified,
#   indicates that support for static_assert is required and that the
#   macro should error out if no mode with that support is found.
#   If specified 'optional', then configuration proceeds regardless,
#   after defining HAVE_CXX_STAIC_ASSERT if and only if a supporting mode
#   is found.
#
# LICENSE
#
#   Copyright (c) 2008 Benjamin Kosnik <bkoz@redhat.com>
#   Copyright (c) 2012 Zack Weinberg <zackw@panix.com>
#   Copyright (c) 2013 Roy Stogner <roystgnr@ices.utexas.edu>
#   Copyright (c) 2015 Paul T. Bauman <pbauman@buffalo.edu>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

AC_DEFUN([AX_CXX_STATIC_ASSERT], [dnl
  m4_if([$1], [], [ax_cxx_static_assert_required=true],
        [$1], [mandatory], [ax_cxx_static_assert_required=true],
        [$1], [optional], [ax_cxx_static_assert_required=false],
        [m4_fatal([invalid second argument `$1' to AX_CXX_STATIC_ASSERT])])dnl

  AC_LANG_PUSH([C++])dnl
  ac_success=no

  AC_CACHE_CHECK([for C++11 static_assert support],
                 ax_cv_cxx_static_assert,
                 [AC_COMPILE_IFELSE([AC_LANG_SOURCE([static_assert( 1==1, "Testing static_assert!");])],
                 [ax_cv_cxx_static_assert=yes],
                 [ax_cv_cxx_static_assert=no])])

  if test x$ax_cv_cxx_static_assert = xyes; then
    ac_success=yes
  fi

  AC_LANG_POP([C++])

  if test x$ax_cxx_static_assert_required = xtrue; then
    if test x$ac_success = xno; then
      AC_MSG_ERROR([*** A compiler supporting static_assert is required.])
    fi
  else
    if test x$ac_success = xno; then
      HAVE_CXX_STATIC_ASSERT=0
      AC_MSG_NOTICE([No compiler supporting static_assert was found])
    else
      HAVE_CXX_STATIC_ASSERT=1
      AC_DEFINE(HAVE_CXX_STATIC_ASSERT,1,
                [define if the compiler supports static_assert])
    fi

    AC_SUBST(HAVE_CXX_STATIC_ASSERT)
  fi
])
