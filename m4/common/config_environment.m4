# SYNOPSIS
#
#   Queries configuration environment.
#
#   AX_SUMMARIZE_ENV([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Queries compile environment and SVN revision for use in configure summary
#   and pre-processing macros.
#
# LAST MODIFICATION
#
#   $Id: config_environment.m4 34934 2012-11-28 18:15:04Z benkirk $
#
# COPYLEFT
#
#   Copyright (c) 2010 Karl W. Schulz <karl@ices.utexas.edu>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_SUMMARIZE_ENV],
[

AC_CANONICAL_HOST

BUILD_USER=${USER}
BUILD_ARCH=${host}
BUILD_HOST=${ac_hostname}
BUILD_DATE=`date +'%F %H:%M'`

# Determine method for querying Source code revisioning (assumes git)

AC_PATH_PROG(gitquery,git)

if (test "x${gitquery}" = "x" -o ! -d $srcdir/.git) ; then
   GIT_REVISION="external" #"cat $srcdir/dist_version"
   GIT_CHECKOUT=false
   BUILD_DEVSTATUS="External Release"
else

   #echo "git revision=`cd $srcdir ; ${gitquery} rev-parse HEAD`"
   #echo "git branch=`cd $srcdir ; ${gitquery} branch | awk '/^\\* / { print $2 }'`"

   GIT_REVISION="`cd $srcdir ; ${gitquery} rev-parse HEAD`"
   GIT_CHECKOUT=true
   BUILD_DEVSTATUS="Development Build"
fi


AC_SUBST(GIT_REVISION)
AC_SUBST(BUILD_DEVSTATUS)
AM_CONDITIONAL(GIT_CHECKOUT,test x${GIT_CHECKOUT} = xtrue )

# Query current version.

BUILD_VERSION=${GIT_REVISION}

# Versioning info - check local developer version (if checked out)

AC_DEFINE_UNQUOTED([BUILD_USER],     "${BUILD_USER}",     [The fine user who built the package])
AC_DEFINE_UNQUOTED([BUILD_ARCH],     "${BUILD_ARCH}",     [Architecture of the build host])
AC_DEFINE_UNQUOTED([BUILD_HOST],     "${BUILD_HOST}",     [Build host name])
AC_DEFINE_UNQUOTED([BUILD_VERSION],  "${BUILD_VERSION}",  [SVN revision])
AC_DEFINE_UNQUOTED([BUILD_DEVSTATUS],"${BUILD_DEVSTATUS}",[Dev/Release build])
AC_DEFINE(         [BUILD_DATE],     __DATE__ " " __TIME__, [Build date])

AC_SUBST(BUILD_USER)
AC_SUBST(BUILD_ARCH)
AC_SUBST(BUILD_HOST)
#AC_SUBST(BUILD_DATE)
AC_SUBST(BUILD_VERSION)

])
