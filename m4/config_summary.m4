# SYNOPSIS
#
#   Summarizes configuration settings.
#
#   AX_SUMMARIZE_CONFIG([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Outputs a summary of relevant configuration settings.
#
# LAST MODIFICATION
#
#   2010-03-24
#

AC_DEFUN([AX_SUMMARIZE_CONFIG],
[

echo
echo '----------------------------------- SUMMARY -----------------------------------'
echo
echo Package version............... : $PACKAGE-$VERSION
echo
echo C++ compiler.................. : $CXX
echo C++ compiler flags............ : $CXXFLAGS
echo Install dir................... : $prefix 
echo Build user.................... : $USER
echo Build host.................... : $BUILD_HOST
echo Configure date................ : $BUILD_DATE
echo Build architecture............ : $BUILD_ARCH
echo Revision id................... : $BUILD_VERSION
echo
echo Testing Options:
echo '  'Number of tuples............ : $n_tuples
if test "x$HAVE_EIGEN" = "x1"; then
  echo '  'Eigen....................... : yes
else
  echo '  'Eigen....................... : no
fi
if test "x$HAVE_METAPHYSICL" = "x1"; then
  echo '  'MetaPhysicL................. : yes
  echo '    'METAPHYSICL_CPPFLAGS...... : $METAPHYSICL_CPPFLAGS
  echo '    'METAPHYSICL_LDFLAGS....... : $METAPHYSICL_LDFLAGS
  echo '    'METAPHYSICL_LIBS.......... : $METAPHYSICL_LIBS
else
  echo '  'MetaPhysicL................. : no
fi
if test "x$HAVE_VEXCL" = "x1"; then
  echo '  'VexCL....................... : yes
  echo '    'VEXCL_CPPFLAGS............ : $VEXCL_CPPFLAGS
  echo '    'VEXCL_LDFLAGS............. : $VEXCL_LDFLAGS
  echo '    'VEXCL_LIBS................ : $VEXCL_LIBS
else
  echo '  'VexCL....................... : no
fi
if test "x$HAVE_VIENNACL" = "x1"; then
  echo '  'ViennaCL.................... : yes
  echo '    'VIENNACL_CPPFLAGS......... : $VIENNACL_CPPFLAGS
  echo '    'VIENNACL_LDFLAGS.......... : $VIENNACL_LDFLAGS
  echo '    'VIENNACL_LIBS............. : $VIENNACL_LIBS
else
  echo '  'ViennaCL.................... : no
fi
if test "x$HAVE_GRVY" = "x1"; then
  echo '  'GRVY........................ : yes
  echo '    'GRVY_CFLAGS............... : $GRVY_CFLAGS
  echo '    'GRVY_LIBS................. : $GRVY_LIBS
else
  echo '  'GRVY........................ : no
fi
if test "x$HAVE_GSL" = "x1"; then
  echo '  'GSL......................... : yes
  echo '    'GSL_CPPFLAGS.............. : $GSL_CPPFLAGS
  echo '    'GSL_LDFLAGS............... : $GSL_LDFLAGS
  echo '    'GSL_LIBS.................. : $GSL_LIBS
else
  echo '  'GLS......................... : no
fi
echo
echo '-------------------------------------------------------------------------------'

echo
echo Configure complete, now type \'make\' and then \'make install\'.
echo

])
