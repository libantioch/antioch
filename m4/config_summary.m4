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
echo SVN revision number........... : $BUILD_VERSION
echo
echo Optional Packages for Testing:
echo '  'Eigen....................... : $enableeigen
if test "x$HAVE_METAPHYSICL" = "x1"; then
  echo '  'MetaPhysicL................. : yes
  echo '  'METAPHYSICL_CPPFLAGS ....... : $METAPHYSICL_CPPFLAGS
  echo '  'METAPHYSICL_LDFLAGS ........ : $METAPHYSICL_LDFLAGS
  echo '  'METAPHYSICL_LIBS ........... : $METAPHYSICL_LIBS
else
  echo '  'MetaPhysicL................. : no
fi
echo
echo '-------------------------------------------------------------------------------'

echo
echo Configure complete, now type \'make\' and then \'make install\'.
echo

])
