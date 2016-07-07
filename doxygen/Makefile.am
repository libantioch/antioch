AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = antioch.sty antioch.dox antioch.page txt_common fig_common

DISTCLEANFILES  = antioch.dox
DISTCLEANFILES += Makefile

dist-hook:
	rm -rf `find $(distdir)/ -name .svn`
	rm -rf `find $(distdir)/ -name .deps`

