[![Build Status](https://travis-ci.org/libantioch/antioch.png)](https://travis-ci.org/libantioch/antioch)

Antioch
=======

A New Templated Implementation Of Chemistry for Hydrodynamics (Antioch) was initiated
to centralize work by some of the Antioch authors within the realm of hypersonic
aerodynamics, based on the [libMesh](https://github.com/libMesh/libmesh.git) finite
element library. In particular, although there exist C++ chemistry libraries, such
as [Cantera](http://code.google.com/p/cantera/), we had needs for both thread-safety
and high performance. Thus, Antioch was born. Antioch originally lived within
the [PECOS](http://pecos.ices.utexas.edu) center at the Institute for Computational
Engineering and Sciences ([ICES](https://www.ices.utexas.edu))
at [The University of Texas at Austin](https://www.utexas.edu).

Dependencies
============

Requirements
------------

Antioch has no required dependencies other than a reasonably modern C++ compiler.

Optional Packages
-----------------

Antioch has been designed to allow several optional packages to facilitate vectorized
evaluation of thermochemistry quantities. In particular, Antioch currently can support:

1. [VexCL](https://github.com/ddemidov/vexcl.git) - This package will enable GPU offload capabilities
2. [Eigen](http://eigen.tuxfamily.org) - Highly optimized vector and matrix types

Building Antioch
================

Antioch uses an Autotools build system, so typical GNU build commands are used.

1. ./bootstrap (generate configure script)
2. ./configure --prefix=/path/to/install (for more options, do ./configure --help)
3. make (note parallel builds are supported)
4. make check (note parallel-tests are supported)
5. make install

Discussion and Development
==========================

Pull Requests for bug fixes and new features are welcome. Mailing lists have been setup for
user questions and discussion (antioch-users@googlegroups.com) as well as development questions
and discussion (antioch-devel@googlegroups.com).
