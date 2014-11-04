Atom
===

[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT) [![Build Status](https://travis-ci.org/kartikkumar/atom.svg?branch=master)](https://travis-ci.org/kartikkumar/atom)

`Accurate Transfer Orbit Model (Atom)` is a C++ template library containing the `Atom` solver. The `Atom` solver is an analog of the Lambert solver using the SGP4/SDP4 propagation models.

Details of the solver can be found in [Kumar, et al. (2014)](#temp)

Features
------

  - Header-only
  - Atom solver function with fully-configurable optional parameters
  - Cartesian-to-TLE conversion function
  - Full suite of tests

Requirements
------

To install this project, please ensure that you have installed the following (install guides are provided on the respective websites):

  - [Git](http://git-scm.com)
  - A C++ compiler, e.g., [GCC](https://gcc.gnu.org/), [clang](http://clang.llvm.org/), [MinGW](http://www.mingw.org/)
  - [CMake](http://www.cmake.org)
  - [Doxygen](http://www.doxygen.org "Doxygen homepage") (optional)

In addition, `Atom` depends on the following libraries:

  - [SML](https://www.github.com/kartikkumar/sml)
  - [SAM](https://www.github.com/kartikkumar/sam)
  - [GSL](http://www.gnu.org/software/gsl)
  - [SGP4](https://www.github.com/kartikkumar/sgp4deorbit)         
  - [Catch](https://www.github.com/philsquared/Catch)
  - [Eigen](http://eigen.tuxfamily.org/) (optional)

These dependencies will be downloaded and configured automatically if they are not already present locally (requires an internet connection).

Installation
------

Run the following commands to download, build, and install this project.

    git clone https://www.github.com/kartikkumar/atom
    cd atom
    git submodule init && git submodule update
    mkdir build && cd build
    cmake ..
    cmake --build .

To install the header files and libraries, run the following from within the `build` directory:

    make install

Note that dependencies are installed by fetching them online, in case they cannot be detected on your local system. If the build process fails, check the error log given. Typically, building fails due to timeout. Simply run the `cmake --build .` command once more.

Build options
-------------

You can pass the follow command-line options when running `CMake`:

  - `-DBUILD_DOCS=[on|off (default)]`: build the [Doxygen](http://www.doxygen.org "Doxygen homepage") documentation ([LaTeX](http://www.latex-project.org/) must be installed with `amsmath` package)
  - `-DBUILD_TESTS`=[on|off (default)]: build tests (execute tests from build-directory using `make test`)
  - `-DBUILD_WITH_EIGEN=[on|off (default)]`: build tests using [Eigen](http://eigen.tuxfamily.org/)
  - `-DBUILD_SHARED_LIBS=[on|off (default)]`: build shared libraries instead of static
  - `-DCMAKE_INSTALL_PREFIX`: set path prefix for install script (`make install`); if not set, defaults to usual locations
  - `-DFORCE_DEPENDENCIES_BUILD=[on|off (default)]`: force local build of dependencies, instead of first searching system-wide using `find_package()`
  - `-DMYLIB_PATH[=build_dir/lib (default]`: set library path
  - `-DMYBIN_PATH[=build_dir/bin (default]`: set binary path
  - `-DMYTEST_PATH[=build_dir/tests (default]`: set tests path

Contributing
------------

Once you've made your great commits:

1. [Fork](https://github.com/kartikkumar/atom/fork) `Atom`
2. Create a topic branch - `git checkout -b my_branch`
3. Push to your branch - `git push origin my_branch`
4. Create a [Pull Request](http://help.github.com/pull-requests/) from your branch
5. That's it!

License
------

See `LICENSE.md`.

Disclaimer
------

The copyright holders are not liable for any damage(s) incurred due to improper use of `Atom`.

Contact
------

Shoot an [email](mailto:me@kartikkumar.com?subject=Atom) if you have any questions.

TODO
------

  - Improve convergence radius of `Atom` solver (which includes figuring out the best scaling for the residual function)
  - Extend test suite
  - Figure out better way (avoiding code duplication) to build `STL`-based and `Eigen`-based tests in the same build tree
  - Add version detection in `CMake` module so that find_package respects minimum version required.
  - Find a way to provide an option to clean installation. 
