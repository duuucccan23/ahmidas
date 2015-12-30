# Setting up #

  1. The source code for Ahmidas is stored in a Google code SVN repository. So the first requirement is a working SVN-client. Tarballs may be an option at a later point, but are not made now.

## SVN ##
Getting a client:
  1. Unix/Linux: GUI variants exist in many forms, command line interfaces are usually installed by default, but always available via package managers such as apt. Note that there is a distinction between a client and a server, and only a client is needed.
  1. Windows: A very good client is [tortoisesvn](http://tortoisesvn.tigris.org/).
  1. Mac OS X: Subversion clients can be installed through various means, presumably one which keeps itself up to date is most desirable. So projects such as [Fink](http://www.finkproject.org/) may be a good choice.

Checkout instructions are available on the project website under Source.

## C++ Compiler ##
The entire code is in C++, so a C++ compiler is necessary, we generally use the GCC here, usually the latest versions.

## CMake ##
  1. The build system is [CMake](http://www.cmake.org), which is more and more available on systems or through systems such as apt and Fink.

## openMPI ##
In order to try parallel programs on a personal computer before running them on a true parallel machine, openmpi can be extremely helpful. This is generally not a preinstalled package on any system, and needs to be installed. Many linux distributions provide the openmpi-dev, openmpi-common and openmpi-bin packages through package managers such as apt, otherwise see: [open-mpi](http://www.open-mpi.org/).

Compilation is then done most easily using the MPI compiler wrapper (usually mpic++).

# Compilation #
  1. We suggest to build the entire code inside the `bin` directory in order to cleanly separate the make and object files. Therefore: run `cmake ..` from within the bin directory, or run `cmake $PATH_TO_MAIN_AHMIDAS_DIRECTORY` from the directory where you want to install everything.
  1. Then `make all` to compile all the tests and all executables.
  1. Preferably run `make test` or the more verbose `ctest -V` to perform all tests.
  1. Note that subdirectories in the source tree will contain their own CMakeList.txt files. These are used while the tree is traversed by CMake, but are not themselves suitable targets for starting a build. The top-level file, located one level _above_ the 'src' subdirectory, sets some basic parameters and starting from a lower level will usually have you encounter compilation errors.

# Documentation #
Finally, the documentation in the `doc` directory is in LaTeX format, so a working TeX distribution is necessary to view the full documentation.