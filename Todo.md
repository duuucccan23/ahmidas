# Introduction #

Use the standard GTD conventions here. No Someday/Maybe items on this list. Only actionable and well defined items.

  * Generate test propagators corresponding to the 8^4 test gauge configuration for use with unit tests.

  * Refactorize functions from some of the longer source files.

  * Define and implement a proper naming convention for the source files. See [NamingConventions](NamingConventions.md).

  * Include proper licenses and comments on all files.

  * Decide about a documentation system (Doxygen?)

  * Properly organize all the libraries (this is a persistent action item, but Mac OS X seems to be very picky about this, so it should serve as a good test.

  * Create a separate shift and multiply function to be used in the staple and square paths.

  * Create a function to create a random gauge field.

  * Implement proper ILDG headers.

  * Gaussian random number generator.

  * Quick and dirty random number generator (something like Base::Random::dirty) for speed.

  * Benchmarking code.

  * Documentation.

  * Rename QCD::reducedTensor to Dirac::Matrix (which is more intuitive).

  * Overload the "<<" operator for output streams in such a way that only one thread does the output.
