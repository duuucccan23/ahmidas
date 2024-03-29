/*
   This is meant to prepare the source for the sequential propagator of PS mesons.
*/

//ahmidas
#include <L0/Ahmidas.h>

// C++ complex library
#include <complex>

// C++ string library
#include <cstring>

// this we need for storing the input we read from an input file
#include <map>
#include <vector>

// this we need for C++ style output
#include <iomanip>
#include <iostream>


/* *** ahmidas interfaces *** */

// representation of Dirac gamma matrices
#include <L0/Dirac/Gamma.h>

// lattice structure
#include <L0/Core/Field.h>

// propagator structure
#include <L0/Core/Propagator.h>
// // underlying Tensor structure
// #include <L0/QCD/Tensor.h>

// IO interface
#include <L1/Tool/IO.h>

// input file reader interface
#include <L2/Input/FileReader.h>



int main(int argc, char **argv)
{
  Ahmidas my_ahmidas(&argc, &argv);

  // lattice size, to be read from input file (and therefore initialized to zero)
  size_t L(0);
  size_t T(0);


  /* those are the containers to be filled by the input file reader */

  // this one contains floating point variables like kappa, mu, ...
  std::map< std::string, double > floats;
  // this one contains the file names (themselves stored in containers) of propagators (or sources, gauge fields, ...)
  std::vector< std::vector< std::string > > files;




  /* ****************************************** */
  /* ****** reading the input file ************ */
  /* ****************************************** */

  // create input file reader, the name of the input file has to be passed as a parameter
  Input::FileReader reader("../example/prova_input.xml");

  // get input parameters
  // note: this is how to invoke a member function of an object in C++:
  // <name of object>.<name of function>(<parameter list>)
  reader.initializeParameters(L, T, files, floats);

  /* ****************************************** */

  // this is needed if we want to have the output (i.e. to the standard output) done by only
  // one process in the parallel version
  Base::Weave weave(L, T);



  // that's how writing to the standard output works in C++
  // note: the "<<" operator also works for most of the ahmidas objects like SU3::Spinor or QCD::Tensor
  if(weave.isRoot())
    std::cout << "\nLattice size: " << L << "x" << L << "x" << L << "x" << T << std::endl;

  // that's how one can access the values in the map "floats"
  double kappa = floats["kappa"];
  double mu    = floats["mu"];

  if(weave.isRoot())
    std::cout << "kappa = " << kappa << ", mu = " << mu << std::endl;

  size_t t_src = size_t(floats["timesliceSource"]); 
  size_t t_srcsnk = size_t(floats["timeseparationSourceSink"]); 
  size_t t_insertion = (t_srcsnk + t_src) % T;

  if(weave.isRoot())
    std::cout << "\nThe following files are going to be read:" << std::endl;

  // there should only be one container in files, which can be accessed by files[0]
  // (similar to accessing an object in a C array)
  for (size_t fileIndex=0; fileIndex<files[0].size(); fileIndex++)
  {
    std::cout << (files[0])[fileIndex] << std::endl;
  }

  // create Propagator structure
  // note: we use dynamical memory allocation (indicated by the C++ keyword "new")
  // thus u_propagator actually is a pointer to the dynamically created object
  // the object itself (*u_propagator) can be accessed via the dereferencing operator "*"

  Core::StochasticPropagator<4> stoc_propagator(L,T);

  /* ****************************************** */
  /* ****** reading the propagator ************ */
  /* ****************************************** */

  Tool::IO::load(&stoc_propagator, files[0], Tool::IO::fileSCIDAC);
  // notes:
  // first parameter has to be pointer to Core::Propagator object
  // last parameter is precision, this can be omitted if lime files have a proper header

  /* ****************************************** */

  // here we declare a gamma-matrix
  Dirac::Gamma<5> gamma5;

  Core::StochasticPropagator<4> sequential_source(stoc_propagator);
  // we are only interested in one timeslice
  sequential_source.select_timeslice(t_insertion);

  sequential_source.rightMultiply(gamma5); 

  Tool::IO::save(&sequential_source, files[1], Tool::IO::fileSCIDAC);

  // leave main function
  return EXIT_SUCCESS;
}
