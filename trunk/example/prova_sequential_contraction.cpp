
/*
   This is meant to be the contrction code for 3pts functions of PS mesons.
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

// correlation function structure
#include <L0/Core/Correlator.h>

// this contains the contractions we actually perform directly in the main function here
#include <L2/Contract/Meson.h>

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
  Input::FileReader reader("../example/prova_sequential_contraction_input.xml");

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

//   size_t t_src = size_t(floats["timesliceSource"]); 

  if(weave.isRoot())
    std::cout << "\nThe following files are going to be read:" << std::endl;

  // there should only be one container in files, which can be accessed by files[0]
  // (similar to accessing an object in a C array)
  for (size_t fileIndex=0; fileIndex<files[0].size(); fileIndex++)
  {
    std::cout << (files[0])[fileIndex] << std::endl;
  }

  Core::StochasticPropagator<4> stoc_propagator(L,T);

  Core::StochasticPropagator<4> sequential_propagator(L,T);

  Core::StochasticSource<4> stoc_source(L,T);

  /* ****************************************** */
  /* ****** reading the propagators *********** */
  /* ****************************************** */

  Tool::IO::load(&sequential_propagator, files[0], Tool::IO::fileSCIDAC, 64);

  Tool::IO::load(&stoc_propagator, files[1], Tool::IO::fileSCIDAC);

  Tool::IO::load(reinterpret_cast< Core::StochasticPropagator<4> * >(&stoc_source), files[2], Tool::IO::fileSCIDAC, 64);

  /* ****************************************** */




  /* ****************************************** */
  /* ****** CONTRACTION *********************** */
  /* ****************************************** */


  /* We want to compute the three point function of PS mesons
    (interpolating field is u*gamma5*d_bar)
    C_3 (t) = sum_x'{Tr[gamma5*S_u (x',t+t_0; x_0,t_0)*gamma0*S_d (x_0,t_0 x',t+t_0)]},
            = sum_x'{Tr[gamma5*S_u (x',t+t_0; x_0,t_0)*gamma0*(gamma5*S^dagger_u(x',t+t_0; x_0,t_0)*gamma5)]}
            = sum_x'{Tr[S_u (x',t+t_0; x_0,t_0)*gamma05*S^dagger_u(x',t+t_0; x_0,t_0)]}
     where S is a quark propagator, x' is the sink 3d-position
     and x_0, t_0 are the source position and timeslice.
     The trace is over spin and colour, the "dagger" acts on spin and color indices only.
  */
 // here we declare a gamma-matrix
  Dirac::Gamma<5> gamma5;
//   Dirac::Gamma<45> gamma05;

  /* First of all we calculate gamma_5*S^dagger_u*gamma5 */

  // for that we need a copy of the up quark Propagator
  Core::StochasticPropagator<4> stoc_bar_propagator(stoc_propagator);
  // the following routine does exacly what we want
  // note: this is how you invoke the member function of a C++ object from the pointer to that object
  stoc_bar_propagator.dagger();
  stoc_bar_propagator *= gamma5;

  // multiply sequential_propagator by gamma05
//  sequential_propagator *= gamma05;

  stoc_source.dagger();
  stoc_source.shift(Base::idx_T, Base::dir_UP, 2);

  sequential_propagator.rightMultiply(gamma5);

  // create the Correlator object
  // note: multiplying the two Propagator objects already performs the colour trace,
  // the full Dirac structure is kept for a reason
  Core::Correlator< Dirac::Matrix > ps_3point((dynamic_cast< Core::Propagator & >(sequential_propagator)) * (dynamic_cast< Core::Propagator & >(stoc_source)));

  // this does the zero-momentum projection
  // note: the same function with an argument projects to any momentum
  ps_3point.sumOverSpatialVolume();




  // Essentially everything is computed now, we only have to perform the Dirac trace.
  // We can access the full Dirac structure of the momentum projected two point function
  // of a timeslice like an element of a C array using the access operator "[<timeslice>]".
  // The member function trace() returns the trace of that Dirac structure (this is just a 4x4 matrix).
  if(weave.isRoot())
    std::cout << "\nPS three point function:" << std::endl;

  std::cout << ps_3point;


  // leave main function
  return EXIT_SUCCESS;
}
