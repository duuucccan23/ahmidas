// $Id$
/*
   This is meant to be a tutorial how to do a contraction.
   As an example, a pion contraction is performed here, including all necessary functionality,
   such as reading an input file, loading a propagator and output of the result.
*/

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
  Input::FileReader reader("../example/pion_contraction_input.xml");

  // get input parameters
  // note: this is how to invoke a member function of an object in C++:
  // <name of object>.<name of function>(<parameter list>)
  reader.initializeParameters(L, T, files, floats);


  /* ****************************************** */



  // that's how writing to the standard output works in C++
  // note: the "<<" operator also works for most of the ahmidas objects like SU3::Spinor or QCD::Tensor
  std::cout << "\nLattice size: " << L << "x" << L << "x" << L << "x" << T << std::endl;

  // that's how one can access the values in the map "floats"
  double kappa = floats["kappa"];
  double mu    = floats["mu"];

  std::cout << "kappa = " << kappa << ", mu = " << mu << std::endl;

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

  Core::Propagator *u_propagator = new Core::Propagator(L,T);




  /* ****************************************** */
  /* ****** reading the propagator ************ */
  /* ****************************************** */


  Tool::IO::load(u_propagator, files[0], Tool::IO::fileSCIDAC, 64);
  // notes:
  // first parameter has to be pointer to Core::Propagator object
  // last parameter is precision, this can be omitted if lime files have a proper header


  /* ****************************************** */



  // Core::Correlator is the structure representing a correlation function
  // note: this is just a declaration of a pointer, no object is created yet;
  Core::Correlator *pion_twopoint;





  /* ****************************************** */
  /* ****** CONTRACTION *********************** */
  /* ****************************************** */


  /* We want to compute the (zero-) momentum projected two point function of a pion
    (interpolating field is u*gamma5*d_bar)
    C_2 (t) = sum_x'{Tr[gamma5*S_u (x',t+t_0; x_0,t_0)*gamma5*S_d (x_0,t_0 x',t+t_0)]},
            = sum_x'{Tr[gamma5*S_u (x',t+t_0; x_0,t_0)*gamma5*(gamma_5*S^dagger_u(x',t+t_0; x_0,t_0)*gamma5)]}
            = sum_x'{Tr[S_u (x',t+t_0; x_0,t_0)*S^dagger_u(x',t+t_0; x_0,t_0)]}
     where S is a quark propagator, x' is the sink 3d-position
     and x_0, t_0 are the fixed source position and timeslice.
     The trace is over spin and colour, the "dagger" acts on spin and color indices only.
     Note that although we can leave out the gamma5 we are actually performing the gamma-multiplication
     for educational reasons.
  */


  /* First of all we calculate gamma_5*S^dagger_u*gamma5 */

  // for that we need a copy of the up quark Propagator
  Core::Propagator *d_bar_propagator = new Core::Propagator(*u_propagator);
  // the following routine does exacly what we want
  // note: this is how you invoke the member function of a C++ object from the pointer to that object,
  // instead of using a "." you use "->"
  d_bar_propagator->revert();


  // here we declare a gamma-matrix
  Dirac::Gamma<5> gamma5;

  // multiply u_propagator by gamma5
  // note: since u_propagator is a pointer we have to access
  // the corresponding object using the dereferencing operator "*"
  (*u_propagator) *= gamma5;
  // same for d_bar_propagator
  (*d_bar_propagator) *= gamma5;


  // create the Correlator object
  // note: multiplying the two Propagator objects already performs the colour trace,
  // the full Dirac structure is kept for a reason
  pion_twopoint = new Core::Correlator(L, T, (*u_propagator)*(*d_bar_propagator));


  // this does the zero-momentum projection
  // note: the same function with an argument projects to any momentum
  pion_twopoint->sumOverSpatialVolume();


  // everything allocated dynamically should be deleted when it's not used anymore
  // (actually it's not necessary, but it's a good C++ programming standard and
  //  it makes your programm safe against unwanted memory access)
  delete u_propagator;
  delete d_bar_propagator;


  // Essentially everything is computed now, we only have to perform the Dirac trace.
  // We can access the full Dirac structure of the momentum projected two point function
  // of a timeslice like an element of a C array using the access operator "[<timeslice>]".
  // The member function trace() returns the trace of that Dirac structure (this is just a 4x4 matrix).
  std::cout << "\nPion two point function:" << std::endl;
  for (size_t t=0; t<T; t++)
  {
    // this is the way formatted output works in C++
    std::cout.width(3);
    std::cout << t << "  ";
    std::cout.width(20);
    // since the Correlator stores complex numbers, we can access the real and imaginary parts using
    // the complex class member functions real() and imag()
    std::cout << std::fixed << std::setprecision(10) << std::showpos << ((*pion_twopoint)[t]).trace().real() << "  ";
    std::cout.width(20);
    std::cout << std::fixed << std::setprecision(10) << std::showpos << ((*pion_twopoint)[t]).trace().imag() << std::endl;
  }


  /* ****************************************** */



  // now we probably want to know whether we have obtained the correct result,
  // we therefore print the values of another well-tested contraction routine
  std::cout <<  "\nOutput of reliable contraction routine:\n" << std::endl;
  std::cout <<  " 0  +0.5412652273    +0" << std::endl;
  std::cout <<  " 1  +0.01456410538   +0" << std::endl;
  std::cout <<  " 2  +0.001637160312  +0" << std::endl;
  std::cout <<  " 3  +0.01443586407   +0\n" << std::endl;


  // Correlator no longer needed
  delete pion_twopoint;

  // leave main function
  return EXIT_SUCCESS;
}