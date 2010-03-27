// $Id$

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

 
  /* ****************************************** */
  /* ****** reading the input file ************ */
  /* ****************************************** */
  
  
  /*
      function declarations:
      
      FileReader::FileReader(std::string const file);
      
      void FileReader::initializeParameters(size_t &L, size_t &T,
                                            std::vector< std::vector< std::string > > &filenames,
                                            std::map< std::string, double > &floats) const;
  */


  /* ****************************************** */




  /* ****************************************** */
  /* ****** reading the propagator ************ */
  /* ****************************************** */


  /*
      function declaration:
      
      void Tool::IO::load(Core::Propagator *propagator, std::vector< std::string > const &filenames,
                Tool::IO::filetype type, size_t const precision);
  */


  /* ****************************************** */




  /* ****************************************** */
  /* ****** CONTRACTION *********************** */
  /* ****************************************** */


  /*
      function declarations:
      
      Propagator &Propagator::revert();
      
      template< size_t Index >
      void Propagator::operator*=(Dirac::Gamma< Index > const &gamma);
      
      Correlator::Correlator(size_t const L_, size_t const T_, Field < QCD::reducedTensor > *d_data);
      
      Core::Field< QCD::reducedTensor > *Propagator::operator*(Propagator const &other) const;
      
      void  Correlator::sumOverSpatialVolume();

  */ 



  /* ****************************************** */


//   std::cout <<  "\nOutput of reliable contraction routine:" << std::endl;
//   std::cout <<  " 0  +0.5412652273    +0" << std::endl;
//   std::cout <<  " 1  +0.01456410538   +0" << std::endl;
//   std::cout <<  " 2  +0.001637160312  +0" << std::endl;
//   std::cout <<  " 3  +0.01443586407   +0\n" << std::endl;


  return 0;
}
