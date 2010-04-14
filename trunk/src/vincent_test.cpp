// $Id: main.cpp 393 2010-03-28 16:32:36Z dinter@ifh.de $

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

#include <time.h>

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

// input APE
#include <L1/Smear/APE.h>



int main(int argc, char **argv)
{
  size_t L = 0;
  size_t T = 0;



  /* ****************************************** */
  /* ****** reading the input file ************ */
  /* ****************************************** */
  
  Input::FileReader reader("../test/input.xml");

  std::map< std::string, double > floats;
  std::vector< size_t * > positions;
  std::map< std::string, int > operators;
  std::vector< std::vector< std::string > > files;

  reader.initializeParameters(L, T, files, floats, positions, operators);

  std::cout << "Lattice size: " << L << "x" << L << "x" << L << "x" << T << std::endl;

  double kappa = floats["kappa"];
  double mu = floats["mu"];

  std::cout << "kappa = " << kappa << ", mu = " << mu << std::endl;




  

  std::vector< std::string > const &propfilesU(files[0]);
  std::vector< std::string > const &propfilesD(files[1]);
  std::string  const &confname("../test/conf.9999");

  
   std::cout << confname << std::endl;
    
   Core::Field<QCD::Gauge>  gauge(L, T);
#ifdef __MPI_ARCH__
  if (myid == 0)
#endif  
     std::cout << "\nStart to read the gauge configuration \n"  << std::endl;
  
  Tool::IO::load(&gauge,confname ,Tool::IO::fileILDG);
#ifdef __MPI_ARCH__
  if (myid == 0)
#endif  
     std::cout  << "Reading ended  \n"  << std::endl;
  

  // Number of propagator to read . Have to be specified in the input file.

  size_t const Nprop=3;
  


  // Create a vector which contain the vector of the propagator name . propfile[] should be defined in the input file.

  std::vector< std::string > *propfiles[Nprop];
  propfiles[0]= new  std::vector< std::string >(propfilesU);
  propfiles[1]= new  std::vector< std::string >(propfilesD);
  propfiles[2]= new  std::vector< std::string >(propfilesU);

  

/* ******************************************** */
/* ******  propagator declaration  ************ */
/* ******************************************** */



//Initialize a vector of propagators for u, d , s quarks.
  
  Core::Propagator *FlavourProp[Nprop];



 /* ****************************************** */
 /* ******  reading the propagator **********  */
 /* ****************************************** */



  std::cout << "\nStart to read propagators \n"  << std::endl;
  
  for (size_t i=0; i<Nprop; i++)
  {
     for (int f=0; f<12; f++)
     {
#ifdef __MPI_ARCH__
  if (myid == 0)
#endif
	std::cout<< (*propfiles[i])[f] << std::endl;
     }
     
     FlavourProp[i] = new Core::Propagator(L, T);
     Tool::IO::load( FlavourProp[i], *propfiles[i], Tool::IO::fileSCIDAC, 64);

#ifdef __MPI_ARCH__
  if (myid == 0)
#endif     
     std::cout << "\nquark propagator successfully loaded \n" << std::endl; 
  }

#ifdef __MPI_ARCH__
  if (myid == 0)
#endif
     std::cout << "\nReading ended ..." << Nprop << " quark propagators successfully loaded" << std::endl;


// Change boundary condition at the propagator level  (from Simon )

  size_t timeslice_source(0);
  size_t timeslice_boundary(T-1);

  for (size_t i=0; i<Nprop; i++)  FlavourProp[i]->changeBoundaryConditions_uniformToFixed(timeslice_source, timeslice_boundary);
  
  std::cout << Base::bar_PROTON  << std::endl;
 
 
 

 /* ****************************************** */
 /* ******  APE smearing  ******************** */
 /* ****************************************** */
  


double elapsedTime;
clock_t stopTime;
clock_t startTime = clock();



  double strength = 0.5; // to be read in the input file
  Smear::APE SmearedGauge = Smear::APE(strength); 
  
  size_t times = 10; // to be read in the input file
  std::cout << "APE smearing gauge field " << times << " times with strength " << strength << ".\n";

  SmearedGauge.smear(gauge, times);

  
 

stopTime = clock();
elapsedTime = (stopTime - startTime) / (CLOCKS_PER_SEC / (double) 1000.0);

std::cout << "Ape smearing done in "<< elapsedTime << "ms\n" << std::endl;

startTime = clock();


//Dirac Matrices

 Dirac::Gamma<0> gamma0;
 Dirac::Gamma<1> gamma1;
 Dirac::Gamma<2> gamma2;
 Dirac::Gamma<3> gamma3;
 Dirac::Gamma<5> gamma5;





//  Dirac::Gamma<1> Gamma_[1];
//  std::vector<Dirac::Gamma<>*> bla;



// The basic structure I want to code is the following , denoting X(x,0) Y(x,0) Z(x,0)  three propagators
// Monomial1 <-  eps^abc eps^def X^de(x,0) Gamma  Y^be(x,0) Gamma_tilde Z^cf(x,0)
// with Monomial1 is a 4x4 matrix ( all colour indices are summed)
// There is 9 "Monomials" which are of  all of the same type ( some transpotition of spinor index)
// Gamma and Gamma_tilde are arbitrary gamma matrices related to the Dirac matrices in the interpolating fields of source/sink

for (size_t i=0; i<Nprop; i++)  // flavour index
  {
     for (size_t j=0; j<Nprop; j++) // flavour index
     { 
	for (size_t k=0; k<Nprop; k++)  // flavour index
	{
//	   for (size_t  mu=0; mu<4; mu++)  // Gamma matrices index
	   //   {
	      for (size_t nu=0; nu<4; nu++)  // Gamma matrices index
		 // {
		 //	 int const l=mu; 
		 //int const m=nu; 
		 //Dirac::Gamma<l> bla;
		 //Dirac::Gamma<m> bla2;
		
		 
		 Core::Correlator C(L, T,(*FlavourProp[i])*gamma0*gamma1*gamma5*(*FlavourProp[j]));
		 

//	      }// Gamma matrices index
//	   }  // Gamma matrices index
	}// flavour index
     }// flavour index
  }// flavour index

     

stopTime = clock();
elapsedTime = (stopTime - startTime) / (CLOCKS_PER_SEC / (double) 1000.0);

std::cout << "Contraction done in "<< elapsedTime << "ms\n" << std::endl;

  
  //Gaussian smearing of the propagators.
  
  // Rotation physical base.
  
  // Random sources position

  // Contraction octet/decuplet standard 

  // write output 



//clean up


for (size_t i=0; i<Nprop; i++)
 {
   delete FlavourProp[i];
 }

//delete  &gauge;
//delete  FlavourProp;

  
  /* ****************************************** */
  /* ****** reading the input file ************ */
  /* ****************************************** */
  
  
  /*
      function declarations:
      
      FileReader::FileReader(std::string const file);
      
      void FileReader::initializeParameters(size_t &L, size_t &T,
                                            std::vector< std::vector< std::string > > &filenames,
                                            std::map< std::string, double > &floats) const;
         input file:
         ../example/pion_contraction_input.xml

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


#ifdef __MPI_ARCH__
  if (myid == 0)
#endif
  std::cout << "\nprogramm is going to exit normally now\n" << std::endl;

  return EXIT_SUCCESS;

}
