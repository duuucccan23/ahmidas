//This program produce the sequential source for meson correlation
//functions. It takes as a source the S0 already rotated in physical
//basis, insert the appropriate operator (for the moment, only gamma5,
//but in future it will be more general) and rotate everything back to
//the twisted basis

//After using this program, use an inverter over it, and use
//"finalize_meson_sequential_propagator" over its output to obtain the
//correct S1.

//Note: this code is approriate for the charged S1, not for the neutral.

//Explanation of what is done. Let's take the S1 of of a '0' S0

// We want to compute D_up^-1*g5*D_dw^-1 , to obtain it we do: 
// 0.25*(1+ig5)*(D_+)^-1*(1+ig5)*g5*(1-ig5)*(D_-)^-1*(1-ig5)
//
// this is the s0 taken as input:
//                                  (1-ig5)*(D_-)^-1*(1-ig5)/2
// this is the output of this prog:
//                       (1+ig5)*g5*(1-ig5)*(D_-)^-1*(1-ig5)*(1+ig5)/4

#include <cstring>
#include <vector>
#include <map>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L0/QCD/Gauge.h>
#include <L0/Core/Propagator.h>
#include <L2/Contract/Baryon.h>
#include <L1/Tool/IO.h>
#include <L1/Smear.h>
#include <L1/Smear/APE.h>
#include <L1/Smear/Jacobi.h>
#include <L2/Input/FileReader.h>


int main(int narg, char **arg)
{
  size_t L=0;
  size_t T=0;

  if(narg<2)
    {
      std::cerr<<"Use: "<<arg[0]<<" input_file"<<std::endl;
      exit(0);
    }

  Input::FileReader reader(arg[1]);

  std::map< std::string, double > floats;
  std::vector< size_t * > positions;
  std::vector< int > operators;
  std::vector< std::vector< std::string > > files;

  reader.initializeParameters(L,T,files,floats,positions,operators);

  Base::Weave weave(L,T);
  weave.barrier();

  if(weave.isRoot()) std::cout<<"Lattice size: "<<L<<"x"<<L<<"x"<<L<<"x"<<T<<std::endl;

  bool const take_slice=bool(floats["take_slice"]);
  size_t const tslice=size_t(floats["tslice"]);
  if(take_slice) if(weave.isRoot()) std::cout<<"will take: "<<tslice<<" timeslice\n"<<std::endl;  

  size_t const s0_flav=size_t(floats["s0_flav"]);
  if(weave.isRoot()) std::cout<<"flavour of s0: "<<s0_flav<<std::endl;  

  PropagatorType prop(L,T);
  Tool::IO::load(&prop,files[0],Tool::IO::fileSCIDAC);
  if(weave.isRoot()) std::cout<<"quark propagator successfully loaded\n"<<std::endl;

  //this is the multiplication by gamma5, to be generalized
  {
    Dirac::Gamma<5> gamma5;
    prop.rightMultiply(gamma5);
    if(weave.isRoot()) std::cout<<"insertion of operator performed\n"<<std::endl;
  }

  //rotate back to physical basis, according to the s0 flavor performing the opposite rotation
  prop.rotateToPhysicalBasis(!s0_flav);
  if(weave.isRoot()) std::cout<<"rotation performed\n"<<std::endl;  

  //if it have to be specified to take a slice, take it
  if(take_slice)
    {
      prop.select_timeslice(tslice); 
      if(weave.isRoot()) std::cout<<"timeslice "<<tslice<<" selected\n"<<std::endl;
    }

  //output the sequential source
  std::vector< std::string > file_out;
  for(size_t j=0;j<files[0].size();j++) file_out.push_back(files[0][j]+".seq");
  Tool::IO::save(&prop,file_out,Tool::IO::fileSCIDAC);
  std::cout<<"sequential source generated and saved successfully\n"<<std::endl;

  return EXIT_SUCCESS;
}
