#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

#include <L0/Print.h>
#include <L0/Ahmidas.h>
#include <L0/Core/Propagator.h>
#include <L0/Core/Correlator.h>
#include <L1/Tool/IO.h>
#include <L2/Contract/Meson.h>
#include <L2/Input/FileReader.h>

int main(int argc, char **argv)
{
  // Ahmidas start(&argc, &argv);

  size_t L_tmp(0);
  size_t T_tmp(0);

  //read inputs from parameter file
  std::map< std::string, double > floats;
  std::vector< std::vector< std::string > > files;
  std::vector< size_t * > positions; // not actually used
  std::vector< int > operators;

  {
    Input::FileReader reader("./check_meson_3pts_input.xml");
    reader.initializeParameters(L_tmp, T_tmp, files, floats, positions, operators);
  }

  size_t const L = L_tmp;
  size_t const T = T_tmp;

  Base::Weave weave(L, T);

  {
    // write some output concerning what is going to be read
    std::ostringstream out;
    out << "Lattice size: " << L << "x" << L << "x" << L << "x" << T << std::endl;
  }

  std::vector< std::string > const &propFiles1 = files[0];
  std::vector< std::string > const &propFiles2 = files[1];

  PropagatorType prop1(L, T);
  PropagatorType prop2(L, T);

  Tool::IO::load(&prop1, propFiles1, Tool::IO::fileSCIDAC);
  if(weave.isRoot()) std::cout<<"First propagator successfully loaded"<<std::endl;
  Tool::IO::load(&prop2, propFiles2, Tool::IO::fileSCIDAC);
  if(weave.isRoot()) std::cout<<"Source successfully loaded"<<std::endl;

  size_t localIndex;
  size_t nf=files[0].size();

  Dirac::Gamma<5> gamma5;
  prop1.rightMultiply(gamma5);

  prop2.dagger();
  //std::cout<<(prop1*prop2)<<std::endl;
  std::cout<<"Fuffa"<<std::endl;
  for(size_t idx_T = 0; idx_T < T; idx_T++)
    {
      std::complex<double> sum;
      sum=0;
      for(size_t idx_Z = 0; idx_Z < L; idx_Z++)
	for(size_t idx_Y = 0; idx_Y < L; idx_Y++)
	  for(size_t idx_X = 0; idx_X < L; idx_X++)
	    {
	      localIndex = weave.globalCoordToLocalIndex(idx_X, idx_Y, idx_Z, idx_T);
	      if (localIndex != weave.localVolume())
		{
		  for(size_t iD1=0;iD1<4;iD1++)
		    for(size_t iD2=0;iD2<4;iD2++)
		      for(size_t iC1=0;iC1<3;iC1++)
			sum+=prop1[localIndex][iD1*3][iD2][iC1]*prop2[localIndex][iD2*3][iD1][iC1];
		}
	    }
      std::cout<<idx_T<<" "<<sum<<std::endl;
    }
  Core::Correlator<Dirac::Matrix> check(prop1 * prop2);

  check.sumOverSpatialVolume();
  
  if(weave.isRoot())
    {
      std::ofstream fout("correlators.dat");
      fout<<check;
    }

  return EXIT_SUCCESS;
}
