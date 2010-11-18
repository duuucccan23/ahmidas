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
  Ahmidas my_ahmidas(&argc, &argv);

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

#ifdef UltraStocCase
  //This have to be done because ultra-stochastic source are not flavor indipendent
  size_t s0=floats["S0_flav"];
  prop1.rotateToPhysicalBasis(!s0);
#endif
  prop2.dagger();
  
  Core::Correlator<Dirac::Matrix> check(prop1 * prop2);

  check.sumOverSpatialVolume();
  
  if(weave.isRoot())
    {
      std::ofstream fout("correlators.dat");
      fout<<check;
    }

  return EXIT_SUCCESS;
}
