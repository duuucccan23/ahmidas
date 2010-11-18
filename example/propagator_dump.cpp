#include <cstring>
#include <vector>
#include <map>
#include <complex>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>

#include <L0/Ahmidas.h>
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
  Ahmidas my_ahmidas(&narg, &arg);


  if(narg<4)
    {
      std::cerr<<"Error: use "<<arg[0]<<" L T file_base NFILE"<<std::endl;
      exit(0);
    }
  
  size_t L=atoi(arg[1]);
  size_t T=atoi(arg[2]);
  size_t NF=atoi(arg[4]);
 
  std::vector<std::string> files;

  char suff[3];
  for(int i=0;i<NF;i++)
    {
      sprintf(suff,".%02d",i);
      files.push_back(arg[3]+std::string(suff));
    }
  
  Base::Weave weave(L,T);
  weave.barrier();

  if(weave.isRoot()) std::cout<<"Lattice size: "<<L<<"x"<<L<<"x"<<L<<"x"<<T<<std::endl;

  //load the propagator
  if(NF==1)
    {
      Core::StochasticPropagator<1> prop(L,T);
      Tool::IO::load(&prop,files,Tool::IO::fileSCIDAC);
      if(weave.isRoot()) std::cout<<std::endl<<"propagator successfully loaded\n"<<std::endl;
      std::cout<<prop<<std::endl;
      std::cerr<<sqrt(prop.normq())<<std::endl;
    }
  if(NF==4)
    {
      Core::StochasticPropagator<4> prop(L,T);
      Tool::IO::load(&prop,files,Tool::IO::fileSCIDAC);
      if(weave.isRoot()) std::cout<<std::endl<<"propagator successfully loaded\n"<<std::endl;
      std::cout<<prop<<std::endl;
      std::cerr<<sqrt(prop.normq())<<std::endl;
    }
  else if(NF==12)
    {
      Core::StochasticPropagator<12> prop(L,T);
      Tool::IO::load(&prop,files,Tool::IO::fileSCIDAC);
      if(weave.isRoot()) std::cout<<std::endl<<"propagator successfully loaded\n"<<std::endl;
      std::cout<<prop<<std::endl;
      std::cerr<<sqrt(prop.normq())<<std::endl;
    }
      
  return EXIT_SUCCESS;
}
