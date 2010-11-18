//File0: configuration, needed to produce the two propagators from DD^-1
//File1: input propagator

//output propagators file names are created automatically adding ".0"
//or ".1" to filenames. These follow APE naming scheme. The 0 is the
//solution of the Dirac Operator containing: -i*mu*g5, which rotate to
//physical basis as (1+ig5)/sqrt(2). The 1 is the other.
// 0 is sometimes reffered as up in the following.
//
// Note that the Dirac Operator and the relative propagator rotate one
// the opposite than the other

#include <complex>
#include <cstring>
#include <map>
#include <vector>
#include <iomanip>
#include <iostream>
#include <L0/Ahmidas.h>
#include <L0/Base/Weave.h>
#include <L0/Core/Field.h>
#include <L0/Core/Propagator.h>
#include <L1/Tool/IO.h>
#include <L2/Input/FileReader.h>

int main(int narg,char **arg)
{
  size_t L,T;
  std::string suff[2]={".0",".1"}; //suffix for output files (up,down)

  if(narg<2)
    {
      std::cerr<<"Use: "<<arg[0]<<" parameters_file"<<std::endl;
      exit(1);
    }

  //read inputs from parameter file
  std::map<std::string,double> floats;
  std::vector<std::vector<std::string> > files;
  std::vector<std::string> file_out[2];
  Input::FileReader reader(arg[1]);
  reader.initializeParameters(L,T,files,floats);
  double kappa=floats["kappa"];
  double mu=floats["mu"];
  double thetax=floats["thetax"];
  double thetay=floats["thetay"];
  double thetaz=floats["thetaz"];
  double thetat=floats["thetat"];
  bool debug=(bool)floats["debug"]; //0=no,!=0 =yes
  bool phys_base=(bool)floats["phys_base"]; //0=no,!=0 =yes

  //the weave!
  Base::Weave weave(L,T);

  //Prepare output file names
  for(size_t iud=0;iud<2;iud++)
    for(size_t j=0;j<files[1].size();j++)
	file_out[iud].push_back(files[1][j]+suff[iud]);
  
  //write some output concerning what is going to be read
  if(weave.isRoot())
    {
      std::cout<<"Lattice size: "<<L<<"x"<<L<<"x"<<L<<"x"<<T<<std::endl;
      std::cout<<"kappa="<<kappa<<", mu="<<mu<<std::endl;
      std::cout<<"thetax="<<thetax<<", ";
      std::cout<<"thetay="<<thetay<<", ";
      std::cout<<"thetaz="<<thetaz<<", ";
      std::cout<<"thetat="<<thetat<<std::endl;
      std::cout<<"will ";
      if(phys_base!=true) std::cout<<"not ";
      std::cout<<"write output in physical base"<<std::endl;
    }

  if(weave.isRoot() and debug==true)
    {
      std::cout<<"The following files are going to be read:"<<std::endl;
      
      for(size_t i=0;i<files.size();i++)
	for(size_t j=0;j<files[i].size();j++)
	  std::cout<<(files[i])[j]<<std::endl;

      for(size_t iud=0;iud<2;iud++)
	{
	  std::cout<<std::endl<<"Output files for quark "<<iud<<std::endl;
	  for(size_t j=0;j<files[1].size();j++)
	    std::cout<<file_out[iud][j]<<std::endl;
	}
    }

  //read the gauge configuration, needed to create u and d propagators
  Core::Field<QCD::Gauge> gauge_field(L,T);
  Tool::IO::load(&gauge_field,files[0][0],Tool::IO::fileILDG);
  if(weave.isRoot()) std::cout<<std::endl<<"gauge field successfully loaded"<<std::endl<<std::endl;

  //read the solution of (DD^-1)
  PropagatorType DD_prop(L,T);
  Tool::IO::load(&DD_prop,files[1],Tool::IO::fileSCIDAC);
  if(weave.isRoot()) std::cout<<"(D+D-)^-1 prop successfully loaded"<<std::endl<<std::endl;

  //prduce the u and d propagator from the D+D-
  //we have defined 0 as up, which requires the application of Q-
  PropagatorType *ud_prop[2];
  ud_prop[0]=new PropagatorType(L,T);
  ud_prop[1]=new PropagatorType(L,T);

  //The first is the solution of D- (that with -i*mu*g5), the second of D+.
  //This is obtained multiplying bi D+,D-
  //The function reconstruct doublet is the full Dirac Operatorm, which include also the gamma5
  DD_prop.reconstruct_doublet(*(ud_prop[0]),*(ud_prop[1]),gauge_field,kappa,mu,thetat,thetax,thetay,thetaz); 
  
  for(size_t iud=0;iud<2;iud++)
    { 
      if(phys_base)
	{
	  ud_prop[iud]->rotateToPhysicalBasis(iud); //0,1 = (1-+ig5)/sqrt(2) inside the routine
	  //note that this is propagator so it rotate the opposite way of the Dirac Operator

	  if(debug and weave.isRoot()) std::cout<<"quark "<<iud<<" rotated"<<std::endl;
	}
      Tool::IO::save(ud_prop[iud],file_out[iud],Tool::IO::fileSCIDAC);
      if(weave.isRoot()) std::cout<<"quark "<<iud<<" file successfully saved"<<std::endl;

      delete ud_prop[iud];
    }

  if(weave.isRoot()) std::cout<<std::endl<<"everything ok so exiting!"<<std::endl<<std::endl;

  return EXIT_SUCCESS;
}
